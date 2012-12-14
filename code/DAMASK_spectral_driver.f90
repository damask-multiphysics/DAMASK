!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the various spectral solvers
!--------------------------------------------------------------------------------------------------
program DAMASK_spectral_Driver
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use DAMASK_interface, only: &
   DAMASK_interface_init, &
   loadCaseFile, &
   geometryFile, &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   appendToOutFile
 use prec, only: &
   pInt, &
   pReal, &
   DAMASK_NaN
 use IO, only: &
   IO_isBlank, &
   IO_open_file, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_error, &
   IO_lc, &
   IO_read_jobBinaryFile, &
   IO_write_jobBinaryFile, &
   IO_intOut
 use math                                                                                           ! need to include the whole module for FFTW
 use mesh,  only : &
   res, &
   geomdim, &
   mesh_NcpElems
 use CPFEM, only: &
   CPFEM_initAll
 use FEsolving, only: &
   restartWrite, &
   restartInc
 use numerics, only: &
   maxCutBack, &
   rotation_tol, &
   mySpectralSolver, &
   regridMode
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results
 use DAMASK_spectral_Utilities, only: &
   tBoundaryCondition, &
   tSolutionState, &
   debugGeneral, &
   cutBack
 use DAMASK_spectral_SolverBasic
#ifdef PETSc
 use DAMASK_spectral_SolverBasicPETSC
 use DAMASK_spectral_SolverAL
#endif
 
 implicit none
 type tLoadCase
   real(pReal), dimension (3,3) :: rotation               = math_I3                                 !< rotation of BC
   type(tBoundaryCondition) ::     P, &                                                             !< stress BC
                                   deformation                                                      !< deformation BC (Fdot or L)
   real(pReal) ::                  time                   = 0.0_pReal, &                            !< length of increment
                                   temperature            = 300.0_pReal                             !< isothermal starting conditions
   integer(pInt) ::                incs                   = 0_pInt, &                               !< number of increments
                                   outputfrequency        = 1_pInt, &                               !< frequency of result writes
                                   restartfrequency       = 0_pInt, &                               !< frequency of restart writes
                                   logscale               = 0_pInt                                  !< linear/logaritmic time inc flag
   logical ::                      followFormerTrajectory = .true.                                  !< follow trajectory of former loadcase 
 end type tLoadCase

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 real(pReal), dimension(9) :: temp_valueVector = 0.0_pReal                                          !< temporarily from loadcase file when reading in tensors (initialize to 0.0)
 logical,     dimension(9) :: temp_maskVector  = .false.                                            !< temporarily from loadcase file when reading in tensors
 integer(pInt), parameter  :: maxNchunksLoadcase = (1_pInt + 9_pInt)*3_pInt +&                      ! deformation, rotation, and stress
                                                   (1_pInt + 1_pInt)*5_pInt +&                      ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                    1_pInt                                          ! dropguessing
 integer(pInt), dimension(1_pInt + maxNchunksLoadcase*2_pInt) :: positions                          ! this is longer than needed for geometry parsing
 
 integer(pInt) :: &
   N_l    = 0_pInt, &
   N_t    = 0_pInt, &
   N_n    = 0_pInt, &
   N_Fdot = 0_pInt, &                                                                               !
   myUnit = 234_pInt
 character(len=1024) :: &
   line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal), dimension(3,3), parameter :: ones = 1.0_pReal, zeros = 0.0_pReal 
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc = 1.0_pReal, timeinc_old = 0.0_pReal   ! elapsed time, begin of interval, time interval, previous time interval
 logical :: guess
 integer(pInt) :: i, j, k, l, errorID, cutBackLevel = 0_pInt, stepFraction = 0_pInt
 integer(pInt) :: currentLoadcase = 0_pInt, inc, &
                  totalIncsCounter = 0_pInt,&
                  notConvergedCounter = 0_pInt, convergedCounter = 0_pInt, &
                  resUnit = 0_pInt, statUnit = 0_pInt, lastRestartWritten = 0_pInt
 character(len=6)  :: loadcase_string
 character(len=1024)  :: incInfo
 type(tLoadCase), allocatable, dimension(:) ::  loadCases
 type(tSolutionState) solres

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
 call CPFEM_initAll(temperature = 300.0_pReal, element = 1_pInt, IP= 1_pInt)
 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectral_driver init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! reading basic information from load case file and allocate data structure containing load cases
 call IO_open_file(myUnit,trim(loadCaseFile))
 rewind(myUnit)
 do
   read(myUnit,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do i = 1_pInt, maxNchunksLoadcase, 1_pInt                                                        ! reading compulsory parameters for loadcase
       select case (IO_lc(IO_stringValue(line,positions,i)))
            case('l','velocitygrad','velgrad','velocitygradient')
                 N_l = N_l + 1_pInt
            case('fdot','dotf')
                 N_Fdot = N_Fdot + 1_pInt
            case('t','time','delta')
                 N_t = N_t + 1_pInt
            case('n','incs','increments','steps','logincs','logincrements','logsteps')
                 N_n = N_n + 1_pInt
        end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

100 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                                                  ! sanity check
       call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                                ! error message for incomplete loadcase
 allocate (loadCases(N_n))
 loadCases%P%myType='p'

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(myUnit)
 do
   read(myUnit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   currentLoadCase = currentLoadCase + 1_pInt
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do i = 1_pInt,maxNchunksLoadcase
     select case (IO_lc(IO_stringValue(line,positions,i)))
       case('fdot','dotf','l','velocitygrad','velgrad','velocitygradient')                          ! assign values for the deformation BC matrix
         temp_valueVector = 0.0_pReal
         if (IO_lc(IO_stringValue(line,positions,i)) == 'fdot'.or. &                                ! in case of Fdot, set type to fdot
             IO_lc(IO_stringValue(line,positions,i)) == 'dotf') then
           loadCases(currentLoadCase)%deformation%myType = 'fdot'
         else
           loadCases(currentLoadCase)%deformation%myType = 'l'
         endif
         forall (j = 1_pInt:9_pInt) temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         enddo
         loadCases(currentLoadCase)%deformation%maskLogical = transpose(reshape(temp_maskVector,[ 3,3]))
         loadCases(currentLoadCase)%deformation%maskFloat   = merge(ones,zeros,&
                                                        loadCases(currentLoadCase)%deformation%maskLogical)
         loadCases(currentLoadCase)%deformation%values      = math_plain9to33(temp_valueVector)
       case('p','pk1','piolakirchhoff','stress')
         temp_valueVector = 0.0_pReal
         forall (j = 1_pInt:9_pInt) temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         enddo
         loadCases(currentLoadCase)%P%maskLogical = transpose(reshape(temp_maskVector,[ 3,3]))
         loadCases(currentLoadCase)%P%maskFloat   = merge(ones,zeros,&
                                                        loadCases(currentLoadCase)%P%maskLogical)
         loadCases(currentLoadCase)%P%values      = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                                     ! increment time
         loadCases(currentLoadCase)%time = IO_floatValue(line,positions,i+1_pInt)
       case('temp','temperature')                                                                   ! starting temperature
         loadCases(currentLoadCase)%temperature = IO_floatValue(line,positions,i+1_pInt)
       case('n','incs','increments','steps')                                                        ! number of increments
         loadCases(currentLoadCase)%incs = IO_intValue(line,positions,i+1_pInt)
       case('logincs','logincrements','logsteps')                                                   ! number of increments (switch to log time scaling)
         loadCases(currentLoadCase)%incs = IO_intValue(line,positions,i+1_pInt)
         loadCases(currentLoadCase)%logscale = 1_pInt
       case('f','freq','frequency','outputfreq')                                                    ! frequency of result writings
         loadCases(currentLoadCase)%outputfrequency = IO_intValue(line,positions,i+1_pInt)                
       case('r','restart','restartwrite')                                                           ! frequency of writing restart information
         loadCases(currentLoadCase)%restartfrequency = max(0_pInt,IO_intValue(line,positions,i+1_pInt))                
       case('guessreset','dropguessing')
         loadCases(currentLoadCase)%followFormerTrajectory = .false.                                ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of currentLoadCase given in euler angles
         temp_valueVector = 0.0_pReal
         l = 1_pInt                                                                                 ! assuming values given in degrees
         k = 1_pInt                                                                                 ! assuming keyword indicating degree/radians present
         select case (IO_lc(IO_stringValue(line,positions,i+1_pInt)))
           case('deg','degree')
           case('rad','radian')                                                                     ! don't convert from degree to radian           
             l = 0_pInt
           case default   
             k = 0_pInt           
         end select
         forall(j = 1_pInt:3_pInt)  temp_valueVector(j) = IO_floatValue(line,positions,i+k+j)
         if (l == 1_pInt) temp_valueVector(1:3) = temp_valueVector(1:3) * inRad                     ! convert to rad
         loadCases(currentLoadCase)%rotation = math_EulerToR(temp_valueVector(1:3))                 ! convert rad Eulers to rotation matrix
       case('rotation','rot')                                                                       ! assign values for the rotation of currentLoadCase matrix
         temp_valueVector = 0.0_pReal
         forall (j = 1_pInt:9_pInt) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         loadCases(currentLoadCase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
101 close(myUnit)

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 loadCases(1)%followFormerTrajectory = .false.                                                      ! cannot guess along trajectory for first inc of first currentLoadCase
 errorID = 0_pInt
 checkLoadcases: do currentLoadCase = 1_pInt, size(loadCases)
   write (loadcase_string, '(i6)' ) currentLoadCase
   write(6,'(1x,a,i6)') 'load case: ', currentLoadCase

   if (.not. loadCases(currentLoadCase)%followFormerTrajectory) &
     write(6,'(2x,a)') 'drop guessing along trajectory'
   if (loadCases(currentLoadCase)%deformation%myType=='l') then
     do j = 1_pInt, 3_pInt
       if (any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .true.) .and. &
           any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .false.)) &
                                                                  errorID = 832_pInt                ! each row should be either fully or not at all defined
     enddo
     write(6,'(2x,a)') 'velocity gradient:'
   else
     write(6,'(2x,a)') 'deformation gradient rate:'
   endif
   write(6,'(3(3(3x,f12.7,1x)/))',advance='no') &
              merge(math_transpose33(loadCases(currentLoadCase)%deformation%values), &
              reshape(spread(huge(1.0_pReal),1,9),[ 3,3]), &
              transpose(loadCases(currentLoadCase)%deformation%maskLogical))
   if (any(loadCases(currentLoadCase)%P%maskLogical .eqv. &
            loadCases(currentLoadCase)%deformation%maskLogical)) errorID = 831_pInt                 ! exclusive or masking only
   if (any(loadCases(currentLoadCase)%P%maskLogical .and. &                                   
           transpose(loadCases(currentLoadCase)%P%maskLogical) .and. &
           reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
           errorID = 838_pInt                                                                       ! no rotation is allowed by stress BC
   write(6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'stress / GPa:',&
              1e-9_pReal*merge(math_transpose33(loadCases(currentLoadCase)%P%values),&
              reshape(spread(huge(1.0_pReal),1,9),[ 3,3]),&
              transpose(loadCases(currentLoadCase)%P%maskLogical))
  if (any(abs(math_mul33x33(loadCases(currentLoadCase)%rotation, &
              math_transpose33(loadCases(currentLoadCase)%rotation))-math_I3) >&
              reshape(spread(rotation_tol,1,9),[ 3,3]))&
              .or. abs(math_det33(loadCases(currentLoadCase)%rotation)) > &
              1.0_pReal + rotation_tol) errorID = 846_pInt                                          ! given rotation matrix contains strain
   if (any(loadCases(currentLoadCase)%rotation /= math_I3)) &
     write(6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'rotation of loadframe:',&
              math_transpose33(loadCases(currentLoadCase)%rotation)
   write(6,'(2x,a,f12.6)') 'temperature:', loadCases(currentLoadCase)%temperature
   if (loadCases(currentLoadCase)%time < 0.0_pReal)          errorID = 834_pInt                     ! negative time increment
   write(6,'(2x,a,f12.6)') 'time:       ', loadCases(currentLoadCase)%time
   if (loadCases(currentLoadCase)%incs < 1_pInt)             errorID = 835_pInt                     ! non-positive incs count
   write(6,'(2x,a,i5)')    'increments: ', loadCases(currentLoadCase)%incs
   if (loadCases(currentLoadCase)%outputfrequency < 1_pInt)  errorID = 836_pInt                     ! non-positive result frequency
   write(6,'(2x,a,i5)')    'output  frequency:  ', &
              loadCases(currentLoadCase)%outputfrequency
   write(6,'(2x,a,i5,/)')    'restart frequency:  ', &
              loadCases(currentLoadCase)%restartfrequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)
 enddo checkLoadcases

 select case (myspectralsolver)
   case (DAMASK_spectral_SolverBasic_label)
     call basic_init()
#ifdef PETSc
   case (DAMASK_spectral_SolverBasicPETSc_label)
     call basicPETSc_init()
   case (DAMASK_spectral_SolverAL_label)
     call AL_init()
#endif
   case default
      call IO_error(error_ID = 891, ext_msg = trim(myspectralsolver))
 end select 
 
!--------------------------------------------------------------------------------------------------
! write header of output file
 if (appendToOutFile) then
   open(newunit=resUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.spectralOut',form='UNFORMATTED', position='APPEND', status='OLD')
   open(newunit=statUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.sta',form='FORMATTED', position='APPEND', status='OLD')
 else
   open(newunit=resUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.spectralOut',form='UNFORMATTED',status='REPLACE')
   write(resUnit) 'load',       trim(loadCaseFile)
   write(resUnit) 'workingdir', trim(getSolverWorkingDirectoryName())
   write(resUnit) 'geometry',   trim(geometryFile)
   write(resUnit) 'resolution', res
   write(resUnit) 'dimension',  geomdim
   write(resUnit) 'materialpoint_sizeResults', materialpoint_sizeResults
   write(resUnit) 'loadcases',  size(loadCases)
   write(resUnit) 'frequencies', loadCases%outputfrequency                                          ! one entry per currentLoadCase
   write(resUnit) 'times', loadCases%time                                                           ! one entry per currentLoadCase
   write(resUnit) 'logscales',  loadCases%logscale         
   write(resUnit) 'increments', loadCases%incs                                                      ! one entry per currentLoadCase
   write(resUnit) 'startingIncrement', restartInc - 1_pInt                                          ! start with writing out the previous inc
   write(resUnit) 'eoh'                                                                             ! end of header
   write(resUnit) materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:mesh_NcpElems)    ! initial (non-deformed or read-in) results
   open(newunit=statUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.sta',form='FORMATTED',status='REPLACE')
   write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'
   if (debugGeneral) write(6,'(a)') 'Header of result file written out'
 endif
!--------------------------------------------------------------------------------------------------
! loopping over loadcases
 loadCaseLooping: do currentLoadCase = 1_pInt, size(loadCases)
   time0 = time                                                                                     ! currentLoadCase start time                
   if (loadCases(currentLoadCase)%followFormerTrajectory) then
     guess = .true.
   else
     guess = .false.                                                                                ! change of load case, homogeneous guess for the first inc
   endif

!--------------------------------------------------------------------------------------------------
! loop oper incs defined in input file for current currentLoadCase
   incLooping: do inc = 1_pInt, loadCases(currentLoadCase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt                                                 

!--------------------------------------------------------------------------------------------------
! forwarding time
     timeinc_old = timeinc
     if (loadCases(currentLoadCase)%logscale == 0_pInt) then                                        ! linear scale
       timeinc = loadCases(currentLoadCase)%time/loadCases(currentLoadCase)%incs                    ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
     else
       if (currentLoadCase == 1_pInt) then                                                          ! 1st currentLoadCase of logarithmic scale            
         if (inc == 1_pInt) then                                                                    ! 1st inc of 1st currentLoadCase of logarithmic scale
           timeinc = loadCases(1)%time*(2.0_pReal**real(    1_pInt-loadCases(1)%incs ,pReal))       ! assume 1st inc is equal to 2nd 
         else                                                                                       ! not-1st inc of 1st currentLoadCase of logarithmic scale
           timeinc = loadCases(1)%time*(2.0_pReal**real(inc-1_pInt-loadCases(1)%incs ,pReal))
         endif
       else                                                                                         ! not-1st currentLoadCase of logarithmic scale
         timeinc = time0 *( (1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real(          inc,pReal)/&
                                                                real(loadCases(currentLoadCase)%incs ,pReal))&
                           -(1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                                 real(loadCases(currentLoadCase)%incs ,pReal)) )
       endif
     endif
     timeinc = timeinc / 2.0_pReal**real(cutBackLevel,pReal)

     if(totalIncsCounter >= restartInc)  then                                                       ! do calculations (otherwise just forwarding) 
       stepFraction = 0_pInt
!--------------------------------------------------------------------------------------------------
! report begin of new increment
       subIncLooping: do while (stepFraction/2_pInt**cutBackLevel <1_pInt)
         time = time + timeinc 
         stepFraction = stepFraction + 1_pInt 
         write(6,'(1/,a)') '###########################################################################'
         write(6,'(a,es12.5'//&
                 ',a,'//IO_intOut(inc)//',a,'//IO_intOut(loadCases(currentLoadCase)%incs)//&
                 ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(2_pInt**cutBackLevel)//&
                 ',a,'//IO_intOut(currentLoadCase)//',a,'//IO_intOut(size(loadCases))//')') &
                 'Time', time, &
                 's: Increment ', inc, '/', loadCases(currentLoadCase)%incs,&
                 '-', stepFraction, '/', 2_pInt**cutBackLevel,&
                 ' of load case ', currentLoadCase,'/',size(loadCases)
         write(incInfo,'(a,'//IO_intOut(totalIncsCounter)//',a,'//IO_intOut(sum(loadCases(:)%incs))//&
               ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(2_pInt**cutBackLevel)//')') &
               'Inc. ',totalIncsCounter,'/',sum(loadCases(:)%incs),&
               '-',stepFraction, '/', 2_pInt**cutBackLevel
         select case(myspectralsolver)
         
           case (DAMASK_spectral_SolverBasic_label)
             solres = basic_solution (&
                  incInfo, guess,timeinc,timeinc_old, &
                  P_BC              = loadCases(currentLoadCase)%P, &
                  F_BC              = loadCases(currentLoadCase)%deformation, &
                  temperature_bc    = loadCases(currentLoadCase)%temperature, &
                  rotation_BC       = loadCases(currentLoadCase)%rotation)
#ifdef PETSc
           case (DAMASK_spectral_SolverBasicPETSC_label)
             solres = BasicPETSC_solution (&
                 incInfo, guess,timeinc,timeinc_old, &
                 P_BC              = loadCases(currentLoadCase)%P, &
                 F_BC              = loadCases(currentLoadCase)%deformation, &
                 temperature_bc    = loadCases(currentLoadCase)%temperature, &
                 rotation_BC       = loadCases(currentLoadCase)%rotation)
            
           case (DAMASK_spectral_SolverAL_label)
             solres = AL_solution (&
                 incInfo, guess,timeinc,timeinc_old, &
                 P_BC              = loadCases(currentLoadCase)%P, &
                 F_BC              = loadCases(currentLoadCase)%deformation, &
                 temperature_bc    = loadCases(currentLoadCase)%temperature, &
                 rotation_BC       = loadCases(currentLoadCase)%rotation)
#endif
         end select 
         cutBack = .False.
         if(solres%termIll .or. .not. solres%converged) then                                        ! no solution found
           if (cutBackLevel < maxCutBack) then                                                      ! do cut back
             write(6,'(/,a)') 'cut back detected'
             cutBack = .True.
             stepFraction = (stepFraction - 1_pInt) * 2_pInt**cutBackLevel
             cutBackLevel = cutBackLevel + 1_pInt
             time    = time - timeinc
             timeinc_old = timeinc
             timeinc = timeinc/2.0_pReal
           elseif (solres%termIll) then                                                             ! material point model cannot find a solution
             if(regridMode > 0_pInt) call quit(-1*(lastRestartWritten+1))
             call IO_error(850_pInt)
           else
             if(regridMode == 2_pInt) call quit(-1*(lastRestartWritten+1))
             guess = .true.                                                                         ! start guessing after first accepted (not converged) (sub)inc
           endif
         else
           guess = .true.                                                                           ! start guessing after first converged (sub)inc
         endif
       if(guess) &
         write(statUnit,*) inc, time, cutBackLevel, solres%converged, solres%iterationsNeeded
       enddo subIncLooping
       cutBackLevel = max(0_pInt, cutBackLevel - 1_pInt)                                            ! try half subincs next inc
       if(solres%converged) then
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(A,'//IO_intOut(totalIncsCounter)//',A)') &
                                     'increment ', totalIncsCounter, ' converged'
       else
         write(6,'(A,'//IO_intOut(totalIncsCounter)//',A)') &
                                     'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       endif

       if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0_pInt) then                      ! at output frequency
         write(6,'(1/,a)') '... writing results to file ......................................'
         write(resUnit)  materialpoint_results                                                          ! write result to file
       endif
       if( loadCases(currentLoadCase)%restartFrequency > 0_pInt .and. &
                      mod(inc,loadCases(currentLoadCase)%restartFrequency) == 0_pInt) then                        ! at frequency of writing restart information set restart parameter for FEsolving (first call to CPFEM_general will write ToDo: true?) 
         restartWrite = .true.
         lastRestartWritten = inc
       endif 
     else !just time forwarding
       time = time + timeinc
       guess = .true.
     endif                                                                                          ! end calculation/forwarding

    enddo incLooping
 enddo loadCaseLooping
 
 select case (myspectralsolver)
 
   case (DAMASK_spectral_SolverBasic_label)
     call basic_destroy()
#ifdef PETSc
   case (DAMASK_spectral_SolverBasicPETSC_label)
     call BasicPETSC_destroy()
     
   case (DAMASK_spectral_SolverAL_label)
     call AL_destroy()
#endif 
 end select
 
 write(6,'(a)') ''
 write(6,'(a)') '##################################################################'
 write(6,'(i6.6,a,i6.6,a,f5.1,a)') convergedCounter, ' out of ', &
                                   notConvergedCounter + convergedCounter, ' (', &
                                   real(convergedCounter, pReal)/&
                                   real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, &
                                   ' %) increments converged!'
 close(resUnit)
 if (notConvergedCounter > 0_pInt) call quit(3_pInt)
 call quit(0_pInt)

end program DAMASK_spectral_Driver

#include "spectral_quit.f90"
