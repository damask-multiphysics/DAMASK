!--------------------------------------------------------------------------------------------------
!* $Id$
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
   IO_write_jobBinaryFile
      
 use math
 
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
   rotation_tol, &
   mySpectralSolver
   
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results
 
 use DAMASK_spectral_Utilities, only: &
   boundaryCondition, &
   solutionState, &
   debugGeneral
 
 use DAMASK_spectral_SolverBasic
 use DAMASK_spectral_SolverAL
 
 implicit none
 
 type loadCase
   real(pReal), dimension (3,3) :: rotation               = math_I3                                 ! rotation of BC
   type(boundaryCondition) ::      P, &                                                             ! stress BC
                                   deformation                                                      ! deformation BC (Fdot or L)
   real(pReal) ::                  time                   = 0.0_pReal, &                            ! length of increment
                                   temperature            = 300.0_pReal                             ! isothermal starting conditions
   integer(pInt) ::                incs                   = 0_pInt, &                               ! number of increments
                                   outputfrequency        = 1_pInt, &                               ! frequency of result writes
                                   restartfrequency       = 0_pInt, &                               ! frequency of restart writes
                                   logscale               = 0_pInt                                  ! linear/logaritmic time inc flag
   logical ::                      followFormerTrajectory = .true.                                  ! follow trajectory of former loadcase 
 end type loadCase

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 real(pReal), dimension(9) :: temp_valueVector                                                      !> temporarily from loadcase file when reading in tensors
 logical,     dimension(9) :: temp_maskVector                                                       !> temporarily from loadcase file when reading in tensors
 integer(pInt), parameter  :: maxNchunksLoadcase = (1_pInt + 9_pInt)*3_pInt +&                      ! deformation, rotation, and stress
                                                   (1_pInt + 1_pInt)*5_pInt +&                      ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                    1_pInt, &                                       ! dropguessing
                              maxNchunksGeom     = 7_pInt, &                                        ! 4 identifiers, 3 values
                              myUnit             = 234_pInt
 integer(pInt), dimension(1_pInt + maxNchunksLoadcase*2_pInt) :: positions                          ! this is longer than needed for geometry parsing
 
 integer(pInt) :: &
   N_l    = 0_pInt, &
   N_t    = 0_pInt, &
   N_n    = 0_pInt, &
   N_Fdot = 0_pInt                                                                                 ! number of Fourier points

 character(len=1024) :: &
   line


!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal), dimension(3,3), parameter :: ones = 1.0_pReal, zeroes = 0.0_pReal 
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc = 1.0_pReal, timeinc_old = 0.0_pReal   ! elapsed time, begin of interval, time interval, previous time interval
 real(pReal) :: guessmode             
 real(pReal),    dimension(3,3) :: temp33_Real
 integer(pInt) :: i, j, k, l, errorID
 integer(pInt) :: currentLoadcase = 0_pInt, inc, &
                  totalIncsCounter = 0_pInt,&
                  notConvergedCounter = 0_pInt, convergedCounter = 0_pInt
 character(len=6)  :: loadcase_string
 
 type(loadCase), allocatable, dimension(:) ::  loadCases
 type(solutionState) solres

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
 call CPFEM_initAll(temperature = 300.0_pReal, element = 1_pInt, IP= 1_pInt)

 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectral_Driver init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ''
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

100 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                                                     ! sanity check
       call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                               ! error message for incomplete loadcase
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
         if (IO_lc(IO_stringValue(line,positions,i)) == 'l'.or. &                          ! in case of given L, set flag to true
             IO_lc(IO_stringValue(line,positions,i)) == 'velocitygrad'.or.&
             IO_lc(IO_stringValue(line,positions,i)) == 'velgrad'.or.&
             IO_lc(IO_stringValue(line,positions,i)) == 'velocitygradient') then
           loadCases(currentLoadCase)%deformation%myType = 'l'
         else
           loadCases(currentLoadCase)%deformation%myType = 'fdot'
         endif
         forall (j = 1_pInt:9_pInt) temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         enddo
         loadCases(currentLoadCase)%deformation%maskLogical = transpose(reshape(temp_maskVector,[ 3,3]))
         loadCases(currentLoadCase)%deformation%maskFloat   = merge(ones,zeroes,&
                                                        loadCases(currentLoadCase)%deformation%maskLogical)
         loadCases(currentLoadCase)%deformation%values      = math_plain9to33(temp_valueVector)
       case('p','pk1','piolakirchhoff','stress')
         temp_valueVector = 0.0_pReal
         forall (j = 1_pInt:9_pInt) temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         enddo
         loadCases(currentLoadCase)%P%maskLogical = transpose(reshape(temp_maskVector,[ 3,3]))
         loadCases(currentLoadCase)%P%maskFloat   = merge(ones,zeroes,&
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
         loadCases(currentLoadCase)%followFormerTrajectory = .false.                                              ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of currentLoadCase given in euler angles
         l = 0_pInt                                                                                 ! assuming values given in radians
         k = 1_pInt                                                                                 ! assuming keyword indicating degree/radians
         select case (IO_lc(IO_stringValue(line,positions,i+1_pInt)))
           case('deg','degree')
             l = 1_pInt                                                                             ! for conversion from degree to radian           
           case('rad','radian') 
           case default               
             k = 0_pInt                                                                             ! immediately readingk in angles, assuming radians
         end select
         forall(j = 1_pInt:3_pInt)  temp33_Real(j,1) = &
                                        IO_floatValue(line,positions,i+k+j) * real(l,pReal) * inRad
         loadCases(currentLoadCase)%rotation = math_EulerToR(temp33_Real(:,1))
       case('rotation','rot')                                                                       ! assign values for the rotation of currentLoadCase matrix
         temp_valueVector = 0.0_pReal
         forall (j = 1_pInt:9_pInt) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         loadCases(currentLoadCase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
101 close(myUnit)

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 loadCases(1)%followFormerTrajectory = .false.                                                             ! cannot guess along trajectory for first inc of first currentLoadCase
 errorID = 0_pInt
 checkLoadcases: do currentLoadCase = 1_pInt, size(loadCases)
   write (loadcase_string, '(i6)' ) currentLoadCase

   write(6,'(2x,a,i6)') 'load case: ', currentLoadCase

   if (.not. loadCases(currentLoadCase)%followFormerTrajectory) write(6,'(2x,a)') 'drop guessing along trajectory'
   if (loadCases(currentLoadCase)%deformation%myType=='l') then
     do j = 1_pInt, 3_pInt
       if (any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .true.) .and. &
           any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .false.)) errorID = 832_pInt               ! each row should be either fully or not at all defined
     enddo
     write(6,'(2x,a)') 'velocity gradient:'
   else
     write(6,'(2x,a)') 'deformation gradient rate:'
   endif
   write (6,'(3(3(3x,f12.7,1x)/))',advance='no') merge(math_transpose33(loadCases(currentLoadCase)%deformation%values),&
                  reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(loadCases(currentLoadCase)%deformation%maskLogical))
   write (6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'stress / GPa:',&
        1e-9_pReal*merge(math_transpose33(loadCases(currentLoadCase)%P%values),&
                         reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(loadCases(currentLoadCase)%P%maskLogical))
   if (any(loadCases(currentLoadCase)%rotation /= math_I3)) &
     write (6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'rotation of loadframe:',&
                                                          math_transpose33(loadCases(currentLoadCase)%rotation)
   write(6,'(2x,a,f12.6)') 'temperature:', loadCases(currentLoadCase)%temperature
   write(6,'(2x,a,f12.6)') 'time:       ', loadCases(currentLoadCase)%time
   write(6,'(2x,a,i5)')    'increments: ', loadCases(currentLoadCase)%incs
   write(6,'(2x,a,i5)')    'output  frequency:  ', loadCases(currentLoadCase)%outputfrequency
   write(6,'(2x,a,i5)')    'restart frequency:  ', loadCases(currentLoadCase)%restartfrequency

   if (any(loadCases(currentLoadCase)%P%maskLogical .eqv. loadCases(currentLoadCase)%deformation%maskLogical)) errorID = 831_pInt          ! exclusive or masking only
   if (any(loadCases(currentLoadCase)%P%maskLogical .and. transpose(loadCases(currentLoadCase)%P%maskLogical) .and. &
     reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
                                               errorID = 838_pInt                                   ! no rotation is allowed by stress BC
   if (any(abs(math_mul33x33(loadCases(currentLoadCase)%rotation,math_transpose33(loadCases(currentLoadCase)%rotation))&
                                      -math_I3) > reshape(spread(rotation_tol,1,9),[ 3,3]))&
                    .or. abs(math_det33(loadCases(currentLoadCase)%rotation)) > 1.0_pReal + rotation_tol)&
                                               errorID = 846_pInt                                   ! given rotation matrix contains strain
   if (loadCases(currentLoadCase)%time < 0.0_pReal)          errorID = 834_pInt                                   ! negative time increment
   if (loadCases(currentLoadCase)%incs < 1_pInt)             errorID = 835_pInt                                   ! non-positive incs count
   if (loadCases(currentLoadCase)%outputfrequency < 1_pInt)  errorID = 836_pInt                                   ! non-positive result frequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)
 enddo checkLoadcases

 select case (myspectralsolver)
 
   case (DAMASK_spectral_SolverBasic_label)
     call basic_init()
     
   case (DAMASK_spectral_SolverAL_label)
     call AL_init()
     
 end select 
!--------------------------------------------------------------------------------------------------
! write header of output file
 if (appendToOutFile) then
   open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut',&
                                   form='UNFORMATTED', position='APPEND', status='OLD')
 else
   open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut',&
                                   form='UNFORMATTED',status='REPLACE')
   write(538) 'load',       trim(loadCaseFile)
   write(538) 'workingdir', trim(getSolverWorkingDirectoryName())
   write(538) 'geometry',   trim(geometryFile)
   write(538) 'resolution', res
   write(538) 'dimension',  geomdim
   write(538) 'materialpoint_sizeResults', materialpoint_sizeResults
   write(538) 'loadcases',  size(loadCases)
   write(538) 'frequencies', loadCases%outputfrequency                                      ! one entry per currentLoadCase
   write(538) 'times', loadCases%time                                                       ! one entry per currentLoadCase
   write(538) 'logscales',  loadCases%logscale         
   write(538) 'increments', loadCases%incs                                                  ! one entry per currentLoadCase
   write(538) 'startingIncrement', restartInc - 1_pInt                                              ! start with writing out the previous inc
   write(538) 'eoh'                                                                                 ! end of header
   write(538) materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:mesh_NcpElems)              ! initial (non-deformed or read-in) results
   if (debugGeneral) write(6,'(a)') 'Header of result file written out'
 endif


!--------------------------------------------------------------------------------------------------
! loopping over loadcases
 loadCaseLooping: do currentLoadCase = 1_pInt, size(loadCases)
   time0 = time                                                                                     ! currentLoadCase start time                
   if (loadCases(currentLoadCase)%followFormerTrajectory) then
     guessmode = 1.0_pReal
   else
     guessmode = 0.0_pReal                                                                          ! change of load case, homogeneous guess for the first inc
   endif

!--------------------------------------------------------------------------------------------------
! loop oper incs defined in input file for current currentLoadCase
   incLooping: do inc = 1_pInt, loadCases(currentLoadCase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt                                                 

!--------------------------------------------------------------------------------------------------
! forwarding time
     timeinc_old = timeinc
     if (loadCases(currentLoadCase)%logscale == 0_pInt) then                                                      ! linear scale
       timeinc = loadCases(currentLoadCase)%time/loadCases(currentLoadCase)%incs                                                ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
     else
       if (currentLoadCase == 1_pInt) then                                                                 ! 1st currentLoadCase of logarithmic scale            
         if (inc == 1_pInt) then                                                                    ! 1st inc of 1st currentLoadCase of logarithmic scale
           timeinc = loadCases(1)%time*(2.0_pReal**real(    1_pInt-loadCases(1)%incs ,pReal))                     ! assume 1st inc is equal to 2nd 
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
     time = time + timeinc

     if(totalIncsCounter >= restartInc)  then                                                       ! do calculations (otherwise just forwarding) 

!--------------------------------------------------------------------------------------------------
! report begin of new increment
       write(6,'(a)') '##################################################################'
       write(6,'(A,I5.5,A,es12.5)') 'Increment ', totalIncsCounter, ' Time ',time
       
       select case (myspectralsolver)
       
         case (DAMASK_spectral_SolverBasic_label)
           solres = basic_solution (&
               guessmode,timeinc,timeinc_old, &
                P_BC              = loadCases(currentLoadCase)%P, &
                F_BC              = loadCases(currentLoadCase)%deformation, &
                temperature_bc    = loadCases(currentLoadCase)%temperature, &
                rotation_BC       = loadCases(currentLoadCase)%rotation)
           
          case (DAMASK_spectral_SolverAL_label)
            solres = AL_solution (&
               guessmode,timeinc,timeinc_old, &
                P_BC              = loadCases(currentLoadCase)%P, &
                F_BC              = loadCases(currentLoadCase)%deformation, &
                temperature_bc    = loadCases(currentLoadCase)%temperature, &
                rotation_BC       = loadCases(currentLoadCase)%rotation)
           
       end select 
 
       write(6,'(a)') ''
       write(6,'(a)') '=================================================================='
       if(solres%converged) then
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' converged'
       else
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       endif

       if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0_pInt) then                                    ! at output frequency
         write(6,'(a)') ''
         write(6,'(a)') '... writing results to file ......................................'
         write(538)  materialpoint_results       ! write result to file
       endif
       
     endif ! end calculation/forwarding
     guessmode = 1.0_pReal                                                                          ! keep guessing along former trajectory during same currentLoadCase

    enddo incLooping
 enddo loadCaseLooping
 
 select case (myspectralsolver)
 
   case (DAMASK_spectral_SolverBasic_label)
     call basic_destroy()
     
   case (DAMASK_spectral_SolverAL_label)
     call AL_destroy()
     
 end select
 
 write(6,'(a)') ''
 write(6,'(a)') '##################################################################'
 write(6,'(i6.6,a,i6.6,a,f5.1,a)') convergedCounter, ' out of ', &
                                   notConvergedCounter + convergedCounter, ' (', &
                                   real(convergedCounter, pReal)/&
                                   real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, &
                                   ' %) increments converged!'
 close(538)
 if (notConvergedCounter > 0_pInt) call quit(3_pInt)
 call quit(0_pInt)

end program DAMASK_spectral_Driver




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine quit(stop_id)
 use prec, only: &
   pInt
   
 implicit none
 integer(pInt), intent(in) :: stop_id
 integer, dimension(8) :: dateAndTime                                                               ! type default integer

 call date_and_time(values = dateAndTime)
 write(6,'(/,a)') 'DAMASK terminated on:'
 write(6,'(a,2(i2.2,a),i4.4)') 'Date:               ',dateAndTime(3),'/',&
                                                      dateAndTime(2),'/',&
                                                      dateAndTime(1) 
 write(6,'(a,2(i2.2,a),i2.2)') 'Time:               ',dateAndTime(5),':',&
                                                      dateAndTime(6),':',&
                                                      dateAndTime(7)  
 if (stop_id == 0_pInt) stop 0                                                                      ! normal termination
 if (stop_id <  0_pInt) then                                                                        ! trigger regridding
   write(0,'(a,i6)') 'restart at ', stop_id*(-1_pInt)
   stop 2
 endif
 if (stop_id == 3_pInt) stop 3                                                                      ! not all steps converged
 stop 1                                                                                             ! error (message from IO_error)
end subroutine
