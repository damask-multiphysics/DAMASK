!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the various spectral solvers
!> @details doing cutbacking, forwarding in case of restart, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_spectral
 use, intrinsic :: &
   iso_fortran_env                                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use prec, only: &
   pInt, &
   pLongInt, &
   pReal, &
   tol_math_check, &
   dNeq
 use DAMASK_interface, only: &
   DAMASK_interface_init, &
   loadCaseFile, &
   geometryFile, &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   appendToOutFile
 use IO, only: &
   IO_read, &
   IO_isBlank, &
   IO_open_file, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_error, &
   IO_lc, &
   IO_intOut, &
   IO_warning, &
   IO_timeStamp, &
   IO_EOF
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_levelBasic
 use math                                                                                           ! need to include the whole module for FFTW
 use mesh, only: &
   grid, &
   geomSize
 use CPFEM2, only: &
   CPFEM_initAll
 use FEsolving, only: &
   restartWrite, &
   restartInc
 use numerics, only: &
   worldrank, &
   worldsize, &
   stagItMax, &
   maxCutBack, &
   spectral_solver, &
   continueCalculation
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results, &
   materialpoint_postResults
 use material, only: &
   thermal_type, &
   damage_type, &
   THERMAL_conduction_ID, &
   DAMAGE_nonlocal_ID
 use spectral_utilities, only: &
   utilities_init, &
   utilities_destroy, &
   tSolutionState, &
   tLoadCase, &
   cutBack, &
   nActiveFields, &
   FIELD_UNDEFINED_ID, &
   FIELD_MECH_ID, &
   FIELD_THERMAL_ID, &
   FIELD_DAMAGE_ID
 use spectral_mech_Basic
 use spectral_mech_AL
 use spectral_mech_Polarisation
 use spectral_damage
 use spectral_thermal
 

 implicit none

#include <petsc/finclude/petscsys.h>

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 real(pReal), dimension(9) :: temp_valueVector = 0.0_pReal                                          !< temporarily from loadcase file when reading in tensors (initialize to 0.0)
 logical,     dimension(9) :: temp_maskVector  = .false.                                            !< temporarily from loadcase file when reading in tensors
 integer(pInt), parameter  :: FILEUNIT         = 234_pInt                                           !< file unit, DAMASK IO does not support newunit feature
 integer(pInt), allocatable, dimension(:) :: chunkPos
 
 integer(pInt) :: &
   N_t   = 0_pInt, &                                                                                !< # of time indicators found in load case file 
   N_n   = 0_pInt, &                                                                                !< # of increment specifiers found in load case file
   N_def = 0_pInt                                                                                   !< # of rate of deformation specifiers found in load case file
 character(len=65536) :: &
   line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal), dimension(3,3), parameter :: &
   ones  = 1.0_pReal, &
   zeros = 0.0_pReal 
 integer(pInt), parameter :: &
   subStepFactor = 2_pInt                                                                           !< for each substep, divide the last time increment by 2.0
 real(pReal) :: &
   time = 0.0_pReal, &                                                                              !< elapsed time
   time0 = 0.0_pReal, &                                                                             !< begin of interval
   timeinc = 1.0_pReal, &                                                                           !< current time interval
   timeIncOld = 0.0_pReal, &                                                                        !< previous time interval
   remainingLoadCaseTime = 0.0_pReal                                                                !< remaining time of current load case
 logical :: &
   guess, &                                                                                         !< guess along former trajectory
   stagIterate
 integer(pInt) :: &
   i, j, k, l, field, &
   errorID, &
   cutBackLevel = 0_pInt, &                                                                         !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
   stepFraction = 0_pInt                                                                            !< fraction of current time interval
 integer(pInt) :: &
   currentLoadcase = 0_pInt, &                                                                      !< current load case
   inc, &                                                                                           !< current increment in current load case
   totalIncsCounter = 0_pInt, &                                                                     !< total # of increments
   convergedCounter = 0_pInt, &                                                                     !< # of converged increments
   notConvergedCounter = 0_pInt, &                                                                  !< # of non-converged increments
   resUnit = 0_pInt, &                                                                              !< file unit for results writing
   statUnit = 0_pInt, &                                                                             !< file unit for statistics output
   lastRestartWritten = 0_pInt, &                                                                   !< total increment # at which last restart information was written
   stagIter
 character(len=6)  :: loadcase_string
 character(len=1024)  :: incInfo                                                                    !< string parsed to solution with information about current load case
 type(tLoadCase), allocatable, dimension(:) :: loadCases                                            !< array of all load cases
 type(tSolutionState), allocatable, dimension(:) :: solres
 integer(MPI_OFFSET_KIND) :: fileOffset
 integer(MPI_OFFSET_KIND), dimension(:), allocatable :: outputSize
 integer(pInt), parameter :: maxByteOut = 2147483647-4096                                           !< limit of one file output write https://trac.mpich.org/projects/mpich/ticket/1742
 integer(pLongInt), dimension(2) :: outputIndex
 PetscErrorCode :: ierr
 external :: &
   quit, &
   MPI_file_open, &
   MPI_file_close, &
   MPI_file_seek, &
   MPI_file_get_position, &
   MPI_file_write, &
   MPI_abort, &
   MPI_allreduce, &
   PETScFinalize

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
 call CPFEM_initAll(el = 1_pInt, ip = 1_pInt)
 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)')   ' <<<+-  DAMASK_spectral init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
!--------------------------------------------------------------------------------------------------
! initialize field solver information
 nActiveFields = 1
 if (any(thermal_type  == THERMAL_conduction_ID  )) nActiveFields = nActiveFields + 1
 if (any(damage_type   == DAMAGE_nonlocal_ID     )) nActiveFields = nActiveFields + 1
 allocate(solres(nActiveFields))

!--------------------------------------------------------------------------------------------------
! reading basic information from load case file and allocate data structure containing load cases
 call IO_open_file(FILEUNIT,trim(loadCaseFile))
 rewind(FILEUNIT)
 do
   line = IO_read(FILEUNIT)
   if (trim(line) == IO_EOF) exit
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   chunkPos = IO_stringPos(line)
   do i = 1_pInt, chunkPos(1)                                                                       ! reading compulsory parameters for loadcase
     select case (IO_lc(IO_stringValue(line,chunkPos,i)))
       case('l','velocitygrad','velgrad','velocitygradient','fdot','dotf','f')
         N_def = N_def + 1_pInt
       case('t','time','delta')
         N_t = N_t + 1_pInt
       case('n','incs','increments','steps','logincs','logincrements','logsteps')
         N_n = N_n + 1_pInt
     end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

 if ((N_def /= N_n) .or. (N_n /= N_t) .or. N_n < 1_pInt) &                                          ! sanity check
   call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                                    ! error message for incomplete loadcase
 allocate (loadCases(N_n))                                                                          ! array of load cases
 loadCases%P%myType='p'
 
  do i = 1, size(loadCases)
   allocate(loadCases(i)%ID(nActiveFields))
   field = 1
   loadCases(i)%ID(field) = FIELD_MECH_ID           ! mechanical active by default
   if (any(thermal_type  == THERMAL_conduction_ID)) then ! thermal field active
     field = field + 1
     loadCases(i)%ID(field) = FIELD_THERMAL_ID 
   endif  
   if (any(damage_type   == DAMAGE_nonlocal_ID))  then ! damage field active
     field = field + 1
     loadCases(i)%ID(field) = FIELD_DAMAGE_ID
   endif
 enddo

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(FILEUNIT)
 do
   line = IO_read(FILEUNIT)
   if (trim(line) == IO_EOF) exit
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   currentLoadCase = currentLoadCase + 1_pInt
   chunkPos = IO_stringPos(line)
   do i = 1_pInt, chunkPos(1)
     select case (IO_lc(IO_stringValue(line,chunkPos,i)))
       case('fdot','dotf','l','velocitygrad','velgrad','velocitygradient','f')                      ! assign values for the deformation BC matrix
         temp_valueVector = 0.0_pReal
         if (IO_lc(IO_stringValue(line,chunkPos,i)) == 'fdot'.or. &                                 ! in case of Fdot, set type to fdot
             IO_lc(IO_stringValue(line,chunkPos,i)) == 'dotf') then
           loadCases(currentLoadCase)%deformation%myType = 'fdot'
         else if (IO_lc(IO_stringValue(line,chunkPos,i)) == 'f') then
           loadCases(currentLoadCase)%deformation%myType = 'f'
         else
           loadCases(currentLoadCase)%deformation%myType = 'l'
         endif
         do j = 1_pInt, 9_pInt
           temp_maskVector(j) = IO_stringValue(line,chunkPos,i+j) /= '*'                            ! true if not a *
         enddo
         do j = 1_pInt,9_pInt 
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,chunkPos,i+j)           ! read value where applicable
         enddo
         loadCases(currentLoadCase)%deformation%maskLogical = &                                     ! logical mask in 3x3 notation
               transpose(reshape(temp_maskVector,[ 3,3]))  
         loadCases(currentLoadCase)%deformation%maskFloat   = &                                     ! float (1.0/0.0) mask in 3x3 notation
               merge(ones,zeros,loadCases(currentLoadCase)%deformation%maskLogical)
         loadCases(currentLoadCase)%deformation%values = math_plain9to33(temp_valueVector)          ! values in 3x3 notation
       case('p','pk1','piolakirchhoff','stress', 's')
         temp_valueVector = 0.0_pReal
         do j = 1_pInt, 9_pInt
           temp_maskVector(j) = IO_stringValue(line,chunkPos,i+j) /= '*'                            ! true if not an asterisk
         enddo
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,chunkPos,i+j)           ! read value where applicable
         enddo
         loadCases(currentLoadCase)%P%maskLogical = transpose(reshape(temp_maskVector,[ 3,3]))
         loadCases(currentLoadCase)%P%maskFloat   = merge(ones,zeros,&
                                                        loadCases(currentLoadCase)%P%maskLogical)
         loadCases(currentLoadCase)%P%values      = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                                     ! increment time
         loadCases(currentLoadCase)%time = IO_floatValue(line,chunkPos,i+1_pInt)
       case('n','incs','increments','steps')                                                        ! number of increments
         loadCases(currentLoadCase)%incs = IO_intValue(line,chunkPos,i+1_pInt)
       case('logincs','logincrements','logsteps')                                                   ! number of increments (switch to log time scaling)
         loadCases(currentLoadCase)%incs = IO_intValue(line,chunkPos,i+1_pInt)
         loadCases(currentLoadCase)%logscale = 1_pInt
       case('freq','frequency','outputfreq')                                                        ! frequency of result writings
         loadCases(currentLoadCase)%outputfrequency = IO_intValue(line,chunkPos,i+1_pInt)                
       case('r','restart','restartwrite')                                                           ! frequency of writing restart information
         loadCases(currentLoadCase)%restartfrequency = &
               max(0_pInt,IO_intValue(line,chunkPos,i+1_pInt))                
       case('guessreset','dropguessing')
         loadCases(currentLoadCase)%followFormerTrajectory = .false.                                ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of currentLoadCase given in euler angles
         temp_valueVector = 0.0_pReal
         l = 1_pInt                                                                                 ! assuming values given in degrees
         k = 1_pInt                                                                                 ! assuming keyword indicating degree/radians present
         select case (IO_lc(IO_stringValue(line,chunkPos,i+1_pInt)))
           case('deg','degree')
           case('rad','radian')                                                                     ! don't convert from degree to radian           
             l = 0_pInt
           case default   
             k = 0_pInt           
         end select
         do j = 1_pInt, 3_pInt
           temp_valueVector(j) = IO_floatValue(line,chunkPos,i+k+j)
         enddo
         if (l == 1_pInt) temp_valueVector(1:3) = temp_valueVector(1:3) * inRad                     ! convert to rad
         loadCases(currentLoadCase)%rotation = math_EulerToR(temp_valueVector(1:3))                 ! convert rad Eulers to rotation matrix
       case('rotation','rot')                                                                       ! assign values for the rotation of currentLoadCase matrix
         temp_valueVector = 0.0_pReal
         do j = 1_pInt, 9_pInt
           temp_valueVector(j) = IO_floatValue(line,chunkPos,i+j)
         enddo
         loadCases(currentLoadCase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
 close(FILEUNIT) 

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 loadCases(1)%followFormerTrajectory = .false.                                                      ! cannot guess along trajectory for first inc of first currentLoadCase
 errorID = 0_pInt
 if (worldrank == 0) then
   checkLoadcases: do currentLoadCase = 1_pInt, size(loadCases)
     write (loadcase_string, '(i6)' ) currentLoadCase
     write(6,'(1x,a,i6)') 'load case: ', currentLoadCase
     if (.not. loadCases(currentLoadCase)%followFormerTrajectory) &
       write(6,'(2x,a)') 'drop guessing along trajectory'
     if (loadCases(currentLoadCase)%deformation%myType=='l') then
       do j = 1_pInt, 3_pInt
         if (any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .true.) .and. &
             any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .false.)) &
                                                                    errorID = 832_pInt              ! each row should be either fully or not at all defined
       enddo
       write(6,'(2x,a)') 'velocity gradient:'
     else if (loadCases(currentLoadCase)%deformation%myType=='f') then
       write(6,'(2x,a)') 'deformation gradient at end of load case:'
     else
       write(6,'(2x,a)') 'deformation gradient rate:'
     endif
     do i = 1_pInt, 3_pInt; do j = 1_pInt, 3_pInt
       if(loadCases(currentLoadCase)%deformation%maskLogical(i,j)) then
         write(6,'(2x,f12.7)',advance='no') loadCases(currentLoadCase)%deformation%values(i,j)
       else
         write(6,'(2x,12a)',advance='no') '    *       '
         endif
       enddo; write(6,'(/)',advance='no')
     enddo
     if (any(loadCases(currentLoadCase)%P%maskLogical .eqv. &
              loadCases(currentLoadCase)%deformation%maskLogical)) errorID = 831_pInt               ! exclusive or masking only
     if (any(loadCases(currentLoadCase)%P%maskLogical .and. &                                   
             transpose(loadCases(currentLoadCase)%P%maskLogical) .and. &
             reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
             errorID = 838_pInt                                                                     ! no rotation is allowed by stress BC
     write(6,'(2x,a)') 'stress / GPa:'
     do i = 1_pInt, 3_pInt; do j = 1_pInt, 3_pInt
       if(loadCases(currentLoadCase)%P%maskLogical(i,j)) then
         write(6,'(2x,f12.7)',advance='no') loadCases(currentLoadCase)%P%values(i,j)*1e-9_pReal
       else
         write(6,'(2x,12a)',advance='no') '    *       '
       endif
       enddo; write(6,'(/)',advance='no')
     enddo
    if (any(abs(math_mul33x33(loadCases(currentLoadCase)%rotation, &
                math_transpose33(loadCases(currentLoadCase)%rotation))-math_I3) >&
                reshape(spread(tol_math_check,1,9),[ 3,3]))&
                .or. abs(math_det33(loadCases(currentLoadCase)%rotation)) > &
                1.0_pReal + tol_math_check) errorID = 846_pInt                                      ! given rotation matrix contains strain
     if (any(dNeq(loadCases(currentLoadCase)%rotation, math_I3))) &
       write(6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'rotation of loadframe:',&
                math_transpose33(loadCases(currentLoadCase)%rotation)
     if (loadCases(currentLoadCase)%time < 0.0_pReal)          errorID = 834_pInt                   ! negative time increment
     write(6,'(2x,a,f12.6)') 'time:       ', loadCases(currentLoadCase)%time
     if (loadCases(currentLoadCase)%incs < 1_pInt)             errorID = 835_pInt                   ! non-positive incs count
     write(6,'(2x,a,i5)')    'increments: ', loadCases(currentLoadCase)%incs
     if (loadCases(currentLoadCase)%outputfrequency < 1_pInt)  errorID = 836_pInt                   ! non-positive result frequency
     write(6,'(2x,a,i5)')    'output  frequency:  ', &
                loadCases(currentLoadCase)%outputfrequency
     write(6,'(2x,a,i5,/)')    'restart frequency:  ', &
                loadCases(currentLoadCase)%restartfrequency
     if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)             ! exit with error message
   enddo checkLoadcases
 endif

!--------------------------------------------------------------------------------------------------
! doing initialization depending on selected solver 
 call Utilities_init()
 do field = 1, nActiveFields
   select case (loadCases(1)%ID(field))
     case(FIELD_MECH_ID)
       select case (spectral_solver)
         case (DAMASK_spectral_SolverBasicPETSc_label)
           call basicPETSc_init
         case (DAMASK_spectral_SolverAL_label)
           if(iand(debug_level(debug_spectral),debug_levelBasic)/= 0 .and. worldrank == 0_pInt) &
           call IO_warning(42_pInt, ext_msg='debug Divergence')
           call AL_init
         
         case (DAMASK_spectral_SolverPolarisation_label)
           if(iand(debug_level(debug_spectral),debug_levelBasic)/= 0 .and. worldrank == 0_pInt) &
           call IO_warning(42_pInt, ext_msg='debug Divergence')
           call Polarisation_init
         
         case default
           call IO_error(error_ID = 891, ext_msg = trim(spectral_solver))
       
       end select 
     
      case(FIELD_THERMAL_ID)
       call spectral_thermal_init
 
     case(FIELD_DAMAGE_ID)
       call spectral_damage_init()

   end select
 enddo
 
!--------------------------------------------------------------------------------------------------
! write header of output file
 if (worldrank == 0) then
   if (.not. appendToOutFile) then                                                                    ! after restart, append to existing results file
     open(newunit=resUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                                 '.spectralOut',form='UNFORMATTED',status='REPLACE')
     write(resUnit) 'load:',       trim(loadCaseFile)                                                 ! ... and write header
     write(resUnit) 'workingdir:', trim(getSolverWorkingDirectoryName())
     write(resUnit) 'geometry:',   trim(geometryFile)
     write(resUnit) 'grid:',       grid
     write(resUnit) 'size:',       geomSize
     write(resUnit) 'materialpoint_sizeResults:', materialpoint_sizeResults
     write(resUnit) 'loadcases:',  size(loadCases)
     write(resUnit) 'frequencies:', loadCases%outputfrequency                                         ! one entry per LoadCase
     write(resUnit) 'times:',      loadCases%time                                                     ! one entry per LoadCase
     write(resUnit) 'logscales:',  loadCases%logscale
     write(resUnit) 'increments:', loadCases%incs                                                     ! one entry per LoadCase
     write(resUnit) 'startingIncrement:', restartInc - 1_pInt                                         ! start with writing out the previous inc
     write(resUnit) 'eoh'    
     close(resUnit)                                                                                   ! end of header
     open(newunit=statUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                                 '.sta',form='FORMATTED',status='REPLACE')
     write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                 ! statistics file
     if (iand(debug_level(debug_spectral),debug_levelBasic) /= 0) &
       write(6,'(/,a)') ' header of result and statistics file written out'
     flush(6)
   else                                                                                             ! open new files ...
     open(newunit=statUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                                 '.sta',form='FORMATTED', position='APPEND', status='OLD')
   endif
 endif

!--------------------------------------------------------------------------------------------------
! prepare MPI parallel out (including opening of file)
 allocate(outputSize(worldsize), source = 0_MPI_OFFSET_KIND)
 outputSize(worldrank+1) = size(materialpoint_results,kind=MPI_OFFSET_KIND)*int(pReal,MPI_OFFSET_KIND)
 call MPI_allreduce(MPI_IN_PLACE,outputSize,worldsize,MPI_LONG,MPI_SUM,PETSC_COMM_WORLD,ierr)       ! get total output size over each process
 if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_allreduce')
 call MPI_file_open(PETSC_COMM_WORLD, &
                    trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut', &
                    MPI_MODE_WRONLY + MPI_MODE_APPEND, &
                    MPI_INFO_NULL, &
                    resUnit, &
                    ierr)
 if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_file_open')
 call MPI_file_get_position(resUnit,fileOffset,ierr)                                                ! get offset from header
 if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_file_get_position')
 fileOffset = fileOffset + sum(outputSize(1:worldrank))                                             ! offset of my process in file (header + processes before me)
 call MPI_file_seek (resUnit,fileOffset,MPI_SEEK_SET,ierr)
 if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_file_seek')

 if (.not. appendToOutFile) then                                                                    ! if not restarting, write 0th increment
   do i=1, size(materialpoint_results,3)/(maxByteOut/(materialpoint_sizeResults*pReal))+1           ! slice the output of my process in chunks not exceeding the limit for one output
     outputIndex=int([(i-1_pInt)*((maxByteOut/pReal)/materialpoint_sizeResults)+1_pInt, &
                      min(i*((maxByteOut/pReal)/materialpoint_sizeResults),size(materialpoint_results,3))],pLongInt)
     call MPI_file_write(resUnit,reshape(materialpoint_results(:,:,outputIndex(1):outputIndex(2)),&
                                   [(outputIndex(2)-outputIndex(1)+1)*materialpoint_sizeResults]), &
                         (outputIndex(2)-outputIndex(1)+1)*materialpoint_sizeResults,&
                         MPI_DOUBLE, MPI_STATUS_IGNORE, ierr)
     if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_file_write')
   enddo
   fileOffset = fileOffset + sum(outputSize)                                                        ! forward to current file position
   if (worldrank == 0) &
     write(6,'(1/,a)') ' ... writing initial configuration to file ........................'
 endif
!--------------------------------------------------------------------------------------------------
! loopping over loadcases
 loadCaseLooping: do currentLoadCase = 1_pInt, size(loadCases)
   time0 = time                                                                                     ! currentLoadCase start time                
   guess = loadCases(currentLoadCase)%followFormerTrajectory                                        ! change of load case? homogeneous guess for the first inc

!--------------------------------------------------------------------------------------------------
! loop oper incs defined in input file for current currentLoadCase
   incLooping: do inc = 1_pInt, loadCases(currentLoadCase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt

!--------------------------------------------------------------------------------------------------
! forwarding time
     timeIncOld = timeinc
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
         timeinc = time0 * &
              ( (1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real(          inc,pReal)/&
                                                    real(loadCases(currentLoadCase)%incs ,pReal))&
               -(1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                     real(loadCases(currentLoadCase)%incs ,pReal)))
       endif
     endif
     timeinc = timeinc / 2.0_pReal**real(cutBackLevel,pReal)                                        ! depending on cut back level, decrease time step

     forwarding: if(totalIncsCounter >= restartInc) then
       stepFraction = 0_pInt

!--------------------------------------------------------------------------------------------------
! loop over sub incs 
       subIncLooping: do while (stepFraction/subStepFactor**cutBackLevel <1_pInt)
         time = time + timeinc                                                                      ! forward time
         stepFraction = stepFraction + 1_pInt 
         remainingLoadCaseTime = time0 - time + loadCases(currentLoadCase)%time + timeInc
           
!--------------------------------------------------------------------------------------------------
! report begin of new increment
         if (worldrank == 0) then
           write(6,'(/,a)') ' ###########################################################################'
           write(6,'(1x,a,es12.5'//&
                   ',a,'//IO_intOut(inc)//',a,'//IO_intOut(loadCases(currentLoadCase)%incs)//&
                   ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(subStepFactor**cutBackLevel)//&
                   ',a,'//IO_intOut(currentLoadCase)//',a,'//IO_intOut(size(loadCases))//')') &
                   'Time', time, &
                   's: Increment ', inc, '/', loadCases(currentLoadCase)%incs,&
                   '-', stepFraction, '/', subStepFactor**cutBackLevel,&
                   ' of load case ', currentLoadCase,'/',size(loadCases)
           flush(6)
           write(incInfo,'(a,'//IO_intOut(totalIncsCounter)//',a,'//IO_intOut(sum(loadCases%incs))//&
                 ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(subStepFactor**cutBackLevel)//')') &
                 'Increment ',totalIncsCounter,'/',sum(loadCases%incs),&
                 '-',stepFraction, '/', subStepFactor**cutBackLevel
         endif     

!--------------------------------------------------------------------------------------------------
! forward fields
         do field = 1, nActiveFields
           select case(loadCases(currentLoadCase)%ID(field))
             case(FIELD_MECH_ID)
               select case (spectral_solver)
                 case (DAMASK_spectral_SolverBasicPETSc_label)
                   call BasicPETSc_forward (&
                       guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                       F_BC               = loadCases(currentLoadCase)%deformation, &
                       P_BC               = loadCases(currentLoadCase)%P, &
                       rotation_BC        = loadCases(currentLoadCase)%rotation)
                 case (DAMASK_spectral_SolverAL_label)
                   call AL_forward (&
                       guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                       F_BC               = loadCases(currentLoadCase)%deformation, &
                       P_BC               = loadCases(currentLoadCase)%P, &
                       rotation_BC        = loadCases(currentLoadCase)%rotation)
                 case (DAMASK_spectral_SolverPolarisation_label)
                   call Polarisation_forward (&
                       guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                       F_BC               = loadCases(currentLoadCase)%deformation, &
                       P_BC               = loadCases(currentLoadCase)%P, &
                       rotation_BC        = loadCases(currentLoadCase)%rotation)
               end select 
     
           case(FIELD_THERMAL_ID)
               call spectral_thermal_forward (&
                   guess,timeinc,timeIncOld,remainingLoadCaseTime)
                   
           case(FIELD_DAMAGE_ID)
               call spectral_damage_forward (&
                   guess,timeinc,timeIncOld,remainingLoadCaseTime)
           end select
         enddo       
           
!--------------------------------------------------------------------------------------------------
! solve fields
         stagIter = 0_pInt
         stagIterate = .true.
         do while (stagIterate)
           do field = 1, nActiveFields
             select case(loadCases(currentLoadCase)%ID(field))
               case(FIELD_MECH_ID)
                 select case (spectral_solver)
                   case (DAMASK_spectral_SolverBasicPETSc_label)
                     solres(field) = BasicPETSC_solution (&
                         incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                         P_BC               = loadCases(currentLoadCase)%P, &
                         F_BC               = loadCases(currentLoadCase)%deformation, &
                         rotation_BC        = loadCases(currentLoadCase)%rotation)
         
                   case (DAMASK_spectral_SolverAL_label)
                     solres(field) = AL_solution (&
                         incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                         P_BC               = loadCases(currentLoadCase)%P, &
                         F_BC               = loadCases(currentLoadCase)%deformation, &
                         rotation_BC        = loadCases(currentLoadCase)%rotation)
         
                   case (DAMASK_spectral_SolverPolarisation_label)
                     solres(field) = Polarisation_solution (&
                         incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                         P_BC               = loadCases(currentLoadCase)%P, &
                         F_BC               = loadCases(currentLoadCase)%deformation, &
                         rotation_BC        = loadCases(currentLoadCase)%rotation)
       
                 end select 
     
               case(FIELD_THERMAL_ID)
                 solres(field) = spectral_thermal_solution (&
                     guess,timeinc,timeIncOld,remainingLoadCaseTime)
 
               case(FIELD_DAMAGE_ID)
                 solres(field) = spectral_damage_solution (&
                     guess,timeinc,timeIncOld,remainingLoadCaseTime)

             end select
             if(.not. solres(field)%converged) exit                                                ! no solution found
           enddo
           stagIter = stagIter + 1_pInt
           stagIterate = stagIter < stagItMax .and. &
                         all(solres(:)%converged) .and. &
                         .not. all(solres(:)%stagConverged)
         enddo     

!--------------------------------------------------------------------------------------------------
! check solution 
         cutBack = .False.                                                                   
         if(solres(1)%termIll .or. .not. all(solres(:)%converged .and. solres(:)%stagConverged)) then ! no solution found
           if (cutBackLevel < maxCutBack) then                                                      ! do cut back
             if (worldrank == 0) write(6,'(/,a)') ' cut back detected'
             cutBack = .True.
             stepFraction = (stepFraction - 1_pInt) * subStepFactor                                 ! adjust to new denominator
             cutBackLevel = cutBackLevel + 1_pInt
             time    = time - timeinc                                                               ! rewind time
             timeinc = timeinc/2.0_pReal
           elseif (solres(1)%termIll) then                                                          ! material point model cannot find a solution, exit in any casy
             call IO_warning(850_pInt)
             call quit(-1_pInt*(lastRestartWritten+1_pInt))                                         ! quit and provide information about last restart inc written (e.g. for regridding)
           elseif (continueCalculation == 1_pInt)  then
             guess = .true.                                                                         ! accept non converged BVP solution   
           else                                                                                     ! default behavior, exit if spectral solver does not converge                                
             call IO_warning(850_pInt)
             call quit(-1_pInt*(lastRestartWritten+1_pInt))                                         ! quit and provide information about last restart inc written (e.g. for regridding)
           endif
         else
           guess = .true.                                                                           ! start guessing after first converged (sub)inc
         endif
         if (.not. cutBack) then
           if (worldrank == 0) then
             write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                               solres%converged, solres%iterationsNeeded                            ! write statistics about accepted solution
             flush(statUnit)
           endif 
         endif  
       enddo subIncLooping
       cutBackLevel = max(0_pInt, cutBackLevel - 1_pInt)                                            ! try half number of subincs next inc
       if(all(solres(:)%converged)) then                                                            ! report converged inc
         convergedCounter = convergedCounter + 1_pInt
         if (worldrank == 0) &
           write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &
                                     ' increment ', totalIncsCounter, ' converged'
       else
         if (worldrank == 0) &
           write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &                                   ! report non-converged inc
                                     ' increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       endif; flush(6)
       if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0_pInt) then                      ! at output frequency
         if (worldrank == 0) &
           write(6,'(1/,a)') ' ... writing results to file ......................................'
         call materialpoint_postResults()
         call MPI_file_seek (resUnit,fileOffset,MPI_SEEK_SET,ierr)
         if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_file_seek')
         do i=1, size(materialpoint_results,3)/(maxByteOut/(materialpoint_sizeResults*pReal))+1     ! slice the output of my process in chunks not exceeding the limit for one output
           outputIndex=int([(i-1_pInt)*((maxByteOut/pReal)/materialpoint_sizeResults)+1_pInt, &
                      min(i*((maxByteOut/pReal)/materialpoint_sizeResults),size(materialpoint_results,3))],pLongInt)
           call MPI_file_write(resUnit,reshape(materialpoint_results(:,:,outputIndex(1):outputIndex(2)),&
                                         [(outputIndex(2)-outputIndex(1)+1)*materialpoint_sizeResults]), &
                               (outputIndex(2)-outputIndex(1)+1)*materialpoint_sizeResults,&
                               MPI_DOUBLE, MPI_STATUS_IGNORE, ierr)
           if(ierr /=0_pInt) call IO_error(894_pInt, ext_msg='MPI_file_write')
         enddo
         fileOffset = fileOffset + sum(outputSize)                                                  ! forward to current file position
       endif
       if( loadCases(currentLoadCase)%restartFrequency > 0_pInt .and. &                             ! at frequency of writing restart information set restart parameter for FEsolving 
                      mod(inc,loadCases(currentLoadCase)%restartFrequency) == 0_pInt) then          ! first call to CPFEM_general will write? 
         restartWrite = .true.
         lastRestartWritten = inc
       endif 
     else forwarding
       time = time + timeinc
       guess = .true.
     endif forwarding

    enddo incLooping
 enddo loadCaseLooping

!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
 if (worldrank == 0) then
   write(6,'(/,a)') ' ###########################################################################'
   write(6,'(1x,i6.6,a,i6.6,a,f5.1,a)') convergedCounter, ' out of ', &
                                     notConvergedCounter + convergedCounter, ' (', &
                                     real(convergedCounter, pReal)/&
                                     real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, &
                                     ' %) increments converged!'
 endif
 call MPI_file_close(resUnit,ierr)
 close(statUnit)

 do field = 1, nActiveFields
   select case(loadCases(1)%ID(field))
     case(FIELD_MECH_ID)
       select case (spectral_solver)
         case (DAMASK_spectral_SolverBasicPETSc_label)
           call BasicPETSC_destroy()
         case (DAMASK_spectral_SolverAL_label)
           call AL_destroy()
         case (DAMASK_spectral_SolverPolarisation_label)
           call Polarisation_destroy()
       end select 
     case(FIELD_THERMAL_ID)
       call spectral_thermal_destroy()
     case(FIELD_DAMAGE_ID)
       call spectral_damage_destroy()
   end select
 enddo
 call utilities_destroy()

 call PETScFinalize(ierr); CHKERRQ(ierr)

 if (notConvergedCounter > 0_pInt) call quit(3_pInt)                                                ! error if some are not converged
 call quit(0_pInt)                                                                                  ! no complains ;)

end program DAMASK_spectral


!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief quit subroutine to mimic behavior of FEM solvers
!> @details exits the Spectral solver and reports time and duration. Exit code 0 signals
!> everything went fine. Exit code 1 signals an error, message according to IO_error. Exit code 
!> 2 signals request for regridding, increment of last saved restart information is written to
!> stderr. Exit code 3 signals no severe problems, but some increments did not converge
!--------------------------------------------------------------------------------------------------
subroutine quit(stop_id)
 use prec, only: &
   pInt
 use numerics, only: &
   worldrank  

 implicit none
 integer(pInt), intent(in) :: stop_id
 integer, dimension(8) :: dateAndTime                                                               ! type default integer

 if (worldrank == 0_pInt) then
   call date_and_time(values = dateAndTime)
   write(6,'(/,a)') 'DAMASK terminated on:'
   write(6,'(a,2(i2.2,a),i4.4)') 'Date:               ',dateAndTime(3),'/',&
                                                        dateAndTime(2),'/',&
                                                        dateAndTime(1)
   write(6,'(a,2(i2.2,a),i2.2)') 'Time:               ',dateAndTime(5),':',&
                                                        dateAndTime(6),':',&
                                                        dateAndTime(7)
 endif
 
 if (stop_id == 0_pInt) stop 0                                                                      ! normal termination
 if (stop_id <  0_pInt) then                                                                        ! trigger regridding
   if (worldrank == 0_pInt) &
     write(0,'(a,i6)') 'restart information available at ', stop_id*(-1_pInt)
   stop 2
 endif
 if (stop_id == 3_pInt) stop 3                                                                      ! not all incs converged
 stop 1                                                                                             ! error (message from IO_error)

end subroutine quit
