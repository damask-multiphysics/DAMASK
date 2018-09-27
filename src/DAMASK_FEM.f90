!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the FEM solver
!> @details doing cutbacking, forwarding in case of restart, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_FEM 
 use, intrinsic :: &
   iso_fortran_env                                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use prec, only: &
   pInt, &
   pReal, &
   tol_math_check
 use DAMASK_interface, only: &
   DAMASK_interface_init, &
   loadCaseFile, &
   getSolverJobName
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
 use CPFEM2, only: &
   CPFEM_initAll
 use FEsolving, only: &
   restartWrite, &
   restartInc
 use numerics, only: &
   worldrank, &
   maxCutBack, &
   stagItMax
 use mesh, only: &
   mesh_Nboundaries, &
   mesh_boundaries, &   
   geomMesh
 use FEM_Utilities, only: &
   utilities_init, &
   tSolutionState, &
   tLoadCase, &
   cutBack, &
   maxFields, &
   nActiveFields, &
   FIELD_MECH_ID, &
   COMPONENT_MECH_X_ID, &
   COMPONENT_MECH_Y_ID, &
   COMPONENT_MECH_Z_ID, &
   FIELD_MECH_label
 use FEM_mech
 
 implicit none
#include <petsc/finclude/petsc.h>

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 integer(pInt), parameter  :: FILEUNIT           = 234_pInt                                         !< file unit, DAMASK IO does not support newunit feature
 integer(pInt), allocatable, dimension(:) :: chunkPos                                               ! this is longer than needed for geometry parsing
 
 integer(pInt) :: &
   N_def = 0_pInt                                                                                   !< # of rate of deformation specifiers found in load case file
 character(len=65536) :: &
   line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.

 integer(pInt), parameter :: &
   subStepFactor = 2_pInt                                                                           !< for each substep, divide the last time increment by 2.0
 real(pReal) :: &
   time = 0.0_pReal, &                                                                              !< elapsed time
   time0 = 0.0_pReal, &                                                                             !< begin of interval
   timeinc = 0.0_pReal, &                                                                           !< current time interval
   timeIncOld = 0.0_pReal, &                                                                        !< previous time interval
   remainingLoadCaseTime = 0.0_pReal                                                                !< remaining time of current load case
 logical :: &
   guess                                                                                            !< guess along former trajectory
 integer(pInt) :: &
   i, &
   errorID, &
   cutBackLevel = 0_pInt, &                                                                         !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
   stepFraction = 0_pInt                                                                            !< fraction of current time interval
 integer(pInt) :: &
   currentLoadcase = 0_pInt, &                                                                      !< current load case
   currentFace = 0_pInt, &
   inc, &                                                                                           !< current increment in current load case
   totalIncsCounter = 0_pInt, &                                                                     !< total No. of increments
   convergedCounter = 0_pInt, &                                                                     !< No. of converged increments
   notConvergedCounter = 0_pInt, &                                                                  !< No. of non-converged increments
   statUnit = 0_pInt, &                                                                             !< file unit for statistics output
   lastRestartWritten = 0_pInt                                                                      !< total increment No. at which last restart information was written
 integer(pInt) :: &
   stagIter, &
   component
 logical :: &
   stagIterate  
 character(len=6)  :: loadcase_string
 character(len=1024)  :: incInfo                                                                    !< string parsed to solution with information about current load case
 type(tLoadCase), allocatable, dimension(:) :: loadCases                                            !< array of all load cases
 type(tSolutionState), allocatable, dimension(:) :: solres
 PetscInt :: faceSet, currentFaceSet
 PetscInt :: field, dimPlex
 PetscErrorCode :: ierr

 external :: &
   MPI_abort, &
   quit

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
 call CPFEM_initAll(el = 1_pInt, ip = 1_pInt)
   write(6,'(/,a)')   ' <<<+-  DAMASK_FEM init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
! reading basic information from load case file and allocate data structure containing load cases
 call DMGetDimension(geomMesh,dimPlex,ierr)! CHKERRQ(ierr)                                            !< dimension of mesh (2D or 3D)
 nActiveFields = 1
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
   do i = 1_pInt, chunkPos(1)                                                                      ! reading compulsory parameters for loadcase
     select case (IO_lc(IO_stringValue(line,chunkPos,i)))
       case('$loadcase')
         N_def = N_def + 1_pInt
     end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

 allocate (loadCases(N_def))         

 do i = 1, size(loadCases)
   allocate(loadCases(i)%fieldBC(nActiveFields))
   field = 1
   loadCases(i)%fieldBC(field)%ID = FIELD_MECH_ID
 enddo

 do i = 1, size(loadCases)
   do field = 1, nActiveFields
     select case (loadCases(i)%fieldBC(field)%ID)
       case(FIELD_MECH_ID)
         loadCases(i)%fieldBC(field)%nComponents = dimPlex                                        !< X, Y (, Z) displacements
         allocate(loadCases(i)%fieldBC(field)%componentBC(loadCases(i)%fieldBC(field)%nComponents))
         do component = 1, loadCases(i)%fieldBC(field)%nComponents
           select case (component)
             case (1)
               loadCases(i)%fieldBC(field)%componentBC(component)%ID = COMPONENT_MECH_X_ID
             case (2)
               loadCases(i)%fieldBC(field)%componentBC(component)%ID = COMPONENT_MECH_Y_ID
             case (3)
               loadCases(i)%fieldBC(field)%componentBC(component)%ID = COMPONENT_MECH_Z_ID
           end select
         enddo  
     end select
     do component = 1, loadCases(i)%fieldBC(field)%nComponents
       allocate(loadCases(i)%fieldBC(field)%componentBC(component)%Value(mesh_Nboundaries), source = 0.0_pReal)
       allocate(loadCases(i)%fieldBC(field)%componentBC(component)%Mask (mesh_Nboundaries), source = .false.)
     enddo
   enddo       
 enddo   

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(FILEUNIT)
 do
   line = IO_read(FILEUNIT)
   if (trim(line) == IO_EOF) exit
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   chunkPos = IO_stringPos(line)
   do i = 1_pInt, chunkPos(1)
     select case (IO_lc(IO_stringValue(line,chunkPos,i)))
!--------------------------------------------------------------------------------------------------
! loadcase information
       case('$loadcase')
         currentLoadCase = IO_intValue(line,chunkPos,i+1_pInt)
       case('face')
         currentFace = IO_intValue(line,chunkPos,i+1_pInt)
         currentFaceSet = -1_pInt
         do faceSet = 1, mesh_Nboundaries
           if (mesh_boundaries(faceSet) == currentFace) currentFaceSet = faceSet
         enddo
         if (currentFaceSet < 0_pInt) call IO_error(error_ID = errorID, ext_msg = 'invalid BC')  
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

!--------------------------------------------------------------------------------------------------
! boundary condition information
       case('x')                                                                                    ! X displacement field  
         do field = 1, nActiveFields
           if (loadCases(currentLoadCase)%fieldBC(field)%ID == FIELD_MECH_ID) then
             do component = 1, loadcases(currentLoadCase)%fieldBC(field)%nComponents
               if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%ID == COMPONENT_MECH_X_ID) then
                 loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask (currentFaceSet) = &
                     .true.
                 loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Value(currentFaceSet) = &
                     IO_floatValue(line,chunkPos,i+1_pInt)
               endif
             enddo
           endif
         enddo  
       case('y')                                                                                    ! Y displacement field                                                                            
         do field = 1, nActiveFields
           if (loadCases(currentLoadCase)%fieldBC(field)%ID == FIELD_MECH_ID) then
             do component = 1, loadcases(currentLoadCase)%fieldBC(field)%nComponents
               if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%ID == COMPONENT_MECH_Y_ID) then
                 loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask (currentFaceSet) = &
                     .true.
                 loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Value(currentFaceSet) = &
                     IO_floatValue(line,chunkPos,i+1_pInt)
               endif
             enddo
           endif
         enddo  
       case('z')                                                                                    ! Z displacement field                                                                       
         do field = 1, nActiveFields
           if (loadCases(currentLoadCase)%fieldBC(field)%ID == FIELD_MECH_ID) then
             do component = 1, loadcases(currentLoadCase)%fieldBC(field)%nComponents
               if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%ID == COMPONENT_MECH_Z_ID) then
                 loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask (currentFaceSet) = &
                     .true.
                 loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Value(currentFaceSet) = &
                     IO_floatValue(line,chunkPos,i+1_pInt)
               endif
             enddo
           endif
         enddo  
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
     do field = 1_pInt, nActiveFields
       select case (loadCases(currentLoadCase)%fieldBC(field)%ID)
         case(FIELD_MECH_ID)
           write(6,'(2x,a)') 'Field '//trim(FIELD_MECH_label)
       
       end select
       do faceSet = 1_pInt, mesh_Nboundaries
          do component = 1_pInt, loadCases(currentLoadCase)%fieldBC(field)%nComponents
            if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask(faceSet)) &
              write(6,'(4x,a,i2,a,i2,a,f12.7)') 'Face  ', mesh_boundaries(faceSet), &
                                                ' Component ', component, & 
                                                ' Value ', loadCases(currentLoadCase)%fieldBC(field)% &
                                                             componentBC(component)%Value(faceSet)
          enddo
        enddo       
     enddo
     write(6,'(2x,a,f12.6)') 'time:       ', loadCases(currentLoadCase)%time
     if (loadCases(currentLoadCase)%incs < 1_pInt)             errorID = 835_pInt                     ! non-positive incs count
     write(6,'(2x,a,i5)')    'increments: ', loadCases(currentLoadCase)%incs
     if (loadCases(currentLoadCase)%outputfrequency < 1_pInt)  errorID = 836_pInt                     ! non-positive result frequency
     write(6,'(2x,a,i5)')    'output  frequency:  ', &
                loadCases(currentLoadCase)%outputfrequency
     write(6,'(2x,a,i5,/)')    'restart frequency:  ', &
                loadCases(currentLoadCase)%restartfrequency
     if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)               ! exit with error message
   enddo checkLoadcases
 endif

!--------------------------------------------------------------------------------------------------
! doing initialization depending on selected solver 
 call Utilities_init()
 do field = 1, nActiveFields
   select case (loadCases(1)%fieldBC(field)%ID)
     case(FIELD_MECH_ID)
       call FEM_mech_init(loadCases(1)%fieldBC(field))
   end select
 enddo   

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
           select case (loadCases(currentLoadCase)%fieldBC(field)%ID)
             case(FIELD_MECH_ID)
               call FEM_mech_forward (&
                   guess,timeinc,timeIncOld,loadCases(currentLoadCase)%fieldBC(field))

          end select
         enddo       
           
!--------------------------------------------------------------------------------------------------
! solve fields
         stagIter = 0_pInt
         stagIterate = .true.
         do while (stagIterate)
           do field = 1, nActiveFields
             select case (loadCases(currentLoadCase)%fieldBC(field)%ID)
               case(FIELD_MECH_ID)
                 solres(field) = FEM_mech_solution (&
                       incInfo,timeinc,timeIncOld,loadCases(currentLoadCase)%fieldBC(field))

             end select
             if(.not. solres(field)%converged) exit                                                ! no solution found
           enddo
           stagIter = stagIter + 1_pInt
           stagIterate = stagIter < stagItMax .and. &
                         all(solres(:)%converged) .and. &
                         .not. all(solres(:)%stagConverged)
         enddo     
           
! check solution 
         cutBack = .False.                                                                   
         if(.not. all(solres(:)%converged .and. solres(:)%stagConverged)) then                      ! no solution found
           if (cutBackLevel < maxCutBack) then                                                   ! do cut back
             if (worldrank == 0) &
               write(6,'(/,a)') ' cut back detected'
             cutBack = .True.
             stepFraction = (stepFraction - 1_pInt) * subStepFactor                                 ! adjust to new denominator
             cutBackLevel = cutBackLevel + 1_pInt
             time    = time - timeinc                                                               ! rewind time
             timeinc = timeinc/2.0_pReal
           else                                                                                     ! default behavior, exit if spectral solver does not converge                                
             call IO_warning(850_pInt)
             call quit(-1_pInt*(lastRestartWritten+1_pInt))                                         ! quit and provide information about last restart inc written (e.g. for regridding)                                                                                       ! continue from non-converged solution and start guessing after accepted (sub)inc
           endif
         else
           guess = .true.                                                                           ! start guessing after first converged (sub)inc
           timeIncOld = timeinc
         endif
         if (.not. cutBack) then
           if (worldrank == 0)  write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                             solres%converged, solres%iterationsNeeded                              ! write statistics about accepted solution
         endif
       enddo subIncLooping
       cutBackLevel = max(0_pInt, cutBackLevel - 1_pInt)                                            ! try half number of subincs next inc
       if(all(solres(:)%converged)) then                                                                    ! report converged inc
         convergedCounter = convergedCounter + 1_pInt
         if (worldrank == 0) then
         write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &
                                     ' increment ', totalIncsCounter, ' converged'
         endif
       else
         if (worldrank == 0) then
         write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &                                     ! report non-converged inc
                                     ' increment ', totalIncsCounter, ' NOT converged'
         endif
         notConvergedCounter = notConvergedCounter + 1_pInt
       endif; flush(6)
       if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0_pInt) then                      ! at output frequency
         if (worldrank == 0) then
         write(6,'(1/,a)') ' ... writing results to file ......................................'
         endif
       endif
       if( loadCases(currentLoadCase)%restartFrequency > 0_pInt .and. &                             ! at frequency of writing restart information set restart parameter for FEsolving 
                      mod(inc,loadCases(currentLoadCase)%restartFrequency) == 0_pInt) then          ! ToDo first call to CPFEM_general will write? 
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
 if (notConvergedCounter > 0_pInt) call quit(2_pInt)                                                ! error if some are not converged
 call quit(0_pInt)                                                                                  ! no complains ;)

end program DAMASK_FEM
