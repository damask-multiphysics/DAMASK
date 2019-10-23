!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the FEM solver
!> @details doing cutbacking, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_FEM 
#include <petsc/finclude/petscsys.h>
  use PetscDM
  use prec
  use DAMASK_interface
  use IO
  use math
  use CPFEM2
  use FEsolving
  use numerics
  use mesh
  use FEM_Utilities
  use FEM_mech
  
  implicit none

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
  integer, allocatable, dimension(:) :: chunkPos                                                    ! this is longer than needed for geometry parsing
  integer :: &
    N_def = 0                                                                                       !< # of rate of deformation specifiers found in load case file
  character(len=65536) :: &
    line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
  integer, parameter :: &
    subStepFactor = 2                                                                               !< for each substep, divide the last time increment by 2.0
  real(pReal) :: &
    time = 0.0_pReal, &                                                                              !< elapsed time
    time0 = 0.0_pReal, &                                                                             !< begin of interval
    timeinc = 0.0_pReal, &                                                                           !< current time interval
    timeIncOld = 0.0_pReal, &                                                                        !< previous time interval
    remainingLoadCaseTime = 0.0_pReal                                                                !< remaining time of current load case
  logical :: &
    guess, &                                                                                         !< guess along former trajectory
    stagIterate
  integer :: &
    i, &
    errorID, &
    cutBackLevel = 0, &                                                                             !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
    stepFraction = 0                                                                                !< fraction of current time interval
  integer :: &
    currentLoadcase = 0, &                                                                          !< current load case
    currentFace = 0, &
    inc, &                                                                                          !< current increment in current load case
    totalIncsCounter = 0, &                                                                         !< total # of increments
    convergedCounter = 0, &                                                                         !< # of converged increments
    notConvergedCounter = 0, &                                                                      !< # of non-converged increments
    fileUnit = 0, &                                                                                 !< file unit for reading load case and writing results
    myStat, &
    statUnit = 0, &                                                                                 !< file unit for statistics output
    stagIter, &
    component  
  character(len=6)  :: loadcase_string
  character(len=1024) :: &
    incInfo
  type(tLoadCase), allocatable, dimension(:) :: loadCases                                            !< array of all load cases
  type(tSolutionState), allocatable, dimension(:) :: solres
  PetscInt :: faceSet, currentFaceSet
  PetscInt :: field, dimPlex
  PetscErrorCode :: ierr
 
  external :: &
    quit

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
  call CPFEM_initAll
  write(6,'(/,a)')   ' <<<+-  DAMASK_FEM init  -+>>>'
 
! reading basic information from load case file and allocate data structure containing load cases
  call DMGetDimension(geomMesh,dimPlex,ierr); CHKERRA(ierr)                                         !< dimension of mesh (2D or 3D)
  nActiveFields = 1
  allocate(solres(nActiveFields))

!--------------------------------------------------------------------------------------------------
! reading basic information from load case file and allocate data structure containing load cases
  open(newunit=fileunit,iostat=myStat,file=trim(loadCaseFile),action='read')
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=trim(loadCaseFile))
  do
    read(fileUnit, '(A)', iostat=myStat) line
    if ( myStat /= 0) exit
    if (IO_isBlank(line)) cycle                                                                     ! skip empty lines
 
    chunkPos = IO_stringPos(line)
    do i = 1, chunkPos(1)                                                                           ! reading compulsory parameters for loadcase
      select case (IO_lc(IO_stringValue(line,chunkPos,i)))
        case('$loadcase')
          N_def = N_def + 1
      end select
    enddo                                                                                           ! count all identifiers to allocate memory and do sanity check
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
          loadCases(i)%fieldBC(field)%nComponents = dimPlex                                         !< X, Y (, Z) displacements
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
  rewind(fileUnit)
  do
    read(fileUnit, '(A)', iostat=myStat) line
    if ( myStat /= 0) exit
    if (IO_isBlank(line)) cycle                                                                     ! skip empty lines
 
    chunkPos = IO_stringPos(line)
    do i = 1, chunkPos(1)
      select case (IO_lc(IO_stringValue(line,chunkPos,i)))
!--------------------------------------------------------------------------------------------------
! loadcase information
        case('$loadcase')
          currentLoadCase = IO_intValue(line,chunkPos,i+1)
        case('face')
          currentFace = IO_intValue(line,chunkPos,i+1)
          currentFaceSet = -1
          do faceSet = 1, mesh_Nboundaries
            if (mesh_boundaries(faceSet) == currentFace) currentFaceSet = faceSet
          enddo
          if (currentFaceSet < 0) call IO_error(error_ID = errorID, ext_msg = 'invalid BC')  
        case('t','time','delta')                                                                    ! increment time
          loadCases(currentLoadCase)%time = IO_floatValue(line,chunkPos,i+1)
        case('n','incs','increments','steps')                                                       ! number of increments
          loadCases(currentLoadCase)%incs = IO_intValue(line,chunkPos,i+1)
        case('logincs','logincrements','logsteps')                                                  ! number of increments (switch to log time scaling)
          loadCases(currentLoadCase)%incs = IO_intValue(line,chunkPos,i+1)
          loadCases(currentLoadCase)%logscale = 1
        case('freq','frequency','outputfreq')                                                       ! frequency of result writings
          loadCases(currentLoadCase)%outputfrequency = IO_intValue(line,chunkPos,i+1)                
        case('guessreset','dropguessing')
          loadCases(currentLoadCase)%followFormerTrajectory = .false.                                ! do not continue to predict deformation along former trajectory

!--------------------------------------------------------------------------------------------------
! boundary condition information
        case('x')                                                                                   ! X displacement field  
          do field = 1, nActiveFields
            if (loadCases(currentLoadCase)%fieldBC(field)%ID == FIELD_MECH_ID) then
              do component = 1, loadcases(currentLoadCase)%fieldBC(field)%nComponents
                if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%ID == COMPONENT_MECH_X_ID) then
                  loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask (currentFaceSet) = &
                      .true.
                  loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Value(currentFaceSet) = &
                      IO_floatValue(line,chunkPos,i+1)
                endif
              enddo
            endif
          enddo  
        case('y')                                                                                   ! Y displacement field                                                                            
          do field = 1, nActiveFields
            if (loadCases(currentLoadCase)%fieldBC(field)%ID == FIELD_MECH_ID) then
              do component = 1, loadcases(currentLoadCase)%fieldBC(field)%nComponents
                if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%ID == COMPONENT_MECH_Y_ID) then
                  loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask (currentFaceSet) = &
                      .true.
                  loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Value(currentFaceSet) = &
                      IO_floatValue(line,chunkPos,i+1)
                endif
              enddo
            endif
          enddo  
        case('z')                                                                                   ! Z displacement field                                                                       
          do field = 1, nActiveFields
            if (loadCases(currentLoadCase)%fieldBC(field)%ID == FIELD_MECH_ID) then
              do component = 1, loadcases(currentLoadCase)%fieldBC(field)%nComponents
                if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%ID == COMPONENT_MECH_Z_ID) then
                  loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask (currentFaceSet) = &
                      .true.
                  loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Value(currentFaceSet) = &
                      IO_floatValue(line,chunkPos,i+1)
                endif
              enddo
            endif
          enddo  
      end select
  enddo; enddo
  close(fileUnit)
 
!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
  loadCases(1)%followFormerTrajectory = .false.                                                      ! cannot guess along trajectory for first inc of first currentLoadCase
  errorID = 0
  checkLoadcases: do currentLoadCase = 1, size(loadCases)
      write (loadcase_string, '(i6)' ) currentLoadCase
      write(6,'(1x,a,i6)') 'load case: ', currentLoadCase
      if (.not. loadCases(currentLoadCase)%followFormerTrajectory) &
        write(6,'(2x,a)') 'drop guessing along trajectory'
      do field = 1, nActiveFields
        select case (loadCases(currentLoadCase)%fieldBC(field)%ID)
          case(FIELD_MECH_ID)
            write(6,'(2x,a)') 'Field '//trim(FIELD_MECH_label)
        
        end select
        do faceSet = 1, mesh_Nboundaries
           do component = 1, loadCases(currentLoadCase)%fieldBC(field)%nComponents
             if (loadCases(currentLoadCase)%fieldBC(field)%componentBC(component)%Mask(faceSet)) &
               write(6,'(4x,a,i2,a,i2,a,f12.7)') 'Face  ', mesh_boundaries(faceSet), &
                                                 ' Component ', component, & 
                                                 ' Value ', loadCases(currentLoadCase)%fieldBC(field)% &
                                                              componentBC(component)%Value(faceSet)
           enddo
         enddo       
      enddo
      write(6,'(2x,a,f12.6)') 'time:       ', loadCases(currentLoadCase)%time
      if (loadCases(currentLoadCase)%incs < 1)             errorID = 835                            ! non-positive incs count
      write(6,'(2x,a,i5)')    'increments: ', loadCases(currentLoadCase)%incs
      if (loadCases(currentLoadCase)%outputfrequency < 1)  errorID = 836                            ! non-positive result frequency
      write(6,'(2x,a,i5)')    'output  frequency:  ', &
                 loadCases(currentLoadCase)%outputfrequency
      if (errorID > 0) call IO_error(error_ID = errorID, ext_msg = loadcase_string)                 ! exit with error message
  enddo checkLoadcases

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call Utilities_init
  do field = 1, nActiveFields
    select case (loadCases(1)%fieldBC(field)%ID)
      case(FIELD_MECH_ID)
        call FEM_mech_init(loadCases(1)%fieldBC(field))
    end select
  enddo   

  open(newunit=statUnit,file=trim(getSolverJobName())//'.sta',form='FORMATTED',status='REPLACE')
  write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                 ! statistics file

  loadCaseLooping: do currentLoadCase = 1, size(loadCases)
    time0 = time                                                                                    ! load case start time
    guess = loadCases(currentLoadCase)%followFormerTrajectory                                       ! change of load case? homogeneous guess for the first inc
 
    incLooping: do inc = 1, loadCases(currentLoadCase)%incs
      totalIncsCounter = totalIncsCounter + 1                                                 

!--------------------------------------------------------------------------------------------------
! forwarding time
      timeIncOld = timeinc                                                                          ! last timeinc that brought former inc to an end
      if (loadCases(currentLoadCase)%logscale == 0) then                                            ! linear scale
        timeinc = loadCases(currentLoadCase)%time/real(loadCases(currentLoadCase)%incs,pReal)
      else
        if (currentLoadCase == 1) then                                                              ! 1st load case of logarithmic scale
          if (inc == 1) then                                                                        ! 1st inc of 1st load case of logarithmic scale
            timeinc = loadCases(1)%time*(2.0_pReal**real(    1-loadCases(1)%incs ,pReal))           ! assume 1st inc is equal to 2nd 
          else                                                                                      ! not-1st inc of 1st load case of logarithmic scale
            timeinc = loadCases(1)%time*(2.0_pReal**real(inc-1-loadCases(1)%incs ,pReal))
          endif
        else                                                                                        ! not-1st load case of logarithmic scale
          timeinc = time0 * &
               ( (1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real(          inc,pReal)/&
                                                     real(loadCases(currentLoadCase)%incs ,pReal))&
                -(1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real( inc-1  ,pReal)/&
                                                      real(loadCases(currentLoadCase)%incs ,pReal)))
        endif
      endif
      timeinc = timeinc * real(subStepFactor,pReal)**real(-cutBackLevel,pReal)                      ! depending on cut back level, decrease time step


      stepFraction = 0                                                                            ! fraction scaled by stepFactor**cutLevel

      subStepLooping: do while (stepFraction < subStepFactor**cutBackLevel)
        remainingLoadCaseTime = loadCases(currentLoadCase)%time+time0 - time
        time = time + timeinc                                                                     ! forward target time
        stepFraction = stepFraction + 1                                                           ! count step
         
!--------------------------------------------------------------------------------------------------
! report begin of new step
        write(6,'(/,a)') ' ###########################################################################'
        write(6,'(1x,a,es12.5'//&
                ',a,'//IO_intOut(inc)//',a,'//IO_intOut(loadCases(currentLoadCase)%incs)//&
                ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(subStepFactor**cutBackLevel)//&
                ',a,'//IO_intOut(currentLoadCase)//',a,'//IO_intOut(size(loadCases))//')') &
                'Time', time, &
                's: Increment ', inc, '/', loadCases(currentLoadCase)%incs,&
                '-', stepFraction, '/', subStepFactor**cutBackLevel,&
                ' of load case ', currentLoadCase,'/',size(loadCases)
        write(incInfo,&
               '(a,'//IO_intOut(totalIncsCounter)//&
               ',a,'//IO_intOut(sum(loadCases%incs))//&
               ',a,'//IO_intOut(stepFraction)//&
               ',a,'//IO_intOut(subStepFactor**cutBackLevel)//')') &
               'Increment ',totalIncsCounter,'/',sum(loadCases%incs),&
               '-',stepFraction, '/', subStepFactor**cutBackLevel
        flush(6)

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
        stagIter = 0
        stagIterate = .true.
        do while (stagIterate)
          do field = 1, nActiveFields
            select case (loadCases(currentLoadCase)%fieldBC(field)%ID)
              case(FIELD_MECH_ID)
                solres(field) = FEM_mech_solution (&
                      incInfo,timeinc,timeIncOld,loadCases(currentLoadCase)%fieldBC(field))

            end select

            if(.not. solres(field)%converged) exit                                                  ! no solution found

          enddo
          stagIter = stagIter + 1
          stagIterate =            stagIter < stagItMax &
                       .and.       all(solres(:)%converged) &
                       .and. .not. all(solres(:)%stagConverged)                                     ! stationary with respect to staggered iteration
        enddo     
         
! check solution 
        cutBack = .False.                                                                   
        if(.not. all(solres(:)%converged .and. solres(:)%stagConverged)) then                       ! no solution found
          if (cutBackLevel < maxCutBack) then                                                       ! do cut back
            write(6,'(/,a)') ' cut back detected'
            cutBack = .True.
            stepFraction = (stepFraction - 1) * subStepFactor                                       ! adjust to new denominator
            cutBackLevel = cutBackLevel + 1
            time    = time - timeinc                                                                ! rewind time
            timeinc = timeinc/2.0_pReal
          else                                                                                      ! default behavior, exit if spectral solver does not converge                                
            call IO_warning(850)
            call quit(1)                                                                            ! quit
          endif
        else
          guess = .true.                                                                            ! start guessing after first converged (sub)inc
          timeIncOld = timeinc
        endif
        if (.not. cutBack) then
          if (worldrank == 0)  write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                            solres%converged, solres%iterationsNeeded                               ! write statistics about accepted solution
        endif
      enddo subStepLooping

      cutBackLevel = max(0, cutBackLevel - 1)                                                       ! try half number of subincs next inc

      if (all(solres(:)%converged)) then
        convergedCounter = convergedCounter + 1
        write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &                                      ! report converged inc
                                  ' increment ', totalIncsCounter, ' converged'
      else
        notConvergedCounter = notConvergedCounter + 1
        write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &                                      ! report non-converged inc
                                  ' increment ', totalIncsCounter, ' NOT converged'
      endif; flush(6)

      if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0) then                            ! at output frequency
        write(6,'(1/,a)') ' ... writing results to file ......................................'
        call CPFEM_results(totalIncsCounter,time)
      endif


    enddo incLooping

  enddo loadCaseLooping
 
 
!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  write(6,'(/,a)') ' ###########################################################################'
  write(6,'(1x,'//IO_intOut(convergedCounter)//',a,'//IO_intOut(notConvergedCounter + convergedCounter)//',a,f5.1,a)') &
    convergedCounter, ' out of ', &
    notConvergedCounter + convergedCounter, ' (', &
    real(convergedCounter, pReal)/&
    real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, ' %) increments converged!'
  flush(6)
  close(statUnit)

  if (notConvergedCounter > 0) call quit(2)                                                         ! error if some are not converged
  call quit(0)                                                                                      ! no complains ;)

end program DAMASK_FEM
