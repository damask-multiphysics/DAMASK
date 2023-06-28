!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the FEM solver
!> @details doing cutbacking, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_mesh
#include <petsc/finclude/petscsys.h>
  use PetscDM
  use prec
  use CLI
  use parallelization
  use IO
  use math
  use materialpoint
  use config
  use discretization_mesh
  use FEM_Utilities
  use mesh_mechanical_FEM

  implicit none(type,external)

  type :: tLoadCase
    real(pREAL)  :: time                   = 0.0_pREAL                                              !< length of increment
    integer      :: incs                   = 0, &                                                   !< number of increments
                    outputfrequency        = 1                                                      !< frequency of result writes
    logical      :: followFormerTrajectory = .true.                                                 !< follow trajectory of former loadcase
    integer,        allocatable, dimension(:) :: faceID
    type(tFieldBC), allocatable, dimension(:) :: fieldBC
  end type tLoadCase

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
  integer, allocatable, dimension(:) :: chunkPos                                                    ! this is longer than needed for geometry parsing
  integer :: &
    N_def = 0                                                                                       !< # of rate of deformation specifiers found in load case file
  character(len=:), allocatable :: &
    line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
  integer, parameter :: &
    subStepFactor = 2                                                                               !< for each substep, divide the last time increment by 2.0
  real(pREAL) :: &
    time = 0.0_pREAL, &                                                                             !< elapsed time
    time0 = 0.0_pREAL, &                                                                            !< begin of interval
    timeinc = 0.0_pREAL, &                                                                          !< current time interval
    timeIncOld = 0.0_pREAL, &                                                                       !< previous time interval
    remainingLoadCaseTime = 0.0_pREAL                                                               !< remaining time of current load case
  logical :: &
    guess, &                                                                                        !< guess along former trajectory
    stagIterate
  integer :: &
    l, &
    i, &
    errorID, &
    cutBackLevel = 0, &                                                                             !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
    stepFraction = 0, &                                                                             !< fraction of current time interval
    currentLoadcase = 0, &                                                                          !< current load case
    currentFace = 0, &
    inc, &                                                                                          !< current increment in current load case
    totalIncsCounter = 0, &                                                                         !< total # of increments
    statUnit = 0, &                                                                                 !< file unit for statistics output
    stagIter, &
    component
  type(tDict), pointer :: &
    num_mesh
  character(len=pSTRLEN), dimension(:), allocatable :: fileContent
  character(len=pSTRLEN) :: &
    incInfo, &
    loadcase_string
  integer :: &
    stagItMax, &                                                                                    !< max number of field level staggered iterations
    maxCutBack                                                                                      !< max number of cutbacks

  type(tLoadCase), allocatable, dimension(:) :: loadCases                                           !< array of all load cases
  type(tSolutionState), allocatable, dimension(:) :: solres
  PetscInt :: faceSet, currentFaceSet, dimPlex
  PetscErrorCode :: err_PETSc
  integer(kind(COMPONENT_UNDEFINED_ID)) :: ID
  external :: &
    quit

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
  call materialpoint_initAll()
  print'(/,1x,a)', '<<<+-  DAMASK_mesh init  -+>>>'; flush(IO_STDOUT)

!---------------------------------------------------------------------
! reading field information from numerics file and do sanity checks
  num_mesh => config_numerics%get_dict('mesh', defaultVal=emptyDict)
  stagItMax  = num_mesh%get_asInt('maxStaggeredIter',defaultVal=10)
  maxCutBack = num_mesh%get_asInt('maxCutBack',defaultVal=3)

  if (stagItMax < 0)  call IO_error(301,ext_msg='maxStaggeredIter')
  if (maxCutBack < 0) call IO_error(301,ext_msg='maxCutBack')

! reading basic information from load case file and allocate data structure containing load cases
  call DMGetDimension(geomMesh,dimPlex,err_PETSc)                                                   !< dimension of mesh (2D or 3D)
  CHKERRA(err_PETSc)
  allocate(solres(1))

!--------------------------------------------------------------------------------------------------
! reading basic information from load case file and allocate data structure containing load cases
  fileContent = IO_readlines(trim(CLI_loadFile))
  do l = 1, size(fileContent)
    line = fileContent(l)
    if (IO_isBlank(line)) cycle                                                                     ! skip empty lines

    chunkPos = IO_strPos(line)
    do i = 1, chunkPos(1)                                                                           ! reading compulsory parameters for loadcase
      select case (IO_strValue(line,chunkPos,i))
        case('$Loadcase')
          N_def = N_def + 1
      end select
    end do                                                                                           ! count all identifiers to allocate memory and do sanity check
  end do

  if (N_def < 1) call IO_error(error_ID = 837)
  allocate(loadCases(N_def))

  do i = 1, size(loadCases)
    allocate(loadCases(i)%fieldBC(1))
    loadCases(i)%fieldBC(1)%ID = FIELD_MECH_ID
  end do

  do i = 1, size(loadCases)
    loadCases(i)%fieldBC(1)%nComponents = dimPlex                                                   !< X, Y (, Z) displacements
    allocate(loadCases(i)%fieldBC(1)%componentBC(loadCases(i)%fieldBC(1)%nComponents))
    do component = 1, loadCases(i)%fieldBC(1)%nComponents
      select case (component)
        case (1)
          loadCases(i)%fieldBC(1)%componentBC(component)%ID = COMPONENT_MECH_X_ID
        case (2)
          loadCases(i)%fieldBC(1)%componentBC(component)%ID = COMPONENT_MECH_Y_ID
        case (3)
          loadCases(i)%fieldBC(1)%componentBC(component)%ID = COMPONENT_MECH_Z_ID
      end select
    end do
    do component = 1, loadCases(i)%fieldBC(1)%nComponents
      allocate(loadCases(i)%fieldBC(1)%componentBC(component)%Value(mesh_Nboundaries), source = 0.0_pREAL)
      allocate(loadCases(i)%fieldBC(1)%componentBC(component)%Mask (mesh_Nboundaries), source = .false.)
    end do
  end do

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
  do l = 1, size(fileContent)
    line = fileContent(l)
    if (IO_isBlank(line)) cycle                                                                     ! skip empty lines

    chunkPos = IO_strPos(line)
    do i = 1, chunkPos(1)
      select case (IO_strValue(line,chunkPos,i))
!--------------------------------------------------------------------------------------------------
! loadcase information
        case('$Loadcase')
          currentLoadCase = IO_intValue(line,chunkPos,i+1)
        case('Face')
          currentFace = IO_intValue(line,chunkPos,i+1)
          currentFaceSet = -1
          do faceSet = 1, mesh_Nboundaries
            if (mesh_boundaries(faceSet) == currentFace) currentFaceSet = faceSet
          end do
          if (currentFaceSet < 0) call IO_error(error_ID = 837, ext_msg = 'invalid BC')
        case('t')
          loadCases(currentLoadCase)%time = IO_realValue(line,chunkPos,i+1)
        case('N')
          loadCases(currentLoadCase)%incs = IO_intValue(line,chunkPos,i+1)
        case('f_out')
          loadCases(currentLoadCase)%outputfrequency = IO_intValue(line,chunkPos,i+1)
        case('estimate_rate')
          loadCases(currentLoadCase)%followFormerTrajectory = .false.                                ! do not continue to predict deformation along former trajectory

!--------------------------------------------------------------------------------------------------
! boundary condition information
        case('X','Y','Z')
          select case(IO_strValue(line,chunkPos,i))
            case('X')
              ID = COMPONENT_MECH_X_ID
            case('Y')
              ID = COMPONENT_MECH_Y_ID
            case('Z')
              ID = COMPONENT_MECH_Z_ID
           end select

           do component = 1, loadcases(currentLoadCase)%fieldBC(1)%nComponents
             if (loadCases(currentLoadCase)%fieldBC(1)%componentBC(component)%ID == ID) then
               loadCases(currentLoadCase)%fieldBC(1)%componentBC(component)%Mask (currentFaceSet) = &
                   .true.
               loadCases(currentLoadCase)%fieldBC(1)%componentBC(component)%Value(currentFaceSet) = &
                   IO_realValue(line,chunkPos,i+1)
             end if
           end do
      end select
    end do
  end do

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
  loadCases(1)%followFormerTrajectory = .false.                                                      ! cannot guess along trajectory for first inc of first currentLoadCase
  errorID = 0
  checkLoadcases: do currentLoadCase = 1, size(loadCases)
    write (loadcase_string, '(i0)' ) currentLoadCase
    print'(/,1x,a,1x,i0)', 'load case:', currentLoadCase
    if (.not. loadCases(currentLoadCase)%followFormerTrajectory) &
      print'(2x,a)', 'drop guessing along trajectory'
    print'(2x,a)', 'Field '//trim(FIELD_MECH_label)

    do faceSet = 1, mesh_Nboundaries
       do component = 1, loadCases(currentLoadCase)%fieldBC(1)%nComponents
         if (loadCases(currentLoadCase)%fieldBC(1)%componentBC(component)%Mask(faceSet)) &
           print'(a,i2,a,i2,a,f12.7)', &
           '    Face ', mesh_boundaries(faceSet), &
           ' Component ', component, &
           ' Value ', loadCases(currentLoadCase)%fieldBC(1)%componentBC(component)%Value(faceSet)
      end do
    end do
    print'(2x,a,f12.6)', 'time:       ', loadCases(currentLoadCase)%time
    if (loadCases(currentLoadCase)%incs < 1)             errorID = 835                            ! non-positive incs count
    print'(2x,a,i5)',    'increments: ', loadCases(currentLoadCase)%incs
    if (loadCases(currentLoadCase)%outputfrequency < 1)  errorID = 836                            ! non-positive result frequency
    print'(2x,a,i5)',    'output frequency: ', &
               loadCases(currentLoadCase)%outputfrequency
    if (errorID > 0) call IO_error(error_ID = errorID, ext_msg = loadcase_string)                 ! exit with error message
  end do checkLoadcases

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call FEM_Utilities_init()
  call FEM_mechanical_init(loadCases(1)%fieldBC(1))
  call config_numerics_deallocate()

  if (worldrank == 0) then
    open(newunit=statUnit,file=trim(getSolverJobName())//'.sta',form='FORMATTED',status='REPLACE')
    write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                  ! statistics file
  end if

  print'(/,1x,a)', '... saving initial configuration ..........................................'
  flush(IO_STDOUT)
  call materialpoint_result(0,0.0_pREAL)

  loadCaseLooping: do currentLoadCase = 1, size(loadCases)
    time0 = time                                                                                    ! load case start time
    guess = loadCases(currentLoadCase)%followFormerTrajectory                                       ! change of load case? homogeneous guess for the first inc

    incLooping: do inc = 1, loadCases(currentLoadCase)%incs
      totalIncsCounter = totalIncsCounter + 1

!--------------------------------------------------------------------------------------------------
! forwarding time
      timeIncOld = timeinc                                                                          ! last timeinc that brought former inc to an end
      timeinc = loadCases(currentLoadCase)%time/real(loadCases(currentLoadCase)%incs,pREAL)
      timeinc = timeinc * real(subStepFactor,pREAL)**real(-cutBackLevel,pREAL)                      ! depending on cut back level, decrease time step
      stepFraction = 0                                                                              ! fraction scaled by stepFactor**cutLevel

      subStepLooping: do while (stepFraction < subStepFactor**cutBackLevel)
        remainingLoadCaseTime = loadCases(currentLoadCase)%time+time0 - time
        time = time + timeinc                                                                       ! forward target time
        stepFraction = stepFraction + 1                                                             ! count step

!--------------------------------------------------------------------------------------------------
! report begin of new step
        print'(/,1x,a)', '###########################################################################'
        print'(1x,a,es12.5,6(a,i0))',&
                'Time', time, &
                's: Increment ', inc, '/', loadCases(currentLoadCase)%incs,&
                '-', stepFraction, '/', subStepFactor**cutBackLevel,&
                ' of load case ', currentLoadCase,'/',size(loadCases)
        write(incInfo,'(4(a,i0))') &
               'Increment ',totalIncsCounter,'/',sum(loadCases%incs),&
               '-',stepFraction, '/', subStepFactor**cutBackLevel
        flush(IO_STDOUT)

        call FEM_mechanical_forward(guess,timeinc,timeIncOld,loadCases(currentLoadCase)%fieldBC(1))

!--------------------------------------------------------------------------------------------------
! solve fields
        stagIter = 0
        stagIterate = .true.
        do while (stagIterate)
          solres(1) = FEM_mechanical_solution(incInfo,timeinc,timeIncOld,loadCases(currentLoadCase)%fieldBC(1))
          if (.not. solres(1)%converged) exit

          stagIter = stagIter + 1
          stagIterate =            stagIter < stagItMax &
                       .and.       all(solres(:)%converged) &
                       .and. .not. all(solres(:)%stagConverged)                                     ! stationary with respect to staggered iteration
        end do

! check solution
        cutBack = .False.
        if (.not. all(solres(:)%converged .and. solres(:)%stagConverged)) then                       ! no solution found
          if (cutBackLevel < maxCutBack) then                                                       ! do cut back
            cutBack = .True.
            stepFraction = (stepFraction - 1) * subStepFactor                                       ! adjust to new denominator
            cutBackLevel = cutBackLevel + 1
            time    = time - timeinc                                                                ! rewind time
            timeinc = timeinc/2.0_pREAL
            print'(/,1x,a)', 'cutting back'
          else                                                                                      ! default behavior, exit if spectral solver does not converge
            if (worldrank == 0) close(statUnit)
            call IO_error(950)
          end if
        else
          guess = .true.                                                                            ! start guessing after first converged (sub)inc
          timeIncOld = timeinc
        end if
        if (.not. cutBack .and. worldrank == 0) then
          write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                            solres%converged, solres%iterationsNeeded                               ! write statistics about accepted solution
          flush(statUnit)
        end if
      end do subStepLooping

      cutBackLevel = max(0, cutBackLevel - 1)                                                       ! try half number of subincs next inc

      if (all(solres(:)%converged)) then
        print'(/,1x,a,1x,i0,1x,a)', 'increment', totalIncsCounter, 'converged'
      else
        print'(/,1x,a,1x,i0,1x,a)', 'increment', totalIncsCounter, 'NOT converged'
      end if; flush(IO_STDOUT)

      if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0) then                            ! at output frequency
        print'(/,1x,a)', '... saving results ........................................................'
        call FEM_mechanical_updateCoords()
        call materialpoint_result(totalIncsCounter,time)
      end if


    end do incLooping

  end do loadCaseLooping


!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  print'(/,1x,a)', '###########################################################################'
  if (worldrank == 0) close(statUnit)

  call quit(0)                                                                                      ! no complains ;)

end program DAMASK_mesh
