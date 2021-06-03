!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the various spectral solvers
!> @details doing cutbacking, forwarding in case of restart, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_grid
#include <petsc/finclude/petscsys.h>
  use PETScsys
  use prec
  use parallelization
  use DAMASK_interface
  use IO
  use config
  use math
  use CPFEM2
  use material
  use spectral_utilities
  use grid_mechanical_spectral_basic
  use grid_mechanical_spectral_polarisation
  use grid_mechanical_FEM
  use grid_damage_spectral
  use grid_thermal_spectral
  use results

  implicit none

  type :: tLoadCase
    type(rotation)           :: rot                                                                 !< rotation of BC
    type(tBoundaryCondition) :: stress, &                                                           !< stress BC
                                deformation                                                         !< deformation BC (dot_F, F, or L)
    real(pReal) ::              t, &                                                                !< length of increment
                                r                                                                   !< ratio of geometric progression
    integer ::                  N, &                                                                !< number of increments
                                f_out, &                                                            !< frequency of result writes
                                f_restart                                                           !< frequency of restart writes
    logical ::                  estimate_rate                                                       !< follow trajectory of former loadcase
  end type tLoadCase

  integer(kind(FIELD_UNDEFINED_ID)), allocatable :: ID(:)

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
  real(pReal), dimension(3,3), parameter :: &
    ones  = 1.0_pReal, &
    zeros = 0.0_pReal
  integer, parameter :: &
    subStepFactor = 2                                                                               !< for each substep, divide the last time increment by 2.0
  real(pReal) :: &
    T_0 = 300.0_pReal, &
    time = 0.0_pReal, &                                                                             !< elapsed time
    time0 = 0.0_pReal, &                                                                            !< begin of interval
    timeinc = 1.0_pReal, &                                                                          !< current time interval
    timeIncOld = 0.0_pReal, &                                                                       !< previous time interval
    remainingLoadCaseTime = 0.0_pReal                                                               !< remaining time of current load case
  logical :: &
    guess, &                                                                                        !< guess along former trajectory
    stagIterate, &
    cutBack = .false.,&
    signal
  integer :: &
    i, j, m, field, &
    errorID = 0, &
    ierr,&
    cutBackLevel = 0, &                                                                             !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
    stepFraction = 0, &                                                                             !< fraction of current time interval
    l = 0, &                                                                                        !< current load case
    inc, &                                                                                          !< current increment in current load case
    totalIncsCounter = 0, &                                                                         !< total # of increments
    statUnit = 0, &                                                                                 !< file unit for statistics output
    stagIter, &
    nActiveFields = 0, &
    maxCutBack, &                                                                                   !< max number of cut backs
    stagItMax                                                                                       !< max number of field level staggered iterations
  character(len=pStringLen) :: &
    incInfo

  type(tLoadCase), allocatable, dimension(:) :: loadCases                                           !< array of all load cases
  type(tSolutionState), allocatable, dimension(:) :: solres
  procedure(grid_mechanical_spectral_basic_init), pointer :: &
    mechanical_init
  procedure(grid_mechanical_spectral_basic_forward), pointer :: &
    mechanical_forward
  procedure(grid_mechanical_spectral_basic_solution), pointer :: &
    mechanical_solution
  procedure(grid_mechanical_spectral_basic_updateCoords), pointer :: &
    mechanical_updateCoords
  procedure(grid_mechanical_spectral_basic_restartWrite), pointer :: &
    mechanical_restartWrite

  external :: &
    quit
  class (tNode), pointer :: &
    num_grid, &
    config_load, &
    load_steps, &
    load_step, &
    solver, &
    initial_conditions, &
    thermal, &
    step_bc, &
    step_mech, &
    step_discretization

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)

  call CPFEM_initAll
  print'(/,a)',   ' <<<+-  DAMASK_grid init  -+>>>'; flush(IO_STDOUT)

  print*, 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print*, 'https://doi.org/10.1007/978-981-10-6855-3_80'


!-------------------------------------------------------------------------------------------------
! reading field paramters from numerics file and do sanity checks
  num_grid => config_numerics%get('grid', defaultVal=emptyDict)
  stagItMax  = num_grid%get_asInt('maxStaggeredIter',defaultVal=10)
  maxCutBack = num_grid%get_asInt('maxCutBack',defaultVal=3)

  if (stagItMax < 0)    call IO_error(301,ext_msg='maxStaggeredIter')
  if (maxCutBack < 0)   call IO_error(301,ext_msg='maxCutBack')

  config_load => YAML_parse_file(trim(interface_loadFile))
  solver => config_load%get('solver')

!--------------------------------------------------------------------------------------------------
! assign mechanics solver depending on selected type

  nActiveFields = 1
  select case (solver%get_asString('mechanical'))
    case ('spectral_basic')
      mechanical_init         => grid_mechanical_spectral_basic_init
      mechanical_forward      => grid_mechanical_spectral_basic_forward
      mechanical_solution     => grid_mechanical_spectral_basic_solution
      mechanical_updateCoords => grid_mechanical_spectral_basic_updateCoords
      mechanical_restartWrite => grid_mechanical_spectral_basic_restartWrite

    case ('spectral_polarization')
      mechanical_init         => grid_mechanical_spectral_polarisation_init
      mechanical_forward      => grid_mechanical_spectral_polarisation_forward
      mechanical_solution     => grid_mechanical_spectral_polarisation_solution
      mechanical_updateCoords => grid_mechanical_spectral_polarisation_updateCoords
      mechanical_restartWrite => grid_mechanical_spectral_polarisation_restartWrite

    case ('FEM')
      mechanical_init         => grid_mechanical_FEM_init
      mechanical_forward      => grid_mechanical_FEM_forward
      mechanical_solution     => grid_mechanical_FEM_solution
      mechanical_updateCoords => grid_mechanical_FEM_updateCoords
      mechanical_restartWrite => grid_mechanical_FEM_restartWrite

    case default
      call IO_error(error_ID = 891, ext_msg = trim(solver%get_asString('mechanical')))

  end select

!--------------------------------------------------------------------------------------------------
! initialize field solver information
  if (solver%get_asString('thermal',defaultVal = 'n/a') == 'spectral') nActiveFields = nActiveFields + 1
  if (solver%get_asString('damage', defaultVal = 'n/a') == 'spectral') nActiveFields = nActiveFields + 1

  allocate(solres(nActiveFields))
  allocate(    ID(nActiveFields))

  field = 1
  ID(field) = FIELD_MECH_ID                                                                         ! mechanical active by default
  thermalActive: if (solver%get_asString('thermal',defaultVal = 'n/a') == 'spectral') then
    field = field + 1
    ID(field) = FIELD_THERMAL_ID
  endif thermalActive
  damageActive: if (solver%get_asString('damage',defaultVal = 'n/a') == 'spectral') then
    field = field + 1
    ID(field) = FIELD_DAMAGE_ID
  endif damageActive


!--------------------------------------------------------------------------------------------------
  load_steps => config_load%get('loadstep')
  allocate(loadCases(load_steps%length))                                                            ! array of load cases

  do l = 1, load_steps%length

    load_step => load_steps%get(l)
    step_bc   => load_step%get('boundary_conditions')
    step_mech => step_bc%get('mechanical')
    loadCases(l)%stress%myType=''
    readMech: do m = 1, step_mech%length
      select case (step_mech%getKey(m))
        case ('L','dot_F','F')                                                                      ! assign values for the deformation BC matrix
          loadCases(l)%deformation%myType = step_mech%getKey(m)
          call getMaskedTensor(loadCases(l)%deformation%values,loadCases(l)%deformation%mask,step_mech%get(m))
        case ('dot_P','P')
          loadCases(l)%stress%myType = step_mech%getKey(m)
          call getMaskedTensor(loadCases(l)%stress%values,loadCases(l)%stress%mask,step_mech%get(m))
      end select
      call loadCases(l)%rot%fromAxisAngle(step_mech%get_as1dFloat('R',defaultVal = real([0.0,0.0,1.0,0.0],pReal)),degrees=.true.)
    enddo readMech
    if (.not. allocated(loadCases(l)%deformation%myType)) call IO_error(error_ID=837,ext_msg = 'L/dot_F/F missing')

    step_discretization => load_step%get('discretization')
    if (.not. step_discretization%contains('t')) call IO_error(error_ID=837,ext_msg = 't missing')
    if (.not. step_discretization%contains('N')) call IO_error(error_ID=837,ext_msg = 'N missing')
    loadCases(l)%t         = step_discretization%get_asFloat('t')
    loadCases(l)%N         = step_discretization%get_asInt  ('N')
    loadCases(l)%r         = step_discretization%get_asFloat('r',         defaultVal= 1.0_pReal)

    loadCases(l)%f_restart = load_step%get_asInt('f_restart', defaultVal=huge(0))
    loadCases(l)%f_out     = load_step%get_asInt('f_out',     defaultVal=1)
    loadCases(l)%estimate_rate = (load_step%get_asBool('estimate_rate',defaultVal=.true.) .and. l>1)

    reportAndCheck: if (worldrank == 0) then
      print'(/,a,i0)', ' load case: ', l
      print*, ' estimate_rate:', loadCases(l)%estimate_rate
      if (loadCases(l)%deformation%myType == 'L') then
        do j = 1, 3
          if (any(loadCases(l)%deformation%mask(j,1:3) .eqv. .true.) .and. &
              any(loadCases(l)%deformation%mask(j,1:3) .eqv. .false.)) errorID = 832                ! each row should be either fully or not at all defined
        enddo
      endif
      if (loadCases(l)%deformation%myType == 'F') then
        print*, ' F:'
      else
        print*, ' '//loadCases(l)%deformation%myType//' / 1/s:'
      endif
      do i = 1, 3; do j = 1, 3
        if (loadCases(l)%deformation%mask(i,j)) then
          write(IO_STDOUT,'(2x,f12.7)',advance='no') loadCases(l)%deformation%values(i,j)
        else
          write(IO_STDOUT,'(2x,12a)',advance='no') '     x      '
        endif
        enddo; write(IO_STDOUT,'(/)',advance='no')
      enddo
      if (any(loadCases(l)%stress%mask .eqv. loadCases(l)%deformation%mask)) errorID = 831
      if (any(loadCases(l)%stress%mask .and. transpose(loadCases(l)%stress%mask) .and. (math_I3<1))) &
        errorID = 838                                                                               ! no rotation is allowed by stress BC

      if (loadCases(l)%stress%myType == 'P')     print*, ' P / MPa:'
      if (loadCases(l)%stress%myType == 'dot_P') print*, ' dot_P / MPa/s:'

      if (loadCases(l)%stress%myType /= '') then
        do i = 1, 3; do j = 1, 3
          if (loadCases(l)%stress%mask(i,j)) then
            write(IO_STDOUT,'(2x,f12.4)',advance='no') loadCases(l)%stress%values(i,j)*1e-6_pReal
          else
            write(IO_STDOUT,'(2x,12a)',advance='no') '     x      '
          endif
          enddo; write(IO_STDOUT,'(/)',advance='no')
        enddo
      endif
      if (any(dNeq(loadCases(l)%rot%asMatrix(), math_I3))) &
        write(IO_STDOUT,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'R:',&
                 transpose(loadCases(l)%rot%asMatrix())

      if (loadCases(l)%r <= 0.0)       errorID = 833
      if (loadCases(l)%t < 0.0_pReal)  errorID = 834
      if (loadCases(l)%N < 1)          errorID = 835
      if (loadCases(l)%f_out < 1)      errorID = 836
      if (loadCases(l)%f_restart < 1)  errorID = 839

      if (dEq(loadCases(l)%r,1.0_pReal,1.e-9_pReal)) then
        print'(a)', '  r: 1 (constant step width)'
      else
        print'(a,f0.3)', '  r: ', loadCases(l)%r
      endif
      print'(a,f0.3)',   '  t: ', loadCases(l)%t
      print'(a,i0)',     '  N: ', loadCases(l)%N
      print'(a,i0)',     '  f_out: ', loadCases(l)%f_out
      if (loadCases(l)%f_restart < huge(0)) &
        print'(a,i0)',   '  f_restart: ', loadCases(l)%f_restart

      if (errorID > 0) call IO_error(error_ID = errorID, el = l)

    endif reportAndCheck
  enddo

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call spectral_Utilities_init
  do field = 1, nActiveFields
    select case (ID(field))
      case(FIELD_MECH_ID)
        call mechanical_init

      case(FIELD_THERMAL_ID)
        initial_conditions => config_load%get('initial_conditions',defaultVal=emptyDict)
        thermal            => initial_conditions%get('thermal',defaultVal=emptyDict)
        call grid_thermal_spectral_init(thermal%get_asFloat('T',defaultVal = T_0))

      case(FIELD_DAMAGE_ID)
        call grid_damage_spectral_init

    end select
  enddo

!--------------------------------------------------------------------------------------------------
! write header of output file
  if (worldrank == 0) then
    writeHeader: if (interface_restartInc < 1) then
      open(newunit=statUnit,file=trim(getSolverJobName())//'.sta',form='FORMATTED',status='REPLACE')
      write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                ! statistics file
    else writeHeader
      open(newunit=statUnit,file=trim(getSolverJobName())//&
                                  '.sta',form='FORMATTED', position='APPEND', status='OLD')
    endif writeHeader
  endif

  writeUndeformed: if (interface_restartInc < 1) then
    print'(/,a)', ' ... writing initial configuration to file ........................'
    flush(IO_STDOUT)
    call CPFEM_results(0,0.0_pReal)
  endif writeUndeformed

  loadCaseLooping: do l = 1, size(loadCases)
    time0 = time                                                                                    ! load case start time
    guess = loadCases(l)%estimate_rate                                                              ! change of load case? homogeneous guess for the first inc

    incLooping: do inc = 1, loadCases(l)%N
      totalIncsCounter = totalIncsCounter + 1

!--------------------------------------------------------------------------------------------------
! forwarding time
      timeIncOld = timeinc                                                                          ! last timeinc that brought former inc to an end
      if (dEq(loadCases(l)%r,1.0_pReal,1.e-9_pReal)) then                                           ! linear scale
        timeinc = loadCases(l)%t/real(loadCases(l)%N,pReal)
      else
        timeinc = loadCases(l)%t * (loadCases(l)%r**(inc-1)-loadCases(l)%r**inc) &
                                 / (1.0_pReal-loadCases(l)%r**loadCases(l)%N)
      endif
      timeinc = timeinc * real(subStepFactor,pReal)**real(-cutBackLevel,pReal)                      ! depending on cut back level, decrease time step

      skipping: if (totalIncsCounter <= interface_restartInc) then                                  ! not yet at restart inc?
        time = time + timeinc                                                                       ! just advance time, skip already performed calculation
        guess = .true.                                                                              ! QUESTION:why forced guessing instead of inheriting loadcase preference
      else skipping
        stepFraction = 0                                                                            ! fraction scaled by stepFactor**cutLevel

        subStepLooping: do while (stepFraction < subStepFactor**cutBackLevel)
          remainingLoadCaseTime = loadCases(l)%t+time0 - time
          time = time + timeinc                                                                     ! forward target time
          stepFraction = stepFraction + 1                                                           ! count step

!--------------------------------------------------------------------------------------------------
! report begin of new step
          print'(/,a)', ' ###########################################################################'
          print'(1x,a,es12.5,6(a,i0))', &
                  'Time', time, &
                  's: Increment ', inc,'/',loadCases(l)%N,&
                  '-', stepFraction,'/',subStepFactor**cutBackLevel,&
                  ' of load case ', l,'/',size(loadCases)
          write(incInfo,'(4(a,i0))') &
                  'Increment ',totalIncsCounter,'/',sum(loadCases%N),&
                  '-', stepFraction,'/',subStepFactor**cutBackLevel
          flush(IO_STDOUT)

!--------------------------------------------------------------------------------------------------
! forward fields
          do field = 1, nActiveFields
            select case(ID(field))
              case(FIELD_MECH_ID)
                call mechanical_forward (&
                        cutBack,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                        deformation_BC     = loadCases(l)%deformation, &
                        stress_BC          = loadCases(l)%stress, &
                        rotation_BC        = loadCases(l)%rot)

              case(FIELD_THERMAL_ID); call grid_thermal_spectral_forward(cutBack)
              case(FIELD_DAMAGE_ID);  call grid_damage_spectral_forward(cutBack)
            end select
          enddo
          if(.not. cutBack) call CPFEM_forward

!--------------------------------------------------------------------------------------------------
! solve fields
          stagIter = 0
          stagIterate = .true.
          do while (stagIterate)
            do field = 1, nActiveFields
              select case(ID(field))
                case(FIELD_MECH_ID)
                  solres(field) = mechanical_solution(incInfo)
                case(FIELD_THERMAL_ID)
                  solres(field) = grid_thermal_spectral_solution(timeinc)
                case(FIELD_DAMAGE_ID)
                  solres(field) = grid_damage_spectral_solution(timeinc)
              end select

              if (.not. solres(field)%converged) exit                                               ! no solution found

            enddo
            stagIter = stagIter + 1
            stagIterate =            stagIter < stagItMax &
                         .and.       all(solres(:)%converged) &
                         .and. .not. all(solres(:)%stagConverged)                                   ! stationary with respect to staggered iteration
          enddo

!--------------------------------------------------------------------------------------------------
! check solution for either advance or retry

          if ( (all(solres(:)%converged .and. solres(:)%stagConverged)) &                           ! converged
               .and. .not. solres(1)%termIll) then                                                  ! and acceptable solution found
            call mechanical_updateCoords
            timeIncOld = timeinc
            cutBack = .false.
            guess = .true.                                                                          ! start guessing after first converged (sub)inc
            if (worldrank == 0) then
              write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                                solres(1)%converged, solres(1)%iterationsNeeded
              flush(statUnit)
            endif
          elseif (cutBackLevel < maxCutBack) then                                                   ! further cutbacking tolerated?
            cutBack = .true.
            stepFraction = (stepFraction - 1) * subStepFactor                                       ! adjust to new denominator
            cutBackLevel = cutBackLevel + 1
            time    = time - timeinc                                                                ! rewind time
            timeinc = timeinc/real(subStepFactor,pReal)                                             ! cut timestep
            print'(/,a)', ' cutting back '
          else                                                                                      ! no more options to continue
            if (worldrank == 0) close(statUnit)
            call IO_error(950)
          endif

        enddo subStepLooping

        cutBackLevel = max(0, cutBackLevel - 1)                                                     ! try half number of subincs next inc

        if (all(solres(:)%converged)) then
          print'(/,a,i0,a)', ' increment ', totalIncsCounter, ' converged'
        else
          print'(/,a,i0,a)', ' increment ', totalIncsCounter, ' NOT converged'
        endif; flush(IO_STDOUT)

        call MPI_Allreduce(interface_SIGUSR1,signal,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
        if (ierr /= 0) error stop 'MPI error'
        if (mod(inc,loadCases(l)%f_out) == 0 .or. signal) then
          print'(1/,a)', ' ... writing results to file ......................................'
          flush(IO_STDOUT)
          call CPFEM_results(totalIncsCounter,time)
        endif
        if(signal) call interface_setSIGUSR1(.false.)
        call MPI_Allreduce(interface_SIGUSR2,signal,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
        if (ierr /= 0) error stop 'MPI error'
        if (mod(inc,loadCases(l)%f_restart) == 0 .or. signal) then
          call mechanical_restartWrite
          call CPFEM_restartWrite
        endif
        if(signal) call interface_setSIGUSR2(.false.)
        call MPI_Allreduce(interface_SIGTERM,signal,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
        if (ierr /= 0) error stop 'MPI error'
        if (signal) exit loadCaseLooping
      endif skipping

    enddo incLooping

  enddo loadCaseLooping


!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  print'(/,a)', ' ###########################################################################'
  if (worldrank == 0) close(statUnit)

  call quit(0)                                                                                      ! no complains ;)


contains

subroutine getMaskedTensor(values,mask,tensor)

  real(pReal), intent(out), dimension(3,3) :: values
  logical,     intent(out), dimension(3,3) :: mask
  class (tNode), pointer :: tensor

  class (tNode), pointer :: row
  integer :: i,j


  values = 0.0
  if (tensor%length == 9) then ! temporary support for deprecated 1D tensor
    do i = 1,9
      mask((i-1)/3+1,mod(i-1,3)+1) = tensor%get_asString(i) /= 'x'
      if (mask((i-1)/3+1,mod(i-1,3)+1)) values((i-1)/3+1,mod(i-1,3)+1) = tensor%get_asFloat(i)
    enddo
  else
    do i = 1,3
      row => tensor%get(i)
      do j = 1,3
        mask(i,j) = row%get_asString(j) /= 'x'                      ! ToDo change to np.masked behavior
        if (mask(i,j)) values(i,j) = row%get_asFloat(j)
      enddo
    enddo
  endif

end subroutine

end program DAMASK_grid
