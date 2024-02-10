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
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  use prec
  use parallelization
  use signal
  use CLI
  use IO
  use config
  use math
  use materialpoint
  use material
  use spectral_utilities
  use grid_mech_utilities
  use grid_mechanical_spectral_basic
  use grid_mechanical_spectral_polarization
  use grid_mechanical_FEM
  use grid_damage_spectral
  use grid_thermal_spectral
  use result

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif

  type :: tLoadCase
    type(tRotation)          :: rot                                                                 !< rotation of BC
    type(tBoundaryCondition) :: stress, &                                                           !< stress BC
                                deformation                                                         !< deformation BC (dot_F, F, or L)
    real(pREAL) ::              t, &                                                                !< length of increment
                                r                                                                   !< ratio of geometric progression
    integer ::                  N, &                                                                !< number of increments
                                f_out, &                                                            !< frequency of result writes
                                f_restart                                                           !< frequency of restart writes
    logical ::                  estimate_rate                                                       !< follow trajectory of former loadcase
  end type tLoadCase

!--------------------------------------------------------------------------------------------------
! field labels information
  enum, bind(c); enumerator :: &
    FIELD_UNDEFINED_ID, &
    FIELD_MECH_ID, &
    FIELD_THERMAL_ID, &
    FIELD_DAMAGE_ID
  end enum

  integer(kind(FIELD_UNDEFINED_ID)), allocatable :: ID(:)

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
  integer, parameter :: &
    subStepFactor = 2                                                                               !< for each substep, divide the last time increment by 2.0
  real(pREAL) :: &
    t = 0.0_pREAL, &                                                                                !< elapsed time
    t_0 = 0.0_pREAL, &                                                                              !< begin of interval
    Delta_t = 1.0_pREAL, &                                                                          !< current time interval
    Delta_t_prev = 0.0_pREAL, &                                                                     !< previous time interval
    t_remaining = 0.0_pREAL                                                                         !< remaining time of current load case
  logical :: &
    guess, &                                                                                        !< guess along former trajectory
    stagIterate, &
    cutBack = .false.,&
    sig
  integer :: &
    i, j, field, &
    errorID = 0, &
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
  integer(MPI_INTEGER_KIND) :: err_MPI
  character(len=pSTRLEN) :: &
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
  type(tDict), pointer :: &
    load, &
    num_solver, &
    num_grid, &
    solver
  character(len=:), allocatable :: &
    fileContent, fname


!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)

  call materialpoint_initAll()
  print'(/,1x,a)', '<<<+-  DAMASK_grid init  -+>>>'; flush(IO_STDOUT)

  print'(/,1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'


!-------------------------------------------------------------------------------------------------
! read (and check) field parameters from numerics file

  num_solver => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_grid => num_solver%get_dict('grid',defaultVal=emptyDict)

  stagItMax  = num_grid%get_asInt('N_staggered_iter_max',defaultVal=10)
  maxCutBack = num_grid%get_asInt('N_cutback_max',defaultVal=3)

  if (stagItMax < 0)    call IO_error(301,ext_msg='N_staggered_iter_max')
  if (maxCutBack < 0)   call IO_error(301,ext_msg='N_cutback_max')


  if (worldrank == 0) then
    fileContent = IO_read(CLI_loadFile)
    fname = CLI_loadFile
    if (scan(fname,'/') /= 0) fname = fname(scan(fname,'/',.true.)+1:)
    call result_openJobFile(parallel=.false.)
    call result_addSetupFile(fileContent,fname,'load case definition (grid solver)')
    call result_closeJobFile()
  end if

  call parallelization_bcast_str(fileContent)
  load => YAML_str_asDict(fileContent)
  solver => load%get_dict('solver')

!--------------------------------------------------------------------------------------------------
! assign mechanics solver depending on selected type

  nActiveFields = 1
  select case (solver%get_asStr('mechanical'))
    case ('spectral_basic')
      mechanical_init         => grid_mechanical_spectral_basic_init
      mechanical_forward      => grid_mechanical_spectral_basic_forward
      mechanical_solution     => grid_mechanical_spectral_basic_solution
      mechanical_updateCoords => grid_mechanical_spectral_basic_updateCoords
      mechanical_restartWrite => grid_mechanical_spectral_basic_restartWrite

    case ('spectral_polarization')
      mechanical_init         => grid_mechanical_spectral_polarization_init
      mechanical_forward      => grid_mechanical_spectral_polarization_forward
      mechanical_solution     => grid_mechanical_spectral_polarization_solution
      mechanical_updateCoords => grid_mechanical_spectral_polarization_updateCoords
      mechanical_restartWrite => grid_mechanical_spectral_polarization_restartWrite

    case ('FEM')
      mechanical_init         => grid_mechanical_FEM_init
      mechanical_forward      => grid_mechanical_FEM_forward
      mechanical_solution     => grid_mechanical_FEM_solution
      mechanical_updateCoords => grid_mechanical_FEM_updateCoords
      mechanical_restartWrite => grid_mechanical_FEM_restartWrite

    case default
      call IO_error(error_ID = 891, ext_msg = trim(solver%get_asStr('mechanical')))

  end select

!--------------------------------------------------------------------------------------------------
! initialize field solver information
  if (solver%get_asStr('thermal',defaultVal = 'n/a') == 'spectral') nActiveFields = nActiveFields + 1
  if (solver%get_asStr('damage', defaultVal = 'n/a') == 'spectral') nActiveFields = nActiveFields + 1

  allocate(solres(nActiveFields))
  allocate(    ID(nActiveFields))

  field = 1
  ID(field) = FIELD_MECH_ID                                                                         ! mechanical active by default
  thermalActive: if (solver%get_asStr('thermal',defaultVal = 'n/a') == 'spectral') then
    field = field + 1
    ID(field) = FIELD_THERMAL_ID
  end if thermalActive
  damageActive: if (solver%get_asStr('damage',defaultVal = 'n/a') == 'spectral') then
    field = field + 1
    ID(field) = FIELD_DAMAGE_ID
  end if damageActive

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call spectral_utilities_init()

  do field = 2, nActiveFields
    select case (ID(field))

      case (FIELD_THERMAL_ID)
        call grid_thermal_spectral_init(num_grid)

      case (FIELD_DAMAGE_ID)
        call grid_damage_spectral_init(num_grid)

    end select
  end do

  call mechanical_init(num_grid)
  call config_numerics_deallocate()

!--------------------------------------------------------------------------------------------------
! write header of output file
  if (worldrank == 0) then
    writeHeader: if (CLI_restartInc < 1) then
      open(newunit=statUnit,file=trim(getSolverJobName())//'.sta',form='FORMATTED',status='REPLACE')
      write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded StagIterationsNeeded' ! statistics file
    else writeHeader
      open(newunit=statUnit,file=trim(getSolverJobName())//&
                                  '.sta',form='FORMATTED', position='APPEND', status='OLD')
    end if writeHeader
  end if

  writeUndeformed: if (CLI_restartInc < 1) then
    print'(/,1x,a)', '... saving initial configuration ..........................................'
    flush(IO_STDOUT)
    call materialpoint_result(0,0.0_pREAL)
  end if writeUndeformed

  loadCases = parseLoadSteps(load%get_list('loadstep'))

  loadCaseLooping: do l = 1, size(loadCases)
    t_0 = t                                                                                         ! load case start time
    guess = loadCases(l)%estimate_rate                                                              ! change of load case? homogeneous guess for the first inc

    incLooping: do inc = 1, loadCases(l)%N
      totalIncsCounter = totalIncsCounter + 1

!--------------------------------------------------------------------------------------------------
! forwarding time
      Delta_t_prev = Delta_t                                                                        ! last time intervall that brought former inc to an end
      if (dEq(loadCases(l)%r,1.0_pREAL,1.e-9_pREAL)) then                                           ! linear scale
        Delta_t = loadCases(l)%t/real(loadCases(l)%N,pREAL)
      else
        Delta_t = loadCases(l)%t * (loadCases(l)%r**(inc-1)-loadCases(l)%r**inc) &
                                 / (1.0_pREAL-loadCases(l)%r**loadCases(l)%N)
      end if
      Delta_t = Delta_t * real(subStepFactor,pREAL)**real(-cutBackLevel,pREAL)                      ! depending on cut back level, decrease time step

      skipping: if (totalIncsCounter <= CLI_restartInc) then                                        ! not yet at restart inc?
        t = t + Delta_t                                                                             ! just advance time, skip already performed calculation
        guess = .true.                                                                              ! QUESTION:why forced guessing instead of inheriting loadcase preference
      else skipping
        stepFraction = 0                                                                            ! fraction scaled by stepFactor**cutLevel

        subStepLooping: do while (stepFraction < subStepFactor**cutBackLevel)
          t_remaining = loadCases(l)%t + t_0 - t
          t = t + Delta_t                                                                           ! forward target time
          stepFraction = stepFraction + 1                                                           ! count step

!--------------------------------------------------------------------------------------------------
! report beginning of new step
          print'(/,1x,a)', '###########################################################################'
          print'(1x,a,1x,es12.5,6(a,i0))', &
                  'Time', t, &
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
                        cutBack,guess,Delta_t,Delta_t_prev,t_remaining, &
                        deformation_BC = loadCases(l)%deformation, &
                        stress_BC      = loadCases(l)%stress, &
                        rotation_BC    = loadCases(l)%rot)

              case(FIELD_THERMAL_ID); call grid_thermal_spectral_forward(cutBack)
              case(FIELD_DAMAGE_ID);  call grid_damage_spectral_forward(cutBack)
            end select
          end do
          if (.not. cutBack) call materialpoint_forward

!--------------------------------------------------------------------------------------------------
! solve fields
          stagIter = 1
          stagIterate = .true.
          do while (stagIterate)

            if (nActiveFields > 1) print'(/,1x,a,i0)', 'Staggered Iteration ',stagIter
            do field = 1, nActiveFields
              select case(ID(field))
                case(FIELD_MECH_ID)
                  solres(field) = mechanical_solution(incInfo)
                case(FIELD_THERMAL_ID)
                  solres(field) = grid_thermal_spectral_solution(Delta_t)
                case(FIELD_DAMAGE_ID)
                  solres(field) = grid_damage_spectral_solution(Delta_t)
              end select

              if (.not. solres(field)%converged) exit                                               ! no solution found

            end do
            stagIter = stagIter + 1
            stagIterate =            stagIter <= stagItMax &
                         .and.       all(solres(:)%converged) &
                         .and. .not. all(solres(:)%stagConverged)                                   ! stationary with respect to staggered iteration
          end do

!--------------------------------------------------------------------------------------------------
! check solution and either advance or retry with smaller timestep

          if ( (all(solres(:)%converged .and. solres(:)%stagConverged)) &                           ! converged
               .and. .not. solres(1)%termIll) then                                                  ! and acceptable solution found
            call mechanical_updateCoords()
            Delta_t_prev = Delta_t
            cutBack = .false.
            guess = .true.                                                                          ! start guessing after first converged (sub)inc
            if (worldrank == 0) then
              write(statUnit,*) totalIncsCounter, t, cutBackLevel, &
                                solres(1)%converged, solres(1)%iterationsNeeded, StagIter
              flush(statUnit)
            end if
          elseif (cutBackLevel < maxCutBack) then                                                   ! further cutbacking tolerated?
            cutBack = .true.
            stepFraction = (stepFraction - 1) * subStepFactor                                       ! adjust to new denominator
            cutBackLevel = cutBackLevel + 1
            t = t - Delta_t
            Delta_t = Delta_t/real(subStepFactor,pREAL)                                             ! cut timestep
            print'(/,1x,a)', 'cutting back '
          else                                                                                      ! no more options to continue
            if (worldrank == 0) close(statUnit)
            call IO_error(950)
          end if

        end do subStepLooping

        cutBackLevel = max(0, cutBackLevel - 1)                                                     ! try half number of subincs next inc

        if (all(solres(:)%converged)) then
          print'(/,1x,a,1x,i0,1x,a)', 'increment', totalIncsCounter, 'converged'
        else
          print'(/,1x,a,1x,i0,1x,a)', 'increment', totalIncsCounter, 'NOT converged'
        end if; flush(IO_STDOUT)

        call MPI_Allreduce(signal_SIGUSR1,sig,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,err_MPI)
        call parallelization_chkerr(err_MPI)
        if (mod(inc,loadCases(l)%f_out) == 0 .or. sig) then
          print'(/,1x,a)', '... saving results ........................................................'
          flush(IO_STDOUT)
          call materialpoint_result(totalIncsCounter,t)
        end if
        if (sig) call signal_setSIGUSR1(.false.)
        call MPI_Allreduce(signal_SIGUSR2,sig,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,err_MPI)
        call parallelization_chkerr(err_MPI)
        if (mod(inc,loadCases(l)%f_restart) == 0 .or. sig) then
          do field = 1, nActiveFields
            select case (ID(field))
              case(FIELD_MECH_ID)
                call mechanical_restartWrite()
              case(FIELD_THERMAL_ID)
                call grid_thermal_spectral_restartWrite()
              case(FIELD_DAMAGE_ID)
                call grid_damage_spectral_restartWrite()
            end select
          end do
          call materialpoint_restartWrite()
        end if
        if (sig) call signal_setSIGUSR2(.false.)
        call MPI_Allreduce(signal_SIGINT,sig,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,err_MPI)
        call parallelization_chkerr(err_MPI)
        if (sig) exit loadCaseLooping
      end if skipping

    end do incLooping

  end do loadCaseLooping


!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  print'(/,1x,a)', '###########################################################################'
  if (worldrank == 0) close(statUnit)

  call quit(0)                                                                                      ! no complains ;)


contains

subroutine getMaskedTensor(values,mask,tensor)

  real(pREAL), intent(out), dimension(3,3) :: values
  logical,     intent(out), dimension(3,3) :: mask
  type(tList), pointer :: tensor

  type(tList), pointer :: row
  integer :: i,j


  values = 0.0_pREAL
  do i = 1,3
    row => tensor%get_list(i)
    do j = 1,3
      mask(i,j) = row%get_asStr(j) == 'x'
      if (.not. mask(i,j)) values(i,j) = row%get_asReal(j)
    end do
  end do

end subroutine getMaskedTensor


function parseLoadsteps(load_steps) result(loadCases)

  type(tList), intent(in), target :: load_steps
  type(tLoadCase), allocatable, dimension(:) :: loadCases                                           !< array of all load cases

  integer :: l,m
  type(tDict), pointer :: &
    load_step, &
    step_bc, &
    step_mech, &
    step_discretization
#ifdef __INTEL_LLVM_COMPILER
  type(tList), pointer :: &
    tensor
#endif


  allocate(loadCases(load_steps%length))
  do l = 1, load_steps%length
    load_step => load_steps%get_dict(l)
    step_bc   => load_step%get_dict('boundary_conditions')
    step_mech => step_bc%get_dict('mechanical')
    loadCases(l)%stress%myType=''
    readMech: do m = 1, step_mech%length
      select case (step_mech%key(m))
        case ('L','dot_F','F')                                                                      ! assign values for the deformation BC matrix
          loadCases(l)%deformation%myType = step_mech%key(m)
#ifdef __INTEL_LLVM_COMPILER
          tensor => step_mech%get_list(m)
          call getMaskedTensor(loadCases(l)%deformation%values,loadCases(l)%deformation%mask,tensor)
#else
          call getMaskedTensor(loadCases(l)%deformation%values,loadCases(l)%deformation%mask,step_mech%get_list(m))
#endif
        case ('dot_P','P')
          loadCases(l)%stress%myType = step_mech%key(m)
#ifdef __INTEL_LLVM_COMPILER
          tensor => step_mech%get_list(m)
          call getMaskedTensor(loadCases(l)%stress%values,loadCases(l)%stress%mask,tensor)
#else
          call getMaskedTensor(loadCases(l)%stress%values,loadCases(l)%stress%mask,step_mech%get_list(m))
#endif
      end select
      call loadCases(l)%rot%fromAxisAngle(step_mech%get_as1dReal('R',defaultVal = real([0.0,0.0,1.0,0.0],pREAL)),degrees=.true.)
    end do readMech
    if (.not. allocated(loadCases(l)%deformation%myType)) call IO_error(error_ID=837,ext_msg = 'L/dot_F/F missing')

    step_discretization => load_step%get_dict('discretization')
    loadCases(l)%t = step_discretization%get_asReal('t')
    loadCases(l)%N = step_discretization%get_asInt  ('N')
    loadCases(l)%r = step_discretization%get_asReal('r',defaultVal= 1.0_pREAL)

    loadCases(l)%f_restart = load_step%get_asInt('f_restart', defaultVal=huge(0))
    if (load_step%get_asStr('f_out',defaultVal='n/a') == 'none') then
      loadCases(l)%f_out = huge(0)
    else
      loadCases(l)%f_out = load_step%get_asInt('f_out', defaultVal=1)
    end if
    loadCases(l)%estimate_rate = (load_step%get_asBool('estimate_rate',defaultVal=.true.) .and. l>1)

    reportAndCheck: if (worldrank == 0) then
      print'(/,1x,a,1x,i0)', 'load case:', l
      print'(2x,a,1x,l1)', 'estimate_rate:', loadCases(l)%estimate_rate
      if (loadCases(l)%deformation%myType == 'F') then
        print'(2x,a)', 'F:'
      else
        print'(2x,a)', loadCases(l)%deformation%myType//' / 1/s:'
      end if
      do i = 1, 3; do j = 1, 3
        if (loadCases(l)%deformation%mask(i,j)) then
          write(IO_STDOUT,'(2x,12a)',advance='no') '     x      '
        else
          write(IO_STDOUT,'(2x,f12.7)',advance='no') loadCases(l)%deformation%values(i,j)
        end if
        end do; write(IO_STDOUT,'(/)',advance='no')
      end do
      if (any(loadCases(l)%stress%mask .eqv. loadCases(l)%deformation%mask)) errorID = 831
      if (any(.not.(loadCases(l)%stress%mask .or. transpose(loadCases(l)%stress%mask)) .and. (math_I3<1))) &
        errorID = 838                                                                               ! no rotation is allowed by stress BC

      if (loadCases(l)%stress%myType == 'P')     print'(2x,a)', 'P / MPa:'
      if (loadCases(l)%stress%myType == 'dot_P') print'(2x,a)', 'dot_P / MPa/s:'

      if (loadCases(l)%stress%myType /= '') then
        do i = 1, 3; do j = 1, 3
          if (loadCases(l)%stress%mask(i,j)) then
            write(IO_STDOUT,'(2x,12a)',advance='no') '     x      '
          else
            write(IO_STDOUT,'(2x,f12.4)',advance='no') loadCases(l)%stress%values(i,j)*1e-6_pREAL
          end if
          end do; write(IO_STDOUT,'(/)',advance='no')
        end do
      end if
      if (any(dNeq(loadCases(l)%rot%asMatrix(), math_I3))) &
        write(IO_STDOUT,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'R:',&
                 transpose(loadCases(l)%rot%asMatrix())

      if (loadCases(l)%r <= 0.0_pREAL) errorID = 833
      if (loadCases(l)%t < 0.0_pREAL)  errorID = 834
      if (loadCases(l)%N < 1)          errorID = 835
      if (loadCases(l)%f_out < 1)      errorID = 836
      if (loadCases(l)%f_restart < 1)  errorID = 839

      if (dEq(loadCases(l)%r,1.0_pREAL,1.e-9_pREAL)) then
        print'(2x,a)', 'r: 1 (constant step width)'
      else
        print'(2x,a,1x,f0.3)', 'r:', loadCases(l)%r
      end if
      print'(2x,a,1x,f0.3)',   't:', loadCases(l)%t
      print'(2x,a,1x,i0)',     'N:', loadCases(l)%N
      if (loadCases(l)%f_out < huge(0)) &
        print'(2x,a,1x,i0)',   'f_out:', loadCases(l)%f_out
      if (loadCases(l)%f_restart < huge(0)) &
        print'(2x,a,1x,i0)',   'f_restart:', loadCases(l)%f_restart

      if (errorID > 0) call IO_error(errorID,label1='line',ID1=l)

    end if reportAndCheck
  end do
end function parseLoadsteps

end program DAMASK_grid
