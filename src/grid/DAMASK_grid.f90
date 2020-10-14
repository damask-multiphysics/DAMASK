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
  use grid_mech_spectral_basic
  use grid_mech_spectral_polarisation
  use grid_mech_FEM
  use grid_damage_spectral
  use grid_thermal_spectral
  use results
  
  implicit none

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
  real(pReal), dimension(9) :: temp_valueVector = 0.0_pReal                                         !< temporarily from loadcase file when reading in tensors (initialize to 0.0)
  logical,     dimension(9) :: temp_maskVector  = .false.                                           !< temporarily from loadcase file when reading in tensors
  integer, allocatable, dimension(:) :: chunkPos
  integer :: &
    N_t   = 0, &                                                                                    !< # of time indicators found in load case file
    N_n   = 0, &                                                                                    !< # of increment specifiers found in load case file
    N_def = 0                                                                                       !< # of rate of deformation specifiers found in load case file
  character(len=:), allocatable :: &
    line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
  real(pReal), dimension(3,3), parameter :: &
    ones  = 1.0_pReal, &
    zeros = 0.0_pReal
  integer, parameter :: &
    subStepFactor = 2                                                                               !< for each substep, divide the last time increment by 2.0
  real(pReal) :: &
    time = 0.0_pReal, &                                                                             !< elapsed time
    time0 = 0.0_pReal, &                                                                            !< begin of interval
    timeinc = 1.0_pReal, &                                                                          !< current time interval
    timeIncOld = 0.0_pReal, &                                                                       !< previous time interval
    remainingLoadCaseTime = 0.0_pReal                                                               !< remaining time of current load case
  logical :: &
    guess, &                                                                                        !< guess along former trajectory
    stagIterate, &
    cutBack = .false.
  integer :: &
    i, j, m, field, &
    errorID = 0, &
    cutBackLevel = 0, &                                                                             !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
    stepFraction = 0                                                                                !< fraction of current time interval
  integer :: &
    currentLoadcase = 0, &                                                                          !< current load case
    inc, &                                                                                          !< current increment in current load case
    totalIncsCounter = 0, &                                                                         !< total # of increments
    statUnit = 0, &                                                                                 !< file unit for statistics output
    stagIter, &
    nActiveFields = 0
  character(len=pStringLen) :: &
    incInfo, &
    loadcase_string
  integer :: &
    maxCutBack, &                                                                                   !< max number of cut backs
    stagItMax                                                                                       !< max number of field level staggered iterations
  type(tLoadCase), allocatable, dimension(:) :: loadCases                                           !< array of all load cases
  type(tLoadCase) :: newLoadCase
  type(tSolutionState), allocatable, dimension(:) :: solres
  procedure(grid_mech_spectral_basic_init), pointer :: &
    mech_init
  procedure(grid_mech_spectral_basic_forward), pointer :: &
    mech_forward
  procedure(grid_mech_spectral_basic_solution), pointer :: &
    mech_solution
  procedure(grid_mech_spectral_basic_updateCoords), pointer :: &
    mech_updateCoords
  procedure(grid_mech_spectral_basic_restartWrite), pointer :: &
    mech_restartWrite
  
  external :: &
    quit
  class (tNode), pointer :: &
    num_grid, &
    debug_grid, &                                                                                    ! pointer to grid debug options
    config_load, &
    load_steps, &
    load_step, &
    step_mech, &
    step_discretization, &
    step_deformation, &
    step_stress

!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
 
  call CPFEM_initAll 
  print'(/,a)',   ' <<<+-  DAMASK_spectral init  -+>>>'; flush(IO_STDOUT)
  
  print*, 'Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print*, 'https://doi.org/10.1007/978-981-10-6855-3_80'

!--------------------------------------------------------------------------------------------------
! initialize field solver information
  nActiveFields = 1
  if (any(thermal_type  == THERMAL_conduction_ID  )) nActiveFields = nActiveFields + 1
  if (any(damage_type   == DAMAGE_nonlocal_ID     )) nActiveFields = nActiveFields + 1
  allocate(solres(nActiveFields))
  allocate(newLoadCase%ID(nActiveFields))

!-------------------------------------------------------------------------------------------------
! reading field paramters from numerics file and do sanity checks
  num_grid => config_numerics%get('grid', defaultVal=emptyDict)
  stagItMax  = num_grid%get_asInt('maxStaggeredIter',defaultVal=10)
  maxCutBack = num_grid%get_asInt('maxCutBack',defaultVal=3)

  if (stagItMax < 0)    call IO_error(301,ext_msg='maxStaggeredIter')
  if (maxCutBack < 0)   call IO_error(301,ext_msg='maxCutBack')

!--------------------------------------------------------------------------------------------------
! assign mechanics solver depending on selected type
 
  debug_grid => config_debug%get('grid',defaultVal=emptyList)
  select case (trim(num_grid%get_asString('solver', defaultVal = 'Basic'))) 
    case ('Basic')
      mech_init         => grid_mech_spectral_basic_init
      mech_forward      => grid_mech_spectral_basic_forward
      mech_solution     => grid_mech_spectral_basic_solution
      mech_updateCoords => grid_mech_spectral_basic_updateCoords
      mech_restartWrite => grid_mech_spectral_basic_restartWrite
  
    case ('Polarisation')
     if(debug_grid%contains('basic')) &
        call IO_warning(42, ext_msg='debug Divergence')
      mech_init         => grid_mech_spectral_polarisation_init
      mech_forward      => grid_mech_spectral_polarisation_forward
      mech_solution     => grid_mech_spectral_polarisation_solution
      mech_updateCoords => grid_mech_spectral_polarisation_updateCoords
      mech_restartWrite => grid_mech_spectral_polarisation_restartWrite
  
    case ('FEM')
     if(debug_grid%contains('basic')) &
        call IO_warning(42, ext_msg='debug Divergence')
      mech_init         => grid_mech_FEM_init
      mech_forward      => grid_mech_FEM_forward
      mech_solution     => grid_mech_FEM_solution
      mech_updateCoords => grid_mech_FEM_updateCoords
      mech_restartWrite => grid_mech_FEM_restartWrite
  
    case default
      call IO_error(error_ID = 891, ext_msg = trim(num_grid%get_asString('solver')))
  
  end select
  

!--------------------------------------------------------------------------------------------------
! reading information from load case file and to sanity checks 
  config_load => YAML_parse_file(trim(interface_loadFile)) 
 
  load_steps => config_load%get('step')
  allocate(loadCases(load_steps%length))                                                            ! array of load cases

  do currentLoadCase = 1, load_steps%length
    load_step => load_steps%get(currentLoadCase)
    step_mech => load_step%get('mech')
    step_discretization => load_step%get('discretization')
    do m = 1, step_mech%length
      select case (step_mech%getKey(m))
        case('L','Fdot','F')
          N_def = N_def + 1
      end select
    enddo
    if(step_discretization%contains('t')) N_t = N_t + 1
    if(step_discretization%contains('N')) N_n = N_n + 1
    if ((N_def /= N_n) .or. (N_n /= N_t) .or. N_n < 1) &                                            ! sanity check
      call IO_error(error_ID=837,el=currentLoadCase,ext_msg = trim(interface_loadFile))             ! error message for incomplete loadcase
  
    newLoadCase%stress%myType='p'
    field = 1
    newLoadCase%ID(field) = FIELD_MECH_ID                                                           ! mechanical active by default
    thermalActive: if (any(thermal_type  == THERMAL_conduction_ID)) then
      field = field + 1
      newLoadCase%ID(field) = FIELD_THERMAL_ID
    endif thermalActive
    damageActive: if (any(damage_type   == DAMAGE_nonlocal_ID)) then
      field = field + 1
      newLoadCase%ID(field) = FIELD_DAMAGE_ID
    endif damageActive
 
    readMech: do m = 1, step_mech%length
      select case (step_mech%getKey(m))
        case('Fdot','L','F')                                                                        ! assign values for the deformation BC matrix
          temp_valueVector = 0.0_pReal 
          if (step_mech%getKey(m) == 'Fdot') then                                                   ! in case of Fdot, set type to fdot
            newLoadCase%deformation%myType = 'fdot'
          else if (step_mech%getKey(m) == 'F') then
            newLoadCase%deformation%myType = 'f'
          else
            newLoadCase%deformation%myType = 'l'
          endif
          step_deformation => step_mech%get(m)
          do j = 1, 9
            temp_maskVector(j) = step_deformation%get_asString(j) /= 'x'                            ! true if not a x
            if (temp_maskVector(j)) temp_valueVector(j) = step_deformation%get_asFloat(j)           ! read value where applicable
          enddo
          newLoadCase%deformation%mask   = transpose(reshape(temp_maskVector,[ 3,3]))               ! mask in 3x3 notation
          newLoadCase%deformation%values = math_9to33(temp_valueVector)                             ! values in 3x3 notation
        case('P')
          temp_valueVector = 0.0_pReal
          step_stress => step_mech%get(m)
          do j = 1, 9
            temp_maskVector(j) = step_stress%get_asString(j) /= 'x'                                 ! true if not an asterisk
            if (temp_maskVector(j)) temp_valueVector(j) = step_stress%get_asFloat(j)                ! read value where applicable
          enddo
          newLoadCase%stress%mask   = transpose(reshape(temp_maskVector,[ 3,3]))
          newLoadCase%stress%values = math_9to33(temp_valueVector)
      end select
    enddo readMech

!--------------------------------------------------------------------------------------------------------
! Read discretization parameters
!--------------------------------------------------------------------------------------------------------
    newLoadCase%time             = step_discretization%get_asFloat('t')
    newLoadCase%incs             = step_discretization%get_asFloat('N')
    newLoadCase%logscale         = step_discretization%get_asBool ('log_timestep', defaultVal= .false.)
    newLoadCase%outputfrequency  = step_discretization%get_asInt  ('f_out',        defaultVal=1)
    newLoadCase%restartfrequency = step_discretization%get_asInt  ('f_restart',    defaultVal=huge(0))
    call newLoadCase%rot%fromAxisAngle(step_discretization%get_asFloats('R', &
                                                         defaultVal = real([0.0,0.0,1.0,0.0],pReal)),degrees=.true.)
            
    newLoadCase%followFormerTrajectory = .not. load_step%get_asBool('drop_guessing',defaultVal=.false.) ! do not continue to predict deformation along former trajectory
  
    newLoadCase%followFormerTrajectory = merge(.true.,.false.,currentLoadCase > 1)                  ! by default, guess from previous load case
  
    reportAndCheck: if (worldrank == 0) then
      write (loadcase_string, '(i0)' ) currentLoadCase
      print'(/,a,i0)', ' load case: ', currentLoadCase
      if (.not. newLoadCase%followFormerTrajectory) &
        print*, ' drop guessing along trajectory'
      if (newLoadCase%deformation%myType == 'l') then
        do j = 1, 3
          if (any(newLoadCase%deformation%mask(j,1:3) .eqv. .true.) .and. &
              any(newLoadCase%deformation%mask(j,1:3) .eqv. .false.)) errorID = 832                 ! each row should be either fully or not at all defined
        enddo
        print*, ' velocity gradient:'
      else if (newLoadCase%deformation%myType == 'f') then
        print*, ' deformation gradient at end of load case:'
      else
        print*, ' deformation gradient rate:'
      endif
      do i = 1, 3; do j = 1, 3
        if(newLoadCase%deformation%mask(i,j)) then
          write(IO_STDOUT,'(2x,f12.7)',advance='no') newLoadCase%deformation%values(i,j)
        else
          write(IO_STDOUT,'(2x,12a)',advance='no') '     *      '
          endif
        enddo; write(IO_STDOUT,'(/)',advance='no')
      enddo
      if (any(newLoadCase%stress%mask .eqv. newLoadCase%deformation%mask)) errorID = 831            ! exclusive or masking only
      if (any(newLoadCase%stress%mask .and. transpose(newLoadCase%stress%mask) .and. (math_I3<1))) &
        errorID = 838                                                                               ! no rotation is allowed by stress BC
      print*, ' stress / GPa:'
      do i = 1, 3; do j = 1, 3
        if(newLoadCase%stress%mask(i,j)) then
          write(IO_STDOUT,'(2x,f12.7)',advance='no') newLoadCase%stress%values(i,j)*1e-9_pReal
        else
          write(IO_STDOUT,'(2x,12a)',advance='no') '     *      '
        endif
        enddo; write(IO_STDOUT,'(/)',advance='no')
      enddo
      if (any(abs(matmul(newLoadCase%rot%asMatrix(), &
                         transpose(newLoadCase%rot%asMatrix()))-math_I3) > &
                 reshape(spread(tol_math_check,1,9),[ 3,3]))) errorID = 846                         ! given rotation matrix contains strain
      if (any(dNeq(newLoadCase%rot%asMatrix(), math_I3))) &
        write(IO_STDOUT,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'rotation of loadframe:',&
                 transpose(newLoadCase%rot%asMatrix())
      if (newLoadCase%time < 0.0_pReal) errorID = 834                                               ! negative time increment
      print'(a,f0.3)', ' time: ', newLoadCase%time
      if (newLoadCase%incs < 1)    errorID = 835                                                    ! non-positive incs count
      print'(a,i0)',  ' increments: ', newLoadCase%incs
      if (newLoadCase%outputfrequency < 1)  errorID = 836                                           ! non-positive result frequency
      print'(a,i0)',  ' output frequency: ', newLoadCase%outputfrequency
      if (newLoadCase%restartfrequency < 1)  errorID = 839                                          ! non-positive restart frequency
      if (newLoadCase%restartfrequency < huge(0)) &
        print'(a,i0)',  ' restart frequency: ', newLoadCase%restartfrequency
      if (errorID > 0) call IO_error(error_ID = errorID, ext_msg = loadcase_string)                 ! exit with error message
    endif reportAndCheck
    loadCases(currentLoadCase) = newLoadCase                                                        ! load case is ok, append it
  enddo

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call spectral_Utilities_init
  do field = 1, nActiveFields
    select case (loadCases(1)%ID(field))
      case(FIELD_MECH_ID)
        call mech_init
      
      case(FIELD_THERMAL_ID)
        call grid_thermal_spectral_init
  
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
      if (debug_grid%contains('basic')) print'(/,a)', ' header of statistics file written out'
      flush(IO_STDOUT)
    else writeHeader
      open(newunit=statUnit,file=trim(getSolverJobName())//&
                                  '.sta',form='FORMATTED', position='APPEND', status='OLD')
    endif writeHeader
  endif
  
  writeUndeformed: if (interface_restartInc < 1) then
    print'(/,a)', ' ... writing initial configuration to file ........................'
    call CPFEM_results(0,0.0_pReal)
  endif writeUndeformed
  
  loadCaseLooping: do currentLoadCase = 1, size(loadCases)
    time0 = time                                                                                    ! load case start time
    guess = loadCases(currentLoadCase)%followFormerTrajectory                                       ! change of load case? homogeneous guess for the first inc
  
    incLooping: do inc = 1, loadCases(currentLoadCase)%incs
      totalIncsCounter = totalIncsCounter + 1

!--------------------------------------------------------------------------------------------------
! forwarding time
      timeIncOld = timeinc                                                                          ! last timeinc that brought former inc to an end
      if (.not. loadCases(currentLoadCase)%logscale) then                                           ! linear scale
        timeinc = loadCases(currentLoadCase)%time/real(loadCases(currentLoadCase)%incs,pReal)
      else
        if (currentLoadCase == 1) then                                                              ! 1st load case of logarithmic scale
          timeinc = loadCases(1)%time*(2.0_pReal**real(max(inc-1,1)-loadCases(1)%incs ,pReal))      ! assume 1st inc is equal to 2nd
        else                                                                                        ! not-1st load case of logarithmic scale
          timeinc = time0 * &
               ( (1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real( inc         ,pReal)/&
                                                     real(loadCases(currentLoadCase)%incs ,pReal))&
                -(1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real( inc-1  ,pReal)/&
                                                     real(loadCases(currentLoadCase)%incs ,pReal)))
        endif
      endif
      timeinc = timeinc * real(subStepFactor,pReal)**real(-cutBackLevel,pReal)                      ! depending on cut back level, decrease time step
  
      skipping: if (totalIncsCounter <= interface_restartInc) then                                  ! not yet at restart inc?
        time = time + timeinc                                                                       ! just advance time, skip already performed calculation
        guess = .true.                                                                              ! QUESTION:why forced guessing instead of inheriting loadcase preference
      else skipping
        stepFraction = 0                                                                            ! fraction scaled by stepFactor**cutLevel
  
        subStepLooping: do while (stepFraction < subStepFactor**cutBackLevel)
          remainingLoadCaseTime = loadCases(currentLoadCase)%time+time0 - time
          time = time + timeinc                                                                     ! forward target time
          stepFraction = stepFraction + 1                                                           ! count step

!--------------------------------------------------------------------------------------------------
! report begin of new step
          print'(/,a)', ' ###########################################################################'
          print'(1x,a,es12.5,6(a,i0))', &
                  'Time', time, &
                  's: Increment ', inc,'/',loadCases(currentLoadCase)%incs,&
                  '-', stepFraction,'/',subStepFactor**cutBackLevel,&
                  ' of load case ', currentLoadCase,'/',size(loadCases)
          write(incInfo,'(4(a,i0))') &
                  'Increment ',totalIncsCounter,'/',sum(loadCases%incs),&
                  '-', stepFraction,'/',subStepFactor**cutBackLevel
          flush(IO_STDOUT)

!--------------------------------------------------------------------------------------------------
! forward fields
          do field = 1, nActiveFields
            select case(loadCases(currentLoadCase)%ID(field))
              case(FIELD_MECH_ID)
                call mech_forward (&
                        cutBack,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                        deformation_BC     = loadCases(currentLoadCase)%deformation, &
                        stress_BC          = loadCases(currentLoadCase)%stress, &
                        rotation_BC        = loadCases(currentLoadCase)%rot)
  
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
              select case(loadCases(currentLoadCase)%ID(field))
                case(FIELD_MECH_ID)
                  solres(field) = mech_solution(incInfo)
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
            call mech_updateCoords
            timeIncOld = timeinc
            cutBack = .false.
            guess = .true.                                                                          ! start guessing after first converged (sub)inc
            if (worldrank == 0) then
              write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                                solres%converged, solres%iterationsNeeded
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
            call IO_warning(850)
            if (worldrank == 0) close(statUnit)
            call quit(0)                                                                            ! quit
          endif
  
        enddo subStepLooping
  
        cutBackLevel = max(0, cutBackLevel - 1)                                                     ! try half number of subincs next inc
  
        if (all(solres(:)%converged)) then
          print'(/,a,i0,a)', ' increment ', totalIncsCounter, ' converged'
        else
          print'(/,a,i0,a)', ' increment ', totalIncsCounter, ' NOT converged'
        endif; flush(IO_STDOUT)
  
        if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0) then                          ! at output frequency
          print'(1/,a)', ' ... writing results to file ......................................'
          flush(IO_STDOUT)
          call CPFEM_results(totalIncsCounter,time)
        endif
        if (mod(inc,loadCases(currentLoadCase)%restartFrequency) == 0) then
          call mech_restartWrite
          call CPFEM_restartWrite
        endif
      endif skipping
  
     enddo incLooping
  
  enddo loadCaseLooping
 
 
!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  print'(/,a)', ' ###########################################################################'
  if (worldrank == 0) close(statUnit)
  
  call quit(0)                                                                                      ! no complains ;)
  
end program DAMASK_grid
