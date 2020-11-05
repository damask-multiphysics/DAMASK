!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Chen Zhang, Michigan State University
!> @brief crystallite state integration functions and reporting of results
!--------------------------------------------------------------------------------------------------

module crystallite
  use prec
  use parallelization
  use IO
  use HDF5_utilities
  use DAMASK_interface
  use config
  use rotations
  use math
  use FEsolving
  use material
  use constitutive
  use discretization
  use lattice
  use results

  implicit none
  private

  real(pReal),               dimension(:,:,:),        allocatable, public :: &
    crystallite_dt                                                                                  !< requested time increment of each grain
  real(pReal),               dimension(:,:,:),        allocatable :: &
    crystallite_subdt, &                                                                            !< substepped time increment of each grain
    crystallite_subFrac, &                                                                          !< already calculated fraction of increment
    crystallite_subStep                                                                             !< size of next integration step
  type(rotation),            dimension(:,:,:),        allocatable :: &
    crystallite_orientation                                                                         !< current orientation
  real(pReal),               dimension(:,:,:,:,:),    allocatable :: &
    crystallite_F0, &                                                                               !< def grad at start of FE inc
    crystallite_subF,  &                                                                            !< def grad to be reached at end of crystallite inc
    crystallite_subF0, &                                                                            !< def grad at start of crystallite inc
    !
    crystallite_Fe, &                                                                               !< current "elastic" def grad (end of converged time step)
    !
    crystallite_Fp, &                                                                               !< current plastic def grad (end of converged time step)
    crystallite_Fp0, &                                                                              !< plastic def grad at start of FE inc
    crystallite_partitionedFp0,&                                                                    !< plastic def grad at start of homog inc
    crystallite_subFp0,&                                                                            !< plastic def grad at start of crystallite inc
    !
    crystallite_Fi, &                                                                               !< current intermediate def grad (end of converged time step)
    crystallite_Fi0, &                                                                              !< intermediate def grad at start of FE inc
    crystallite_partitionedFi0,&                                                                    !< intermediate def grad at start of homog inc
    crystallite_subFi0,&                                                                            !< intermediate def grad at start of crystallite inc
    !
    crystallite_Lp0, &                                                                              !< plastic velocitiy grad at start of FE inc
    crystallite_partitionedLp0, &                                                                   !< plastic velocity grad at start of homog inc
    !
    crystallite_Li, &                                                                               !< current intermediate velocitiy grad (end of converged time step)
    crystallite_Li0, &                                                                              !< intermediate velocitiy grad at start of FE inc
    crystallite_partitionedLi0, &                                                                   !< intermediate velocity grad at start of homog inc
    !
    crystallite_S0, &                                                                               !< 2nd Piola-Kirchhoff stress vector at start of FE inc
    crystallite_partitionedS0                                                                       !< 2nd Piola-Kirchhoff stress vector at start of homog inc
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public, protected :: &
    crystallite_P, &                                                                                !< 1st Piola-Kirchhoff stress per grain
    crystallite_Lp, &                                                                               !< current plastic velocitiy grad (end of converged time step)
    crystallite_S, &                                                                                !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
    crystallite_partitionedF0                                                                       !< def grad at start of homog inc
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
    crystallite_partitionedF                                                                        !< def grad to be reached at end of homog inc

  logical,                    dimension(:,:,:),         allocatable, public :: &
    crystallite_requested                                                                           !< used by upper level (homogenization) to request crystallite calculation
  logical,                    dimension(:,:,:),         allocatable :: &
    crystallite_converged                                                                           !< convergence flag

  type :: tOutput                                                                                   !< new requested output (per phase)
    character(len=pStringLen), allocatable, dimension(:) :: &
      label
  end type tOutput
  type(tOutput), allocatable, dimension(:) :: output_constituent

  type :: tNumerics
    integer :: &
      iJacoLpresiduum, &                                                                            !< frequency of Jacobian update of residuum in Lp
      nState, &                                                                                     !< state loop limit
      nStress                                                                                       !< stress loop limit
    real(pReal) :: &
      subStepMinCryst, &                                                                            !< minimum (relative) size of sub-step allowed during cutback
      subStepSizeCryst, &                                                                           !< size of first substep when cutback
      subStepSizeLp, &                                                                              !< size of first substep when cutback in Lp calculation
      subStepSizeLi, &                                                                              !< size of first substep when cutback in Li calculation
      stepIncreaseCryst, &                                                                          !< increase of next substep size when previous substep converged
      rtol_crystalliteState, &                                                                      !< relative tolerance in state loop
      rtol_crystalliteStress, &                                                                     !< relative tolerance in stress loop
      atol_crystalliteStress                                                                        !< absolute tolerance in stress loop
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

  type :: tDebugOptions
    logical :: &
      basic, &
      extensive, &
      selective
    integer :: &
      element, &
      ip, &
      grain
  end type tDebugOptions

  type(tDebugOptions) :: debugCrystallite

  procedure(integrateStateFPI), pointer :: integrateState

  public :: &
    crystallite_init, &
    crystallite_stress, &
    crystallite_stressTangent, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    crystallite_results, &
    crystallite_restartWrite, &
    crystallite_restartRead, &
    crystallite_forward, &
    crystallite_initializeRestorationPoints, &
    crystallite_windForward, &
    crystallite_restore

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init

  logical, dimension(discretization_nIPs,discretization_Nelems) :: devNull
  integer :: &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    cMax, &                                                                                         !< maximum number of  integration point components
    iMax, &                                                                                         !< maximum number of integration points
    eMax                                                                                            !< maximum number of elements


  class(tNode), pointer :: &
    num_crystallite, &
    debug_crystallite, &                                                                            ! pointer to debug options for crystallite
    phases, &
    phase, &
    generic_param

  print'(/,a)', ' <<<+-  crystallite init  -+>>>'

  debug_crystallite => config_debug%get('crystallite', defaultVal=emptyList)
  debugCrystallite%basic     = debug_crystallite%contains('basic')
  debugCrystallite%extensive = debug_crystallite%contains('extensive')
  debugCrystallite%selective = debug_crystallite%contains('selective')
  debugCrystallite%element   = config_debug%get_asInt('element', defaultVal=1)
  debugCrystallite%ip        = config_debug%get_asInt('integrationpoint', defaultVal=1)
  debugCrystallite%grain     = config_debug%get_asInt('grain', defaultVal=1)

  cMax = homogenization_maxNconstituents
  iMax = discretization_nIPs
  eMax = discretization_Nelems

  allocate(crystallite_partitionedF(3,3,cMax,iMax,eMax),source=0.0_pReal)

  allocate(crystallite_S0, &
           crystallite_F0, crystallite_Fi0,crystallite_Fp0, &
                           crystallite_Li0,crystallite_Lp0, &
           crystallite_partitionedS0, &
           crystallite_partitionedF0,crystallite_partitionedFp0,crystallite_partitionedFi0, &
                                   crystallite_partitionedLp0,crystallite_partitionedLi0, &
           crystallite_S,crystallite_P, &
           crystallite_Fe,crystallite_Fi,crystallite_Fp, &
                          crystallite_Li,crystallite_Lp, &
           crystallite_subF,crystallite_subF0, &
           crystallite_subFp0,crystallite_subFi0, &
           source = crystallite_partitionedF)

  allocate(crystallite_dt(cMax,iMax,eMax),source=0.0_pReal)
  allocate(crystallite_subdt,crystallite_subFrac,crystallite_subStep, &
           source = crystallite_dt)

  allocate(crystallite_orientation(cMax,iMax,eMax))

  allocate(crystallite_requested(cMax,iMax,eMax),             source=.false.)
  allocate(crystallite_converged(cMax,iMax,eMax),             source=.true.)

  num_crystallite => config_numerics%get('crystallite',defaultVal=emptyDict)

  num%subStepMinCryst        = num_crystallite%get_asFloat ('subStepMin',       defaultVal=1.0e-3_pReal)
  num%subStepSizeCryst       = num_crystallite%get_asFloat ('subStepSize',      defaultVal=0.25_pReal)
  num%stepIncreaseCryst      = num_crystallite%get_asFloat ('stepIncrease',     defaultVal=1.5_pReal)
  num%subStepSizeLp          = num_crystallite%get_asFloat ('subStepSizeLp',    defaultVal=0.5_pReal)
  num%subStepSizeLi          = num_crystallite%get_asFloat ('subStepSizeLi',    defaultVal=0.5_pReal)
  num%rtol_crystalliteState  = num_crystallite%get_asFloat ('rtol_State',       defaultVal=1.0e-6_pReal)
  num%rtol_crystalliteStress = num_crystallite%get_asFloat ('rtol_Stress',      defaultVal=1.0e-6_pReal)
  num%atol_crystalliteStress = num_crystallite%get_asFloat ('atol_Stress',      defaultVal=1.0e-8_pReal)
  num%iJacoLpresiduum        = num_crystallite%get_asInt   ('iJacoLpresiduum',  defaultVal=1)
  num%nState                 = num_crystallite%get_asInt   ('nState',           defaultVal=20)
  num%nStress                = num_crystallite%get_asInt   ('nStress',          defaultVal=40)

  if(num%subStepMinCryst   <= 0.0_pReal)      call IO_error(301,ext_msg='subStepMinCryst')
  if(num%subStepSizeCryst  <= 0.0_pReal)      call IO_error(301,ext_msg='subStepSizeCryst')
  if(num%stepIncreaseCryst <= 0.0_pReal)      call IO_error(301,ext_msg='stepIncreaseCryst')

  if(num%subStepSizeLp <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLp')
  if(num%subStepSizeLi <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLi')

  if(num%rtol_crystalliteState  <= 0.0_pReal) call IO_error(301,ext_msg='rtol_crystalliteState')
  if(num%rtol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='rtol_crystalliteStress')
  if(num%atol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='atol_crystalliteStress')

  if(num%iJacoLpresiduum < 1)                 call IO_error(301,ext_msg='iJacoLpresiduum')

  if(num%nState < 1)                          call IO_error(301,ext_msg='nState')
  if(num%nStress< 1)                          call IO_error(301,ext_msg='nStress')

  select case(num_crystallite%get_asString('integrator',defaultVal='FPI'))
    case('FPI')
      integrateState => integrateStateFPI
    case('Euler')
      integrateState => integrateStateEuler
    case('AdaptiveEuler')
      integrateState => integrateStateAdaptiveEuler
    case('RK4')
      integrateState => integrateStateRK4
    case('RKCK45')
      integrateState => integrateStateRKCK45
    case default
     call IO_error(301,ext_msg='integrator')
  end select

  phases => config_material%get('phase')

  allocate(output_constituent(phases%length))
  do c = 1, phases%length
    phase => phases%get(c)
    generic_param  => phase%get('generic',defaultVal = emptyDict)
#if defined(__GFORTRAN__)
    output_constituent(c)%label  = output_asStrings(generic_param)
#else
    output_constituent(c)%label  = generic_param%get_asStrings('output',defaultVal=emptyStringArray)
#endif
  enddo


!--------------------------------------------------------------------------------------------------
! initialize
 !$OMP PARALLEL DO PRIVATE(i,c)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1), FEsolving_execIP(2); do c = 1, homogenization_Nconstituents(material_homogenizationAt(e))
      crystallite_Fp0(1:3,1:3,c,i,e) = material_orientation0(c,i,e)%asMatrix()                      ! Fp reflects initial orientation (see 10.1016/j.actamat.2006.01.005)
      crystallite_Fp0(1:3,1:3,c,i,e) = crystallite_Fp0(1:3,1:3,c,i,e) &
                                     / math_det33(crystallite_Fp0(1:3,1:3,c,i,e))**(1.0_pReal/3.0_pReal)
      crystallite_Fi0(1:3,1:3,c,i,e) = constitutive_initialFi(c,i,e)
      crystallite_F0(1:3,1:3,c,i,e)  = math_I3
      crystallite_Fe(1:3,1:3,c,i,e)  = math_inv33(matmul(crystallite_Fi0(1:3,1:3,c,i,e), &
                                                         crystallite_Fp0(1:3,1:3,c,i,e)))           ! assuming that euler angles are given in internal strain free configuration
      crystallite_Fp(1:3,1:3,c,i,e)  = crystallite_Fp0(1:3,1:3,c,i,e)
      crystallite_Fi(1:3,1:3,c,i,e)  = crystallite_Fi0(1:3,1:3,c,i,e)
      crystallite_requested(c,i,e) = .true.
    enddo; enddo
  enddo
  !$OMP END PARALLEL DO


  crystallite_partitionedFp0 = crystallite_Fp0
  crystallite_partitionedFi0 = crystallite_Fi0
  crystallite_partitionedF0  = crystallite_F0
  crystallite_partitionedF   = crystallite_F0

  call crystallite_orientations()

  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2)
      do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
        call constitutive_dependentState(crystallite_partitionedF0(1:3,1:3,c,i,e), &
                                         crystallite_partitionedFp0(1:3,1:3,c,i,e), &
                                         c,i,e)                                                     ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  devNull = crystallite_stress()

#ifdef DEBUG
  if (debugCrystallite%basic) then
    print'(a42,1x,i10)', '    # of elements:                       ', eMax
    print'(a42,1x,i10)', '    # of integration points/element:     ', iMax
    print'(a42,1x,i10)', 'max # of constituents/integration point: ', cMax
    flush(IO_STDOUT)
  endif
#endif

end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P)
!--------------------------------------------------------------------------------------------------
function crystallite_stress()

  logical, dimension(discretization_nIPs,discretization_Nelems) :: crystallite_stress
  real(pReal) :: &
    formerSubStep
  integer :: &
    NiterationCrystallite, &                                                                        ! number of iterations in crystallite loop
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    s
  logical, dimension(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems) :: todo     !ToDo: need to set some values to false for different Ngrains
  real(pReal),               dimension(:,:,:,:,:),    allocatable :: &
    subLp0,&                                                                                        !< plastic velocity grad at start of crystallite inc
    subLi0                                                                                          !< intermediate velocity grad at start of crystallite inc


  todo = .false.

  subLp0 = crystallite_partitionedLp0
  subLi0 = crystallite_partitionedLi0



!--------------------------------------------------------------------------------------------------
! initialize to starting condition
  crystallite_subStep = 0.0_pReal
  !$OMP PARALLEL DO
  elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2); do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
      homogenizationRequestsCalculation: if (crystallite_requested(c,i,e)) then
        plasticState    (material_phaseAt(c,e))%subState0(      :,material_phaseMemberAt(c,i,e)) = &
        plasticState    (material_phaseAt(c,e))%partitionedState0(:,material_phaseMemberAt(c,i,e))

        do s = 1, phase_Nsources(material_phaseAt(c,e))
          sourceState(material_phaseAt(c,e))%p(s)%subState0(      :,material_phaseMemberAt(c,i,e)) = &
          sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phaseMemberAt(c,i,e))
        enddo
        crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_partitionedFp0(1:3,1:3,c,i,e)
        crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_partitionedFi0(1:3,1:3,c,i,e)
        crystallite_subF0(1:3,1:3,c,i,e)  = crystallite_partitionedF0(1:3,1:3,c,i,e)
        crystallite_subFrac(c,i,e) = 0.0_pReal
        crystallite_subStep(c,i,e) = 1.0_pReal/num%subStepSizeCryst
        todo(c,i,e) = .true.
        crystallite_converged(c,i,e) = .false.                                                      ! pretend failed step of 1/subStepSizeCryst
      endif homogenizationRequestsCalculation
    enddo; enddo
  enddo elementLooping1
  !$OMP END PARALLEL DO

  NiterationCrystallite = 0
  cutbackLooping: do while (any(todo(:,FEsolving_execIP(1):FEsolving_execIP(2),FEsolving_execELem(1):FEsolving_execElem(2))))
    NiterationCrystallite = NiterationCrystallite + 1

#ifdef DEBUG
    if (debugCrystallite%extensive) &
      print'(a,i6)', '<< CRYST stress >> crystallite iteration ',NiterationCrystallite
#endif
    !$OMP PARALLEL DO PRIVATE(formerSubStep)
    elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
      do i = FEsolving_execIP(1),FEsolving_execIP(2)
        do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
!--------------------------------------------------------------------------------------------------
!  wind forward
          if (crystallite_converged(c,i,e)) then
            formerSubStep = crystallite_subStep(c,i,e)
            crystallite_subFrac(c,i,e) = crystallite_subFrac(c,i,e) + crystallite_subStep(c,i,e)
            crystallite_subStep(c,i,e) = min(1.0_pReal - crystallite_subFrac(c,i,e), &
                                             num%stepIncreaseCryst * crystallite_subStep(c,i,e))

            todo(c,i,e) = crystallite_subStep(c,i,e) > 0.0_pReal                        ! still time left to integrate on?
            if (todo(c,i,e)) then
              crystallite_subF0 (1:3,1:3,c,i,e) = crystallite_subF(1:3,1:3,c,i,e)
              subLp0(1:3,1:3,c,i,e) = crystallite_Lp  (1:3,1:3,c,i,e)
              subLi0(1:3,1:3,c,i,e) = crystallite_Li  (1:3,1:3,c,i,e)
              crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_Fp  (1:3,1:3,c,i,e)
              crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_Fi  (1:3,1:3,c,i,e)
              plasticState(    material_phaseAt(c,e))%subState0(:,material_phaseMemberAt(c,i,e)) &
                = plasticState(material_phaseAt(c,e))%state(    :,material_phaseMemberAt(c,i,e))
              do s = 1, phase_Nsources(material_phaseAt(c,e))
                sourceState(    material_phaseAt(c,e))%p(s)%subState0(:,material_phaseMemberAt(c,i,e)) &
                  = sourceState(material_phaseAt(c,e))%p(s)%state(    :,material_phaseMemberAt(c,i,e))
              enddo
            endif

!--------------------------------------------------------------------------------------------------
!  cut back (reduced time and restore)
          else
            crystallite_subStep(c,i,e)       = num%subStepSizeCryst * crystallite_subStep(c,i,e)
            crystallite_Fp   (1:3,1:3,c,i,e) =            crystallite_subFp0(1:3,1:3,c,i,e)
            crystallite_Fi   (1:3,1:3,c,i,e) =            crystallite_subFi0(1:3,1:3,c,i,e)
            crystallite_S    (1:3,1:3,c,i,e) =            crystallite_S0    (1:3,1:3,c,i,e)
            if (crystallite_subStep(c,i,e) < 1.0_pReal) then                                        ! actual (not initial) cutback
              crystallite_Lp (1:3,1:3,c,i,e) =            subLp0(1:3,1:3,c,i,e)
              crystallite_Li (1:3,1:3,c,i,e) =            subLi0(1:3,1:3,c,i,e)
            endif
            plasticState    (material_phaseAt(c,e))%state(    :,material_phaseMemberAt(c,i,e)) &
              = plasticState(material_phaseAt(c,e))%subState0(:,material_phaseMemberAt(c,i,e))
            do s = 1, phase_Nsources(material_phaseAt(c,e))
              sourceState(    material_phaseAt(c,e))%p(s)%state(    :,material_phaseMemberAt(c,i,e)) &
                = sourceState(material_phaseAt(c,e))%p(s)%subState0(:,material_phaseMemberAt(c,i,e))
            enddo

                                                                                                    ! cant restore dotState here, since not yet calculated in first cutback after initialization
            todo(c,i,e) = crystallite_subStep(c,i,e) > num%subStepMinCryst                          ! still on track or already done (beyond repair)
          endif

!--------------------------------------------------------------------------------------------------
!  prepare for integration
          if (todo(c,i,e)) then
            crystallite_subF(1:3,1:3,c,i,e) = crystallite_subF0(1:3,1:3,c,i,e) &
                                            + crystallite_subStep(c,i,e) *( crystallite_partitionedF (1:3,1:3,c,i,e) &
                                                                           -crystallite_partitionedF0(1:3,1:3,c,i,e))
            crystallite_Fe(1:3,1:3,c,i,e) = matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                                   math_inv33(matmul(crystallite_Fi(1:3,1:3,c,i,e), &
                                                                     crystallite_Fp(1:3,1:3,c,i,e))))
            crystallite_subdt(c,i,e) = crystallite_subStep(c,i,e) * crystallite_dt(c,i,e)
            crystallite_converged(c,i,e) = .false.
            call integrateState(c,i,e)
          endif

        enddo
      enddo
    enddo elementLooping3
    !$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------
!  integrate --- requires fully defined state array (basic + dependent state)
    where(.not. crystallite_converged .and. crystallite_subStep > num%subStepMinCryst) &            ! do not try non-converged but fully cutbacked any further
      todo = .true.                                                                                 ! TODO: again unroll this into proper elementloop to avoid N^2 for single point evaluation


  enddo cutbackLooping

! return whether converged or not
  crystallite_stress = .false.
  elementLooping5: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2)
      crystallite_stress(i,e) = all(crystallite_converged(:,i,e))
    enddo
  enddo elementLooping5

end function crystallite_stress


!--------------------------------------------------------------------------------------------------
!> @brief Backup data for homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine crystallite_initializeRestorationPoints(i,e)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  integer :: &
    c, &                                                                                            !< constituent number
    s

  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
    crystallite_partitionedFp0(1:3,1:3,c,i,e) = crystallite_Fp0(1:3,1:3,c,i,e)
    crystallite_partitionedLp0(1:3,1:3,c,i,e) = crystallite_Lp0(1:3,1:3,c,i,e)
    crystallite_partitionedFi0(1:3,1:3,c,i,e) = crystallite_Fi0(1:3,1:3,c,i,e)
    crystallite_partitionedLi0(1:3,1:3,c,i,e) = crystallite_Li0(1:3,1:3,c,i,e)
    crystallite_partitionedF0(1:3,1:3,c,i,e)  = crystallite_F0(1:3,1:3,c,i,e)
    crystallite_partitionedS0(1:3,1:3,c,i,e)  = crystallite_S0(1:3,1:3,c,i,e)

    plasticState(material_phaseAt(c,e))%partitionedState0(:,material_phasememberAt(c,i,e)) = &
    plasticState(material_phaseAt(c,e))%state0(         :,material_phasememberAt(c,i,e))
    do s = 1, phase_Nsources(material_phaseAt(c,e))
      sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phasememberAt(c,i,e)) = &
      sourceState(material_phaseAt(c,e))%p(s)%state0(         :,material_phasememberAt(c,i,e))
    enddo
  enddo

end subroutine crystallite_initializeRestorationPoints


!--------------------------------------------------------------------------------------------------
!> @brief Wind homog inc forward.
!--------------------------------------------------------------------------------------------------
subroutine crystallite_windForward(i,e)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  integer :: &
    c, &                                                                                            !< constituent number
    s

  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
    crystallite_partitionedF0 (1:3,1:3,c,i,e) = crystallite_partitionedF(1:3,1:3,c,i,e)
    crystallite_partitionedFp0(1:3,1:3,c,i,e) = crystallite_Fp        (1:3,1:3,c,i,e)
    crystallite_partitionedLp0(1:3,1:3,c,i,e) = crystallite_Lp        (1:3,1:3,c,i,e)
    crystallite_partitionedFi0(1:3,1:3,c,i,e) = crystallite_Fi        (1:3,1:3,c,i,e)
    crystallite_partitionedLi0(1:3,1:3,c,i,e) = crystallite_Li        (1:3,1:3,c,i,e)
    crystallite_partitionedS0 (1:3,1:3,c,i,e) = crystallite_S         (1:3,1:3,c,i,e)

    plasticState    (material_phaseAt(c,e))%partitionedState0(:,material_phasememberAt(c,i,e)) = &
    plasticState    (material_phaseAt(c,e))%state          (:,material_phasememberAt(c,i,e))
    do s = 1, phase_Nsources(material_phaseAt(c,e))
      sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phasememberAt(c,i,e)) = &
      sourceState(material_phaseAt(c,e))%p(s)%state          (:,material_phasememberAt(c,i,e))
    enddo
  enddo

end subroutine crystallite_windForward


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restore(i,e,includeL)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  logical, intent(in) :: &
    includeL                                                                                        !< protect agains fake cutback
  integer :: &
    c, &                                                                                            !< constituent number
    s

  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
    if (includeL) then
      crystallite_Lp(1:3,1:3,c,i,e) = crystallite_partitionedLp0(1:3,1:3,c,i,e)
      crystallite_Li(1:3,1:3,c,i,e) = crystallite_partitionedLi0(1:3,1:3,c,i,e)
    endif                                                                                           ! maybe protecting everything from overwriting makes more sense
    crystallite_Fp(1:3,1:3,c,i,e)   = crystallite_partitionedFp0(1:3,1:3,c,i,e)
    crystallite_Fi(1:3,1:3,c,i,e)   = crystallite_partitionedFi0(1:3,1:3,c,i,e)
    crystallite_S (1:3,1:3,c,i,e)   = crystallite_partitionedS0 (1:3,1:3,c,i,e)

    plasticState    (material_phaseAt(c,e))%state(          :,material_phasememberAt(c,i,e)) = &
    plasticState    (material_phaseAt(c,e))%partitionedState0(:,material_phasememberAt(c,i,e))
    do s = 1, phase_Nsources(material_phaseAt(c,e))
      sourceState(material_phaseAt(c,e))%p(s)%state(          :,material_phasememberAt(c,i,e)) = &
      sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phasememberAt(c,i,e))
    enddo
  enddo

end subroutine crystallite_restore


!--------------------------------------------------------------------------------------------------
!> @brief Calculate tangent (dPdF).
!--------------------------------------------------------------------------------------------------
function crystallite_stressTangent(c,i,e) result(dPdF)

  real(pReal), dimension(3,3,3,3) :: dPdF
  integer, intent(in) :: &
    c, &                                                                                            !< counter in constituent loop
    i, &                                                                                            !< counter in integration point loop
    e                                                                                               !< counter in element loop
  integer :: &
    o, &
    p

  real(pReal), dimension(3,3)     ::   devNull, &
                                       invSubFp0,invSubFi0,invFp,invFi, &
                                       temp_33_1, temp_33_2, temp_33_3, temp_33_4
  real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                       dSdF, &
                                       dSdFi, &
                                       dLidS, &                                                     ! tangent in lattice configuration
                                       dLidFi, &
                                       dLpdS, &
                                       dLpdFi, &
                                       dFidS, &
                                       dFpinvdF, &
                                       rhs_3333, &
                                       lhs_3333, &
                                       temp_3333
  real(pReal), dimension(9,9)::        temp_99
  logical :: error


        call constitutive_SandItsTangents(devNull,dSdFe,dSdFi, &
                                         crystallite_Fe(1:3,1:3,c,i,e), &
                                         crystallite_Fi(1:3,1:3,c,i,e),c,i,e)
        call constitutive_LiAndItsTangents(devNull,dLidS,dLidFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e), &
                                           c,i,e)

        invFp = math_inv33(crystallite_Fp(1:3,1:3,c,i,e))
        invFi = math_inv33(crystallite_Fi(1:3,1:3,c,i,e))
        invSubFp0 = math_inv33(crystallite_subFp0(1:3,1:3,c,i,e))
        invSubFi0 = math_inv33(crystallite_subFi0(1:3,1:3,c,i,e))

        if (sum(abs(dLidS)) < tol_math_check) then
          dFidS = 0.0_pReal
        else
          lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
          do o=1,3; do p=1,3
            lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) &
                                  + crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidFi(1:3,1:3,o,p))
            lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) &
                                  + invFi*invFi(p,o)
            rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) &
                                  - crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidS(1:3,1:3,o,p))
          enddo; enddo
          call math_invert(temp_99,error,math_3333to99(lhs_3333))
          if (error) then
            call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                            ext_msg='inversion error in analytic tangent calculation')
            dFidS = 0.0_pReal
          else
            dFidS = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
          endif
          dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
        endif

        call constitutive_LpAndItsTangents(devNull,dLpdS,dLpdFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                     ! call constitutive law to calculate Lp tangent in lattice configuration
        dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

!--------------------------------------------------------------------------------------------------
! calculate dSdF
        temp_33_1 = transpose(matmul(invFp,invFi))
        temp_33_2 = matmul(crystallite_subF(1:3,1:3,c,i,e),invSubFp0)
        temp_33_3 = matmul(matmul(crystallite_subF(1:3,1:3,c,i,e),invFp), invSubFi0)

        do o=1,3; do p=1,3
          rhs_3333(p,o,1:3,1:3)  = matmul(dSdFe(p,o,1:3,1:3),temp_33_1)
          temp_3333(1:3,1:3,p,o) = matmul(matmul(temp_33_2,dLpdS(1:3,1:3,p,o)), invFi) &
                                 + matmul(temp_33_3,dLidS(1:3,1:3,p,o))
        enddo; enddo
        lhs_3333 = crystallite_subdt(c,i,e)*math_mul3333xx3333(dSdFe,temp_3333) &
                 + math_mul3333xx3333(dSdFi,dFidS)

        call math_invert(temp_99,error,math_eye(9)+math_3333to99(lhs_3333))
        if (error) then
          call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                          ext_msg='inversion error in analytic tangent calculation')
          dSdF = rhs_3333
        else
          dSdF = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
        endif

!--------------------------------------------------------------------------------------------------
! calculate dFpinvdF
        temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
        do o=1,3; do p=1,3
          dFpinvdF(1:3,1:3,p,o) = -crystallite_subdt(c,i,e) &
                                * matmul(invSubFp0, matmul(temp_3333(1:3,1:3,p,o),invFi))
        enddo; enddo

!--------------------------------------------------------------------------------------------------
! assemble dPdF
        temp_33_1 = matmul(crystallite_S(1:3,1:3,c,i,e),transpose(invFp))
        temp_33_2 = matmul(invFp,temp_33_1)
        temp_33_3 = matmul(crystallite_subF(1:3,1:3,c,i,e),invFp)
        temp_33_4 = matmul(temp_33_3,crystallite_S(1:3,1:3,c,i,e))

        dPdF = 0.0_pReal
        do p=1,3
          dPdF(p,1:3,p,1:3) = transpose(temp_33_2)
        enddo
        do o=1,3; do p=1,3
          dPdF(1:3,1:3,p,o) = dPdF(1:3,1:3,p,o) &
                            + matmul(matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                     dFpinvdF(1:3,1:3,p,o)),temp_33_1) &
                            + matmul(matmul(temp_33_3,dSdF(1:3,1:3,p,o)), &
                                     transpose(invFp)) &
                            + matmul(temp_33_4,transpose(dFpinvdF(1:3,1:3,p,o)))
        enddo; enddo

end function crystallite_stressTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations

  integer &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e                                                                                               !< counter in element loop

  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2)
      do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
        call crystallite_orientation(c,i,e)%fromMatrix(transpose(math_rotationalPart(crystallite_Fe(1:3,1:3,c,i,e))))
  enddo; enddo; enddo
  !$OMP END PARALLEL DO

  nonlocalPresent: if (any(plasticState%nonlocal)) then
    !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)
      if (plasticState(material_phaseAt(1,e))%nonlocal) then
        do i = FEsolving_execIP(1),FEsolving_execIP(2)
          call plastic_nonlocal_updateCompatibility(crystallite_orientation, &
                                                    phase_plasticityInstance(material_phaseAt(1,e)),i,e)
        enddo
      endif
    enddo
    !$OMP END PARALLEL DO
  endif nonlocalPresent

end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(ipc,ip,el, tensor33)

  real(pReal), dimension(3,3) :: crystallite_push33ToRef
  real(pReal), dimension(3,3), intent(in) :: tensor33
  real(pReal), dimension(3,3)             :: T
  integer, intent(in):: &
    el, &
    ip, &
    ipc

  T = matmul(material_orientation0(ipc,ip,el)%asMatrix(), &                                         ! ToDo: initial orientation correct?
             transpose(math_inv33(crystallite_subF(1:3,1:3,ipc,ip,el))))
  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief writes crystallite results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine crystallite_results

  integer :: p,o
  real(pReal),    allocatable, dimension(:,:,:) :: selected_tensors
  type(rotation), allocatable, dimension(:)     :: selected_rotations
  character(len=:),allocatable                  :: group,structureLabel

  do p=1,size(material_name_phase)
    group = trim('current/constituent')//'/'//trim(material_name_phase(p))//'/generic'

    call results_closeGroup(results_addGroup(group))

    do o = 1, size(output_constituent(p)%label)
      select case (output_constituent(p)%label(o))
        case('F')
          selected_tensors = select_tensors(crystallite_partitionedF,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'deformation gradient','1')
        case('F_e')
          selected_tensors = select_tensors(crystallite_Fe,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'elastic deformation gradient','1')
        case('F_p')
          selected_tensors = select_tensors(crystallite_Fp,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'plastic deformation gradient','1')
        case('F_i')
          selected_tensors = select_tensors(crystallite_Fi,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'inelastic deformation gradient','1')
        case('L_p')
          selected_tensors = select_tensors(crystallite_Lp,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'plastic velocity gradient','1/s')
        case('L_i')
          selected_tensors = select_tensors(crystallite_Li,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'inelastic velocity gradient','1/s')
        case('P')
          selected_tensors = select_tensors(crystallite_P,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'First Piola-Kirchoff stress','Pa')
        case('S')
          selected_tensors = select_tensors(crystallite_S,p)
          call results_writeDataset(group,selected_tensors,output_constituent(p)%label(o),&
                                   'Second Piola-Kirchoff stress','Pa')
        case('O')
          select case(lattice_structure(p))
            case(lattice_ISO_ID)
              structureLabel = 'iso'
            case(lattice_FCC_ID)
              structureLabel = 'fcc'
            case(lattice_BCC_ID)
              structureLabel = 'bcc'
            case(lattice_BCT_ID)
              structureLabel = 'bct'
            case(lattice_HEX_ID)
              structureLabel = 'hex'
            case(lattice_ORT_ID)
              structureLabel = 'ort'
          end select
          selected_rotations = select_rotations(crystallite_orientation,p)
          call results_writeDataset(group,selected_rotations,output_constituent(p)%label(o),&
                                   'crystal orientation as quaternion',structureLabel)
      end select
    enddo
  enddo

  contains

  !------------------------------------------------------------------------------------------------
  !> @brief select tensors for output
  !------------------------------------------------------------------------------------------------
  function select_tensors(dataset,instance)

    integer, intent(in) :: instance
    real(pReal), dimension(:,:,:,:,:), intent(in) :: dataset
    real(pReal), allocatable, dimension(:,:,:) :: select_tensors
    integer :: e,i,c,j

    allocate(select_tensors(3,3,count(material_phaseAt==instance)*discretization_nIPs))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, discretization_nIPs
        do c = 1, size(material_phaseAt,1)                                                          !ToDo: this needs to be changed for varying Ngrains
          if (material_phaseAt(c,e) == instance) then
            j = j + 1
            select_tensors(1:3,1:3,j) = dataset(1:3,1:3,c,i,e)
          endif
        enddo
      enddo
    enddo

  end function select_tensors


!--------------------------------------------------------------------------------------------------
!> @brief select rotations for output
!--------------------------------------------------------------------------------------------------
  function select_rotations(dataset,instance)

    integer, intent(in) :: instance
    type(rotation), dimension(:,:,:), intent(in) :: dataset
    type(rotation), allocatable, dimension(:) :: select_rotations
    integer :: e,i,c,j

    allocate(select_rotations(count(material_phaseAt==instance)*homogenization_maxNconstituents*discretization_nIPs))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, discretization_nIPs
        do c = 1, size(material_phaseAt,1)                                                          !ToDo: this needs to be changed for varying Ngrains
           if (material_phaseAt(c,e) == instance) then
             j = j + 1
             select_rotations(j) = dataset(c,i,e)
           endif
        enddo
      enddo
   enddo

 end function select_rotations

end subroutine crystallite_results


!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
function integrateStress(ipc,ip,el,timeFraction) result(broken)

  integer, intent(in)::         el, &                                                               ! element index
                                      ip, &                                                         ! integration point index
                                      ipc                                                           ! grain index
  real(pReal), optional, intent(in) :: timeFraction                                                 ! fraction of timestep

  real(pReal), dimension(3,3)::       F, &                                                          ! deformation gradient at end of timestep
                                      Fp_new, &                                                     ! plastic deformation gradient at end of timestep
                                      invFp_new, &                                                  ! inverse of Fp_new
                                      invFp_current, &                                              ! inverse of Fp_current
                                      Lpguess, &                                                    ! current guess for plastic velocity gradient
                                      Lpguess_old, &                                                ! known last good guess for plastic velocity gradient
                                      Lp_constitutive, &                                            ! plastic velocity gradient resulting from constitutive law
                                      residuumLp, &                                                 ! current residuum of plastic velocity gradient
                                      residuumLp_old, &                                             ! last residuum of plastic velocity gradient
                                      deltaLp, &                                                    ! direction of next guess
                                      Fi_new, &                                                     ! gradient of intermediate deformation stages
                                      invFi_new, &
                                      invFi_current, &                                              ! inverse of Fi_current
                                      Liguess, &                                                    ! current guess for intermediate velocity gradient
                                      Liguess_old, &                                                ! known last good guess for intermediate velocity gradient
                                      Li_constitutive, &                                            ! intermediate velocity gradient resulting from constitutive law
                                      residuumLi, &                                                 ! current residuum of intermediate velocity gradient
                                      residuumLi_old, &                                             ! last residuum of intermediate velocity gradient
                                      deltaLi, &                                                    ! direction of next guess
                                      Fe, &                                                         ! elastic deformation gradient
                                      S, &                                                          ! 2nd Piola-Kirchhoff Stress in plastic (lattice) configuration
                                      A, &
                                      B, &
                                      temp_33
  real(pReal), dimension(9) ::        temp_9                                                        ! needed for matrix inversion by LAPACK
  integer,     dimension(9) ::        devNull_9                                                     ! needed for matrix inversion by LAPACK
  real(pReal), dimension(9,9) ::      dRLp_dLp, &                                                   ! partial derivative of residuum (Jacobian for Newton-Raphson scheme)
                                      dRLi_dLi                                                      ! partial derivative of residuumI (Jacobian for Newton-Raphson scheme)
  real(pReal), dimension(3,3,3,3)::   dS_dFe, &                                                     ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dS_dFi, &
                                      dFe_dLp, &                                                    ! partial derivative of elastic deformation gradient
                                      dFe_dLi, &
                                      dFi_dLi, &
                                      dLp_dFi, &
                                      dLi_dFi, &
                                      dLp_dS, &
                                      dLi_dS
  real(pReal)                         steplengthLp, &
                                      steplengthLi, &
                                      dt, &                                                         ! time increment
                                      atol_Lp, &
                                      atol_Li, &
                                      devNull
  integer                             NiterationStressLp, &                                         ! number of stress integrations
                                      NiterationStressLi, &                                         ! number of inner stress integrations
                                      ierr, &                                                       ! error indicator for LAPACK
                                      o, &
                                      p, &
                                      jacoCounterLp, &
                                      jacoCounterLi                                                 ! counters to check for Jacobian update
  logical :: error,broken

  broken = .true.

  if (present(timeFraction)) then
    dt = crystallite_subdt(ipc,ip,el) * timeFraction
    F  = crystallite_subF0(1:3,1:3,ipc,ip,el) &
       + (crystallite_subF(1:3,1:3,ipc,ip,el) - crystallite_subF0(1:3,1:3,ipc,ip,el)) * timeFraction
  else
    dt = crystallite_subdt(ipc,ip,el)
    F  = crystallite_subF(1:3,1:3,ipc,ip,el)
  endif

  call constitutive_dependentState(crystallite_partitionedF(1:3,1:3,ipc,ip,el), &
                                   crystallite_Fp(1:3,1:3,ipc,ip,el),ipc,ip,el)

  Lpguess = crystallite_Lp(1:3,1:3,ipc,ip,el)                                                       ! take as first guess
  Liguess = crystallite_Li(1:3,1:3,ipc,ip,el)                                                       ! take as first guess

  call math_invert33(invFp_current,devNull,error,crystallite_subFp0(1:3,1:3,ipc,ip,el))
  if (error) return ! error
  call math_invert33(invFi_current,devNull,error,crystallite_subFi0(1:3,1:3,ipc,ip,el))
  if (error) return ! error

  A = matmul(F,invFp_current)                                                                       ! intermediate tensor needed later to calculate dFe_dLp

  jacoCounterLi  = 0
  steplengthLi   = 1.0_pReal
  residuumLi_old = 0.0_pReal
  Liguess_old    = Liguess

  NiterationStressLi = 0
  LiLoop: do
    NiterationStressLi = NiterationStressLi + 1
    if (NiterationStressLi>num%nStress) return ! error

    invFi_new = matmul(invFi_current,math_I3 - dt*Liguess)
    Fi_new    = math_inv33(invFi_new)

    jacoCounterLp  = 0
    steplengthLp   = 1.0_pReal
    residuumLp_old = 0.0_pReal
    Lpguess_old    = Lpguess

    NiterationStressLp = 0
    LpLoop: do
      NiterationStressLp = NiterationStressLp + 1
      if (NiterationStressLp>num%nStress) return ! error

      B  = math_I3 - dt*Lpguess
      Fe = matmul(matmul(A,B), invFi_new)
      call constitutive_SandItsTangents(S, dS_dFe, dS_dFi, &
                                        Fe, Fi_new, ipc, ip, el)

      call constitutive_LpAndItsTangents(Lp_constitutive, dLp_dS, dLp_dFi, &
                                         S, Fi_new, ipc, ip, el)

      !* update current residuum and check for convergence of loop
      atol_Lp = max(num%rtol_crystalliteStress * max(norm2(Lpguess),norm2(Lp_constitutive)), &      ! absolute tolerance from largest acceptable relative error
                    num%atol_crystalliteStress)                                                     ! minimum lower cutoff
      residuumLp = Lpguess - Lp_constitutive

      if (any(IEEE_is_NaN(residuumLp))) then
        return ! error
      elseif (norm2(residuumLp) < atol_Lp) then                                                     ! converged if below absolute tolerance
        exit LpLoop
      elseif (NiterationStressLp == 1 .or. norm2(residuumLp) < norm2(residuumLp_old)) then          ! not converged, but improved norm of residuum (always proceed in first iteration)...
        residuumLp_old = residuumLp                                                                 ! ...remember old values and...
        Lpguess_old    = Lpguess
        steplengthLp   = 1.0_pReal                                                                  ! ...proceed with normal step length (calculate new search direction)
      else                                                                                          ! not converged and residuum not improved...
        steplengthLp = num%subStepSizeLp * steplengthLp                                             ! ...try with smaller step length in same direction
        Lpguess      = Lpguess_old &
                     + deltaLp * stepLengthLp
        cycle LpLoop
      endif

      calculateJacobiLi: if (mod(jacoCounterLp, num%iJacoLpresiduum) == 0) then
        jacoCounterLp = jacoCounterLp + 1

        do o=1,3; do p=1,3
          dFe_dLp(o,1:3,p,1:3) = - dt * A(o,p)*transpose(invFi_new)                                 ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
        enddo; enddo
        dRLp_dLp = math_eye(9) &
                 - math_3333to99(math_mul3333xx3333(math_mul3333xx3333(dLp_dS,dS_dFe),dFe_dLp))
        temp_9 = math_33to9(residuumLp)
        call dgesv(9,1,dRLp_dLp,9,devNull_9,temp_9,9,ierr)                                          ! solve dRLp/dLp * delta Lp = -res for delta Lp
        if (ierr /= 0) return ! error
        deltaLp = - math_9to33(temp_9)
      endif calculateJacobiLi

      Lpguess = Lpguess &
              + deltaLp * steplengthLp
    enddo LpLoop

    call constitutive_LiAndItsTangents(Li_constitutive, dLi_dS, dLi_dFi, &
                                       S, Fi_new, ipc, ip, el)

    !* update current residuum and check for convergence of loop
    atol_Li = max(num%rtol_crystalliteStress * max(norm2(Liguess),norm2(Li_constitutive)), &        ! absolute tolerance from largest acceptable relative error
                  num%atol_crystalliteStress)                                                       ! minimum lower cutoff
    residuumLi = Liguess - Li_constitutive
    if (any(IEEE_is_NaN(residuumLi))) then
      return ! error
    elseif (norm2(residuumLi) < atol_Li) then                                                       ! converged if below absolute tolerance
      exit LiLoop
    elseif (NiterationStressLi == 1 .or. norm2(residuumLi) < norm2(residuumLi_old)) then            ! not converged, but improved norm of residuum (always proceed in first iteration)...
      residuumLi_old = residuumLi                                                                   ! ...remember old values and...
      Liguess_old    = Liguess
      steplengthLi   = 1.0_pReal                                                                    ! ...proceed with normal step length (calculate new search direction)
    else                                                                                            ! not converged and residuum not improved...
      steplengthLi = num%subStepSizeLi * steplengthLi                                               ! ...try with smaller step length in same direction
      Liguess      = Liguess_old &
                   + deltaLi * steplengthLi
      cycle LiLoop
    endif

    calculateJacobiLp: if (mod(jacoCounterLi, num%iJacoLpresiduum) == 0) then
      jacoCounterLi = jacoCounterLi + 1

      temp_33 = matmul(matmul(A,B),invFi_current)
      do o=1,3; do p=1,3
        dFe_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*temp_33                                             ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
        dFi_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*invFi_current
      enddo; enddo
      do o=1,3; do p=1,3
        dFi_dLi(1:3,1:3,o,p) = matmul(matmul(Fi_new,dFi_dLi(1:3,1:3,o,p)),Fi_new)
      enddo; enddo
      dRLi_dLi  = math_eye(9) &
                - math_3333to99(math_mul3333xx3333(dLi_dS,  math_mul3333xx3333(dS_dFe, dFe_dLi) &
                                                          + math_mul3333xx3333(dS_dFi, dFi_dLi)))  &
                - math_3333to99(math_mul3333xx3333(dLi_dFi, dFi_dLi))
      temp_9 = math_33to9(residuumLi)
      call dgesv(9,1,dRLi_dLi,9,devNull_9,temp_9,9,ierr)                                            ! solve dRLi/dLp * delta Li = -res for delta Li
      if (ierr /= 0) return ! error
      deltaLi = - math_9to33(temp_9)
    endif calculateJacobiLp

    Liguess = Liguess &
            + deltaLi * steplengthLi
  enddo LiLoop

  invFp_new = matmul(invFp_current,B)
  call math_invert33(Fp_new,devNull,error,invFp_new)
  if (error) return ! error

  crystallite_P    (1:3,1:3,ipc,ip,el) = matmul(matmul(F,invFp_new),matmul(S,transpose(invFp_new)))
  crystallite_S    (1:3,1:3,ipc,ip,el) = S
  crystallite_Lp   (1:3,1:3,ipc,ip,el) = Lpguess
  crystallite_Li   (1:3,1:3,ipc,ip,el) = Liguess
  crystallite_Fp   (1:3,1:3,ipc,ip,el) = Fp_new / math_det33(Fp_new)**(1.0_pReal/3.0_pReal)         ! regularize
  crystallite_Fi   (1:3,1:3,ipc,ip,el) = Fi_new
  crystallite_Fe   (1:3,1:3,ipc,ip,el) = matmul(matmul(F,invFp_new),invFi_new)
  broken = .false.

end function integrateStress


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
subroutine integrateStateFPI(g,i,e)

  integer, intent(in) :: &
    e, &                                                                                            !< element index in element loop
    i, &                                                                                            !< integration point index in ip loop
    g                                                                                               !< grain index in grain loop
  integer :: &
    NiterationState, &                                                                              !< number of iterations in state loop
    p, &
    c, &
    s, &
    size_pl
  integer, dimension(maxval(phase_Nsources)) :: &
    size_so
  real(pReal) :: &
    zeta
  real(pReal), dimension(max(constitutive_plasticity_maxSizeDotState,constitutive_source_maxSizeDotState)) :: &
    r                                                                                               ! state residuum
  real(pReal), dimension(constitutive_plasticity_maxSizeDotState,2) :: &
    plastic_dotState
  real(pReal), dimension(constitutive_source_maxSizeDotState,2,maxval(phase_Nsources)) :: source_dotState
  logical :: &
    broken

  p = material_phaseAt(g,e)
  c = material_phaseMemberAt(g,i,e)

  broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                    crystallite_partitionedF0, &
                                    crystallite_Fi(1:3,1:3,g,i,e), &
                                    crystallite_partitionedFp0, &
                                    crystallite_subdt(g,i,e), g,i,e,p,c)
  if(broken) return

  size_pl = plasticState(p)%sizeDotState
  plasticState(p)%state(1:size_pl,c) = plasticState(p)%subState0(1:size_pl,c) &
                                     + plasticState(p)%dotState (1:size_pl,c) &
                                     * crystallite_subdt(g,i,e)
  plastic_dotState(1:size_pl,2) = 0.0_pReal
  do s = 1, phase_Nsources(p)
    size_so(s) = sourceState(p)%p(s)%sizeDotState
    sourceState(p)%p(s)%state(1:size_so(s),c) = sourceState(p)%p(s)%subState0(1:size_so(s),c) &
                                              + sourceState(p)%p(s)%dotState (1:size_so(s),c) &
                                              * crystallite_subdt(g,i,e)
    source_dotState(1:size_so(s),2,s) = 0.0_pReal
  enddo

  iteration: do NiterationState = 1, num%nState

    if(nIterationState > 1) plastic_dotState(1:size_pl,2) = plastic_dotState(1:size_pl,1)
    plastic_dotState(1:size_pl,1) = plasticState(p)%dotState(:,c)
    do s = 1, phase_Nsources(p)
      if(nIterationState > 1) source_dotState(1:size_so(s),2,s) = source_dotState(1:size_so(s),1,s)
      source_dotState(1:size_so(s),1,s) = sourceState(p)%p(s)%dotState(:,c)
    enddo

    broken = integrateStress(g,i,e)
    if(broken) exit iteration

    broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                          crystallite_partitionedF0, &
                                          crystallite_Fi(1:3,1:3,g,i,e), &
                                          crystallite_partitionedFp0, &
                                          crystallite_subdt(g,i,e), g,i,e,p,c)
    if(broken) exit iteration

    zeta = damper(plasticState(p)%dotState(:,c),plastic_dotState(1:size_pl,1),&
                                                plastic_dotState(1:size_pl,2))
    plasticState(p)%dotState(:,c) = plasticState(p)%dotState(:,c) * zeta &
                                  + plastic_dotState(1:size_pl,1) * (1.0_pReal - zeta)
    r(1:size_pl) = plasticState(p)%state    (1:size_pl,c) &
                 - plasticState(p)%subState0(1:size_pl,c)  &
                 - plasticState(p)%dotState (1:size_pl,c) * crystallite_subdt(g,i,e)
    plasticState(p)%state(1:size_pl,c) = plasticState(p)%state(1:size_pl,c) &
                                       - r(1:size_pl)
    crystallite_converged(g,i,e) = converged(r(1:size_pl), &
                                             plasticState(p)%state(1:size_pl,c), &
                                             plasticState(p)%atol(1:size_pl))
    do s = 1, phase_Nsources(p)
      zeta = damper(sourceState(p)%p(s)%dotState(:,c), &
                    source_dotState(1:size_so(s),1,s),&
                    source_dotState(1:size_so(s),2,s))
      sourceState(p)%p(s)%dotState(:,c) = sourceState(p)%p(s)%dotState(:,c) * zeta &
                                        + source_dotState(1:size_so(s),1,s)* (1.0_pReal - zeta)
      r(1:size_so(s)) = sourceState(p)%p(s)%state    (1:size_so(s),c)  &
                      - sourceState(p)%p(s)%subState0(1:size_so(s),c)  &
                      - sourceState(p)%p(s)%dotState (1:size_so(s),c) * crystallite_subdt(g,i,e)
      sourceState(p)%p(s)%state(1:size_so(s),c) = sourceState(p)%p(s)%state(1:size_so(s),c) &
                                                - r(1:size_so(s))
      crystallite_converged(g,i,e) = &
      crystallite_converged(g,i,e) .and. converged(r(1:size_so(s)), &
                                                   sourceState(p)%p(s)%state(1:size_so(s),c), &
                                                   sourceState(p)%p(s)%atol(1:size_so(s)))
    enddo

    if(crystallite_converged(g,i,e)) then
      broken = constitutive_deltaState(crystallite_S(1:3,1:3,g,i,e), &
                                       crystallite_Fe(1:3,1:3,g,i,e), &
                                       crystallite_Fi(1:3,1:3,g,i,e),g,i,e,p,c)
      exit iteration
    endif

  enddo iteration


  contains

  !--------------------------------------------------------------------------------------------------
  !> @brief calculate the damping for correction of state and dot state
  !--------------------------------------------------------------------------------------------------
  real(pReal) pure function damper(current,previous,previous2)

  real(pReal), dimension(:), intent(in) ::&
    current, previous, previous2

  real(pReal) :: dot_prod12, dot_prod22

  dot_prod12 = dot_product(current  - previous,  previous - previous2)
  dot_prod22 = dot_product(previous - previous2, previous - previous2)
  if ((dot_product(current,previous) < 0.0_pReal .or. dot_prod12 < 0.0_pReal) .and. dot_prod22 > 0.0_pReal) then
    damper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
  else
    damper = 1.0_pReal
  endif

  end function damper

end subroutine integrateStateFPI


!--------------------------------------------------------------------------------------------------
!> @brief integrate state with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
subroutine integrateStateEuler(g,i,e)

  integer, intent(in) :: &
    e, &                                                                                            !< element index in element loop
    i, &                                                                                            !< integration point index in ip loop
    g                                                                                               !< grain index in grain loop
  integer :: &
    p, &
    c, &
    s, &
    sizeDotState
  logical :: &
    broken

  p = material_phaseAt(g,e)
  c = material_phaseMemberAt(g,i,e)

  broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                    crystallite_partitionedF0, &
                                    crystallite_Fi(1:3,1:3,g,i,e), &
                                    crystallite_partitionedFp0, &
                                    crystallite_subdt(g,i,e), g,i,e,p,c)
  if(broken) return

  sizeDotState = plasticState(p)%sizeDotState
  plasticState(p)%state(1:sizeDotState,c) = plasticState(p)%subState0(1:sizeDotState,c) &
                                          + plasticState(p)%dotState (1:sizeDotState,c) &
                                            * crystallite_subdt(g,i,e)
  do s = 1, phase_Nsources(p)
    sizeDotState = sourceState(p)%p(s)%sizeDotState
    sourceState(p)%p(s)%state(1:sizeDotState,c) = sourceState(p)%p(s)%subState0(1:sizeDotState,c) &
                                                + sourceState(p)%p(s)%dotState (1:sizeDotState,c) &
                                                  * crystallite_subdt(g,i,e)
  enddo

  broken = constitutive_deltaState(crystallite_S(1:3,1:3,g,i,e), &
                                   crystallite_Fe(1:3,1:3,g,i,e), &
                                  crystallite_Fi(1:3,1:3,g,i,e),g,i,e,p,c)
  if(broken) return

  broken = integrateStress(g,i,e)
  crystallite_converged(g,i,e) = .not. broken

end subroutine integrateStateEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 1st order Euler method with adaptive step size
!--------------------------------------------------------------------------------------------------
subroutine integrateStateAdaptiveEuler(g,i,e)

  integer, intent(in) :: &
    e, &                                                                                            !< element index in element loop
    i, &                                                                                            !< integration point index in ip loop
    g                                                                                               !< grain index in grain loop
  integer :: &
    p, &
    c, &
    s, &
    sizeDotState
  logical :: &
    broken

  real(pReal), dimension(constitutive_plasticity_maxSizeDotState) :: residuum_plastic
  real(pReal), dimension(constitutive_source_maxSizeDotState,maxval(phase_Nsources)) :: residuum_source


  p = material_phaseAt(g,e)
  c = material_phaseMemberAt(g,i,e)

  broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                        crystallite_partitionedF0, &
                                        crystallite_Fi(1:3,1:3,g,i,e), &
                                        crystallite_partitionedFp0, &
                                        crystallite_subdt(g,i,e), g,i,e,p,c)
  if(broken) return

  sizeDotState = plasticState(p)%sizeDotState

  residuum_plastic(1:sizeDotState) = - plasticState(p)%dotstate(1:sizeDotState,c) * 0.5_pReal * crystallite_subdt(g,i,e)
  plasticState(p)%state(1:sizeDotState,c) = plasticState(p)%subState0(1:sizeDotState,c) &
                                          + plasticState(p)%dotstate(1:sizeDotState,c) * crystallite_subdt(g,i,e)
  do s = 1, phase_Nsources(p)
    sizeDotState = sourceState(p)%p(s)%sizeDotState

    residuum_source(1:sizeDotState,s)  = - sourceState(p)%p(s)%dotstate(1:sizeDotState,c) &
                                       * 0.5_pReal * crystallite_subdt(g,i,e)
    sourceState(p)%p(s)%state(1:sizeDotState,c) = sourceState(p)%p(s)%subState0(1:sizeDotState,c) &
                                                + sourceState(p)%p(s)%dotstate(1:sizeDotState,c) * crystallite_subdt(g,i,e)
  enddo

  broken = constitutive_deltaState(crystallite_S(1:3,1:3,g,i,e), &
                                   crystallite_Fe(1:3,1:3,g,i,e), &
                                   crystallite_Fi(1:3,1:3,g,i,e),g,i,e,p,c)
  if(broken) return

  broken = integrateStress(g,i,e)
  if(broken) return

  broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                        crystallite_partitionedF0, &
                                        crystallite_Fi(1:3,1:3,g,i,e), &
                                        crystallite_partitionedFp0, &
                                        crystallite_subdt(g,i,e), g,i,e,p,c)
  if(broken) return


  sizeDotState = plasticState(p)%sizeDotState
  crystallite_converged(g,i,e) = converged(residuum_plastic(1:sizeDotState) &
                                           + 0.5_pReal * plasticState(p)%dotState(:,c) * crystallite_subdt(g,i,e), &
                                           plasticState(p)%state(1:sizeDotState,c), &
                                           plasticState(p)%atol(1:sizeDotState))

  do s = 1, phase_Nsources(p)
    sizeDotState = sourceState(p)%p(s)%sizeDotState
    crystallite_converged(g,i,e) = &
    crystallite_converged(g,i,e) .and. converged(residuum_source(1:sizeDotState,s) &
                                                 + 0.5_pReal*sourceState(p)%p(s)%dotState(:,c)*crystallite_subdt(g,i,e), &
                                                 sourceState(p)%p(s)%state(1:sizeDotState,c), &
                                                 sourceState(p)%p(s)%atol(1:sizeDotState))
  enddo

end subroutine integrateStateAdaptiveEuler


!---------------------------------------------------------------------------------------------------
!> @brief Integrate state (including stress integration) with the classic Runge Kutta method
!---------------------------------------------------------------------------------------------------
subroutine integrateStateRK4(g,i,e)

  integer, intent(in) :: g,i,e

  real(pReal), dimension(3,3), parameter :: &
    A = reshape([&
      0.5_pReal, 0.0_pReal, 0.0_pReal, &
      0.0_pReal, 0.5_pReal, 0.0_pReal, &
      0.0_pReal, 0.0_pReal, 1.0_pReal],&
      shape(A))
  real(pReal), dimension(3), parameter :: &
    C = [0.5_pReal, 0.5_pReal, 1.0_pReal]
  real(pReal), dimension(4), parameter :: &
    B = [1.0_pReal/6.0_pReal, 1.0_pReal/3.0_pReal, 1.0_pReal/3.0_pReal, 1.0_pReal/6.0_pReal]

  call integrateStateRK(g,i,e,A,B,C)

end subroutine integrateStateRK4


!---------------------------------------------------------------------------------------------------
!> @brief Integrate state (including stress integration) with the Cash-Carp method
!---------------------------------------------------------------------------------------------------
subroutine integrateStateRKCK45(g,i,e)

  integer, intent(in) :: g,i,e

  real(pReal), dimension(5,5), parameter :: &
    A = reshape([&
      1._pReal/5._pReal,       .0_pReal,             .0_pReal,               .0_pReal,                  .0_pReal, &
      3._pReal/40._pReal,      9._pReal/40._pReal,   .0_pReal,               .0_pReal,                  .0_pReal, &
      3_pReal/10._pReal,       -9._pReal/10._pReal,  6._pReal/5._pReal,      .0_pReal,                  .0_pReal, &
      -11._pReal/54._pReal,    5._pReal/2._pReal,    -70.0_pReal/27.0_pReal, 35.0_pReal/27.0_pReal,     .0_pReal, &
      1631._pReal/55296._pReal,175._pReal/512._pReal,575._pReal/13824._pReal,44275._pReal/110592._pReal,253._pReal/4096._pReal],&
      shape(A))
  real(pReal), dimension(5), parameter :: &
    C = [0.2_pReal, 0.3_pReal, 0.6_pReal, 1.0_pReal, 0.875_pReal]
  real(pReal), dimension(6), parameter :: &
    B = &
      [37.0_pReal/378.0_pReal, .0_pReal, 250.0_pReal/621.0_pReal, &
      125.0_pReal/594.0_pReal, .0_pReal, 512.0_pReal/1771.0_pReal], &
    DB = B - &
      [2825.0_pReal/27648.0_pReal,    .0_pReal,                18575.0_pReal/48384.0_pReal,&
      13525.0_pReal/55296.0_pReal, 277.0_pReal/14336.0_pReal,  1._pReal/4._pReal]

  call integrateStateRK(g,i,e,A,B,C,DB)

end subroutine integrateStateRKCK45


!--------------------------------------------------------------------------------------------------
!> @brief Integrate state (including stress integration) with an explicit Runge-Kutta method or an
!! embedded explicit Runge-Kutta method
!--------------------------------------------------------------------------------------------------
subroutine integrateStateRK(g,i,e,A,B,CC,DB)


  real(pReal), dimension(:,:), intent(in) :: A
  real(pReal), dimension(:),   intent(in) :: B, CC
  real(pReal), dimension(:),   intent(in), optional :: DB

  integer, intent(in) :: &
    e, &                                                                                            !< element index in element loop
    i, &                                                                                            !< integration point index in ip loop
    g                                                                                               !< grain index in grain loop
  integer :: &
    stage, &                                                                                        ! stage index in integration stage loop
    n, &
    p, &
    c, &
    s, &
    sizeDotState
  logical :: &
    broken
  real(pReal), dimension(constitutive_source_maxSizeDotState,size(B),maxval(phase_Nsources)) :: source_RKdotState
  real(pReal), dimension(constitutive_plasticity_maxSizeDotState,size(B))                    :: plastic_RKdotState

  p = material_phaseAt(g,e)
  c = material_phaseMemberAt(g,i,e)

  broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                        crystallite_partitionedF0, &
                                        crystallite_Fi(1:3,1:3,g,i,e), &
                                        crystallite_partitionedFp0, &
                                        crystallite_subdt(g,i,e), g,i,e,p,c)
  if(broken) return

  do stage = 1,size(A,1)
    sizeDotState = plasticState(p)%sizeDotState
    plastic_RKdotState(1:sizeDotState,stage) = plasticState(p)%dotState(:,c)
    plasticState(p)%dotState(:,c) = A(1,stage) * plastic_RKdotState(1:sizeDotState,1)
    do s = 1, phase_Nsources(p)
      sizeDotState = sourceState(p)%p(s)%sizeDotState
      source_RKdotState(1:sizeDotState,stage,s) = sourceState(p)%p(s)%dotState(:,c)
      sourceState(p)%p(s)%dotState(:,c) = A(1,stage) * source_RKdotState(1:sizeDotState,1,s)
    enddo

    do n = 2, stage
      sizeDotState = plasticState(p)%sizeDotState
      plasticState(p)%dotState(:,c) = plasticState(p)%dotState(:,c) &
                                    + A(n,stage) * plastic_RKdotState(1:sizeDotState,n)
      do s = 1, phase_Nsources(p)
        sizeDotState = sourceState(p)%p(s)%sizeDotState
        sourceState(p)%p(s)%dotState(:,c) = sourceState(p)%p(s)%dotState(:,c) &
                                          + A(n,stage) * source_RKdotState(1:sizeDotState,n,s)
      enddo
    enddo

    sizeDotState = plasticState(p)%sizeDotState
    plasticState(p)%state(1:sizeDotState,c) = plasticState(p)%subState0(1:sizeDotState,c) &
                                            + plasticState(p)%dotState (1:sizeDotState,c) &
                                              * crystallite_subdt(g,i,e)
    do s = 1, phase_Nsources(p)
      sizeDotState = sourceState(p)%p(s)%sizeDotState
      sourceState(p)%p(s)%state(1:sizeDotState,c) = sourceState(p)%p(s)%subState0(1:sizeDotState,c) &
                                                  + sourceState(p)%p(s)%dotState (1:sizeDotState,c) &
                                                    * crystallite_subdt(g,i,e)
    enddo

    broken = integrateStress(g,i,e,CC(stage))
    if(broken) exit

    broken = constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                          crystallite_partitionedF0, &
                                          crystallite_Fi(1:3,1:3,g,i,e), &
                                          crystallite_partitionedFp0, &
                                          crystallite_subdt(g,i,e)*CC(stage), g,i,e,p,c)
    if(broken) exit

  enddo
  if(broken) return

  sizeDotState = plasticState(p)%sizeDotState

  plastic_RKdotState(1:sizeDotState,size(B)) = plasticState (p)%dotState(:,c)
  plasticState(p)%dotState(:,c) = matmul(plastic_RKdotState(1:sizeDotState,1:size(B)),B)
  plasticState(p)%state(1:sizeDotState,c) = plasticState(p)%subState0(1:sizeDotState,c) &
                                          + plasticState(p)%dotState (1:sizeDotState,c) &
                                            * crystallite_subdt(g,i,e)
  if(present(DB)) &
    broken = .not. converged( matmul(plastic_RKdotState(1:sizeDotState,1:size(DB)),DB) &
                                             * crystallite_subdt(g,i,e), &
                                        plasticState(p)%state(1:sizeDotState,c), &
                                        plasticState(p)%atol(1:sizeDotState))

  do s = 1, phase_Nsources(p)
    sizeDotState = sourceState(p)%p(s)%sizeDotState

    source_RKdotState(1:sizeDotState,size(B),s) = sourceState(p)%p(s)%dotState(:,c)
    sourceState(p)%p(s)%dotState(:,c)  = matmul(source_RKdotState(1:sizeDotState,1:size(B),s),B)
    sourceState(p)%p(s)%state(1:sizeDotState,c) = sourceState(p)%p(s)%subState0(1:sizeDotState,c) &
                                                + sourceState(p)%p(s)%dotState (1:sizeDotState,c) &
                                                  * crystallite_subdt(g,i,e)
    if(present(DB)) &
      broken = broken .or. .not. converged(matmul(source_RKdotState(1:sizeDotState,1:size(DB),s),DB) &
                                                 * crystallite_subdt(g,i,e), &
                                           sourceState(p)%p(s)%state(1:sizeDotState,c), &
                                           sourceState(p)%p(s)%atol(1:sizeDotState))
  enddo
  if(broken) return

  broken = constitutive_deltaState(crystallite_S(1:3,1:3,g,i,e), &
                                   crystallite_Fe(1:3,1:3,g,i,e), &
                                   crystallite_Fi(1:3,1:3,g,i,e),g,i,e,p,c)
  if(broken) return

  broken = integrateStress(g,i,e)
  crystallite_converged(g,i,e) = .not. broken


end subroutine integrateStateRK


!--------------------------------------------------------------------------------------------------
!> @brief determines whether a point is converged
!--------------------------------------------------------------------------------------------------
logical pure function converged(residuum,state,atol)

  real(pReal), intent(in), dimension(:) ::&
    residuum, state, atol
  real(pReal) :: &
    rTol

  rTol = num%rTol_crystalliteState

  converged = all(abs(residuum) <= max(atol, rtol*abs(state)))

end function converged


!--------------------------------------------------------------------------------------------------
!> @brief Write current  restart information (Field and constitutive data) to file.
! ToDo: Merge data into one file for MPI, move state to constitutive and homogenization, respectively
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restartWrite

  integer :: i
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  print*, ' writing field and constitutive data required for restart to file';flush(IO_STDOUT)

  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName,'a')

  call HDF5_write(fileHandle,crystallite_partitionedF,'F')
  call HDF5_write(fileHandle,crystallite_Fp,        'F_p')
  call HDF5_write(fileHandle,crystallite_Fi,        'F_i')
  call HDF5_write(fileHandle,crystallite_Lp,        'L_p')
  call HDF5_write(fileHandle,crystallite_Li,        'L_i')
  call HDF5_write(fileHandle,crystallite_S,           'S')

  groupHandle = HDF5_addGroup(fileHandle,'constituent')
  do i = 1,size(material_name_phase)
    write(datasetName,'(i0,a)') i,'_omega_plastic'
    call HDF5_write(groupHandle,plasticState(i)%state,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_addGroup(fileHandle,'materialpoint')
  do i = 1, size(material_name_homogenization)
    write(datasetName,'(i0,a)') i,'_omega_homogenization'
    call HDF5_write(groupHandle,homogState(i)%state,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  call HDF5_closeFile(fileHandle)

end subroutine crystallite_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read data for restart
! ToDo: Merge data into one file for MPI, move state to constitutive and homogenization, respectively
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restartRead

  integer :: i
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  print'(/,a,i0,a)', ' reading restart information of increment from file'

  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName)

  call HDF5_read(fileHandle,crystallite_F0, 'F')
  call HDF5_read(fileHandle,crystallite_Fp0,'F_p')
  call HDF5_read(fileHandle,crystallite_Fi0,'F_i')
  call HDF5_read(fileHandle,crystallite_Lp0,'L_p')
  call HDF5_read(fileHandle,crystallite_Li0,'L_i')
  call HDF5_read(fileHandle,crystallite_S0, 'S')

  groupHandle = HDF5_openGroup(fileHandle,'constituent')
  do i = 1,size(material_name_phase)
    write(datasetName,'(i0,a)') i,'_omega_plastic'
    call HDF5_read(groupHandle,plasticState(i)%state0,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_openGroup(fileHandle,'materialpoint')
  do i = 1,size(material_name_homogenization)
    write(datasetName,'(i0,a)') i,'_omega_homogenization'
    call HDF5_read(groupHandle,homogState(i)%state0,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  call HDF5_closeFile(fileHandle)

end subroutine crystallite_restartRead


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine crystallite_forward

  integer :: i, j

  crystallite_F0  = crystallite_partitionedF
  crystallite_Fp0 = crystallite_Fp
  crystallite_Lp0 = crystallite_Lp
  crystallite_Fi0 = crystallite_Fi
  crystallite_Li0 = crystallite_Li
  crystallite_S0  = crystallite_S

  do i = 1, size(plasticState)
    plasticState(i)%state0 = plasticState(i)%state
  enddo
  do i = 1, size(sourceState)
    do j = 1,phase_Nsources(i)
      sourceState(i)%p(j)%state0 = sourceState(i)%p(j)%state
  enddo; enddo
  do i = 1,size(material_name_homogenization)
    homogState  (i)%state0 = homogState  (i)%state
    thermalState(i)%state0 = thermalState(i)%state
    damageState (i)%state0 = damageState (i)%state
  enddo

end subroutine crystallite_forward

end module crystallite
