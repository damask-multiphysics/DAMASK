!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, damage & thermal internal microstructure state
!--------------------------------------------------------------------------------------------------
module phase
  use prec
  use constants
  use math
  use rotations
  use polynomials
  use IO
  use config
  use material
  use results
  use lattice
  use discretization
  use parallelization
  use HDF5
  use HDF5_utilities

  implicit none
  private


  character(len=2), allocatable, dimension(:) :: phase_lattice
  real(pReal),      allocatable, dimension(:) :: phase_cOverA
  real(pReal),      allocatable, dimension(:) :: phase_rho

  type(tRotationContainer), dimension(:), allocatable :: &
    phase_O_0, &
    phase_O

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

  type(tPlasticState), allocatable, dimension(:), public :: &
    plasticState
  type(tState),  allocatable, dimension(:), public :: &
    damageState


  interface

! == cleaned:begin =================================================================================
    module subroutine mechanical_init(phases)
      class(tNode), pointer :: phases
    end subroutine mechanical_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine thermal_init(phases)
      class(tNode), pointer :: phases
    end subroutine thermal_init


    module subroutine mechanical_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine mechanical_results

    module subroutine damage_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine damage_results

    module subroutine thermal_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine thermal_results

    module subroutine mechanical_forward()
    end subroutine mechanical_forward

    module subroutine damage_forward()
    end subroutine damage_forward

    module subroutine thermal_forward()
    end subroutine thermal_forward


    module subroutine mechanical_restore(ce,includeL)
      integer, intent(in) :: ce
      logical, intent(in) :: includeL
    end subroutine mechanical_restore

    module subroutine damage_restore(ce)
      integer, intent(in) :: ce
    end subroutine damage_restore


    module function phase_mechanical_dPdF(Delta_t,co,ce) result(dPdF)
      real(pReal), intent(in) :: Delta_t
      integer, intent(in) :: &
        co, &                                                                                       !< counter in constituent loop
        ce
      real(pReal), dimension(3,3,3,3) :: dPdF
    end function phase_mechanical_dPdF

    module subroutine mechanical_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mechanical_restartWrite

    module subroutine thermal_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine thermal_restartWrite

    module subroutine mechanical_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mechanical_restartRead

    module subroutine thermal_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine thermal_restartRead

    module function mechanical_S(ph,en) result(S)
      integer, intent(in) :: ph,en
      real(pReal), dimension(3,3) :: S
    end function mechanical_S

    module function mechanical_L_p(ph,en) result(L_p)
      integer, intent(in) :: ph,en
      real(pReal), dimension(3,3) :: L_p
    end function mechanical_L_p

    module function mechanical_F_e(ph,en) result(F_e)
      integer, intent(in) :: ph,en
      real(pReal), dimension(3,3) :: F_e
    end function mechanical_F_e


    module function phase_F(co,ce) result(F)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: F
    end function phase_F

    module function phase_P(co,ce) result(P)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: P
    end function phase_P

    pure module function thermal_T(ph,en) result(T)
      integer, intent(in) :: ph,en
      real(pReal) :: T
    end function thermal_T

    module function thermal_dot_T(ph,en) result(dot_T)
      integer, intent(in) :: ph,en
      real(pReal) :: dot_T
    end function thermal_dot_T

    module function damage_phi(ph,en) result(phi)
      integer, intent(in) :: ph,en
      real(pReal) :: phi
    end function damage_phi


    module subroutine phase_set_F(F,co,ce)
      real(pReal), dimension(3,3), intent(in) :: F
      integer, intent(in) :: co, ce
    end subroutine phase_set_F

    module subroutine phase_thermal_setField(T,dot_T, co,ce)
      real(pReal), intent(in) :: T, dot_T
      integer, intent(in) :: co, ce
    end subroutine phase_thermal_setField

    module subroutine phase_set_phi(phi,co,ce)
      real(pReal), intent(in) :: phi
      integer, intent(in) :: co, ce
    end subroutine phase_set_phi


    module function phase_mu_phi(co,ce) result(mu)
      integer, intent(in) :: co, ce
      real(pReal) :: mu
    end function phase_mu_phi

    module function phase_K_phi(co,ce) result(K)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: K
    end function phase_K_phi


    module function phase_mu_T(co,ce) result(mu)
      integer, intent(in) :: co, ce
      real(pReal) :: mu
    end function phase_mu_T

    module function phase_K_T(co,ce) result(K)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: K
    end function phase_K_T

! == cleaned:end ===================================================================================

    module function phase_thermal_constitutive(Delta_t,ph,en) result(converged_)

      real(pReal), intent(in) :: Delta_t
      integer, intent(in) :: ph, en
      logical :: converged_

    end function phase_thermal_constitutive

    module function phase_damage_constitutive(Delta_t,co,ce) result(converged_)
      real(pReal), intent(in) :: Delta_t
      integer, intent(in) :: co, ce
      logical :: converged_
    end function phase_damage_constitutive

    module function phase_mechanical_constitutive(Delta_t,co,ce) result(converged_)
      real(pReal), intent(in) :: Delta_t
      integer, intent(in) :: co, ce
      logical :: converged_
    end function phase_mechanical_constitutive

    !ToDo: Merge all the stiffness functions
    module function phase_homogenizedC66(ph,en) result(C)
      integer, intent(in) :: ph, en
      real(pReal), dimension(6,6) :: C
    end function phase_homogenizedC66
    module function phase_damage_C66(C66,ph,en) result(C66_degraded)
      real(pReal), dimension(6,6), intent(in)  :: C66
      integer,                     intent(in)  :: ph,en
      real(pReal), dimension(6,6) :: C66_degraded
    end function phase_damage_C66

    module function phase_f_phi(phi,co,ce) result(f)
      integer, intent(in) :: ce,co
      real(pReal), intent(in) :: &
        phi                                                                                         !< damage parameter
      real(pReal) :: &
        f
    end function phase_f_phi

    module function phase_f_T(ph,en) result(f)
      integer, intent(in) :: ph, en
      real(pReal) :: f
    end function phase_f_T

    module subroutine plastic_nonlocal_updateCompatibility(orientation,ph,ip,el)
      integer, intent(in) :: &
        ph, &
        ip, &
        el
        type(tRotationContainer), dimension(:), intent(in) :: orientation
    end subroutine plastic_nonlocal_updateCompatibility

    module subroutine plastic_dependentState(ph,en)
      integer, intent(in) :: &
        ph, &
        en
    end subroutine plastic_dependentState


    module subroutine damage_anisobrittle_LiAndItsTangent(Ld, dLd_dTstar, S, ph,en)
      integer, intent(in) :: ph, en
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine damage_anisobrittle_LiAndItsTangent

  end interface


  type(tDebugOptions) :: debugConstitutive
#if __INTEL_COMPILER >= 1900
  public :: &
    prec, &
    math, &
    rotations, &
    IO, &
    config, &
    material, &
    results, &
    lattice, &
    discretization, &
    HDF5_utilities
#endif

  public :: &
    phase_init, &
    phase_homogenizedC66, &
    phase_f_phi, &
    phase_f_T, &
    phase_K_phi, &
    phase_K_T, &
    phase_mu_phi, &
    phase_mu_T, &
    phase_results, &
    phase_allocateState, &
    phase_forward, &
    phase_restore, &
    plastic_nonlocal_updateCompatibility, &
    converged, &
    phase_mechanical_constitutive, &
    phase_thermal_constitutive, &
    phase_damage_constitutive, &
    phase_mechanical_dPdF, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    phase_restartWrite, &
    phase_restartRead, &
    phase_thermal_setField, &
    phase_set_phi, &
    phase_P, &
    phase_set_F, &
    phase_F

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize constitutive models for individual physics
!--------------------------------------------------------------------------------------------------
subroutine phase_init

  integer :: &
    ph, ce, co, ma
  class (tNode), pointer :: &
    debug_constitutive, &
    materials, &
    phases, &
    phase


  print'(/,1x,a)', '<<<+-  phase init  -+>>>'; flush(IO_STDOUT)

  debug_constitutive => config_debug%get('phase', defaultVal=emptyList)
  debugConstitutive%basic     = debug_constitutive%contains('basic')
  debugConstitutive%extensive = debug_constitutive%contains('extensive')
  debugConstitutive%selective = debug_constitutive%contains('selective')
  debugConstitutive%element   = config_debug%get_asInt('element',         defaultVal = 1)
  debugConstitutive%ip        = config_debug%get_asInt('integrationpoint',defaultVal = 1)
  debugConstitutive%grain     = config_debug%get_asInt('constituent',     defaultVal = 1)


  materials => config_material%get('material')
  phases    => config_material%get('phase')

  allocate(phase_lattice(phases%length))
  allocate(phase_cOverA(phases%length),source=-1.0_pReal)
  allocate(phase_rho(phases%length))
  allocate(phase_O_0(phases%length))

  do ph = 1,phases%length
    phase => phases%get(ph)
    phase_lattice(ph) = phase%get_asString('lattice')
    if (all(phase_lattice(ph) /= ['cF','cI','hP','tI'])) &
      call IO_error(130,ext_msg='phase_init: '//phase%get_asString('lattice'))
    if (any(phase_lattice(ph) == ['hP','tI'])) &
      phase_cOverA(ph) = phase%get_asFloat('c/a')
    phase_rho(ph) = phase%get_asFloat('rho',defaultVal=0.0_pReal)
    allocate(phase_O_0(ph)%data(count(material_phaseID==ph)))
  end do

  do ce = 1, size(material_phaseID,2)
    ma = discretization_materialAt((ce-1)/discretization_nIPs+1)
    do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
      ph = material_phaseID(co,ce)
      phase_O_0(ph)%data(material_phaseEntry(co,ce)) = material_O_0(ma)%data(co)
    end do
  end do

  allocate(phase_O(phases%length))
  do ph = 1,phases%length
    phase_O(ph)%data = phase_O_0(ph)%data
  end do

  call mechanical_init(phases)
  call damage_init
  call thermal_init(phases)

  call crystallite_init()

end subroutine phase_init


!--------------------------------------------------------------------------------------------------
!> @brief Allocate the components of the state structure for a given phase
!--------------------------------------------------------------------------------------------------
subroutine phase_allocateState(state, &
                               NEntries,sizeState,sizeDotState,sizeDeltaState,offsetDeltaState)

  class(tState), intent(inout) :: &
    state
  integer, intent(in) :: &
    NEntries, &
    sizeState, &
    sizeDotState, &
    sizeDeltaState
  integer, intent(in), optional :: &
    offsetDeltaState

  state%sizeState        = sizeState
  state%sizeDotState     = sizeDotState
  state%sizeDeltaState   = sizeDeltaState
  if (present(offsetDeltaState)) then
    state%offsetDeltaState = offsetDeltaState                                                       ! ToDo: this is a fix for broken nonlocal
  else
    state%offsetDeltaState = sizeState-sizeDeltaState                                               ! deltaState occupies latter part of state by definition
  end if

  allocate(state%atol             (sizeState),          source=0.0_pReal)
  allocate(state%state0           (sizeState,NEntries), source=0.0_pReal)
  allocate(state%state            (sizeState,NEntries), source=0.0_pReal)

  allocate(state%dotState      (sizeDotState,NEntries), source=0.0_pReal)

  allocate(state%deltaState  (sizeDeltaState,NEntries), source=0.0_pReal)
  state%deltaState2 => state%state(state%offsetDeltaState+1: &
                                   state%offsetDeltaState+state%sizeDeltaState,:)

end subroutine phase_allocateState


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine phase_restore(ce,includeL)

  logical, intent(in) :: includeL
  integer, intent(in) :: ce


  call mechanical_restore(ce,includeL)
  call damage_restore(ce)

end subroutine phase_restore


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
!--------------------------------------------------------------------------------------------------
subroutine phase_forward()

  call mechanical_forward()
  call damage_forward()
  call thermal_forward()

end subroutine phase_forward


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine phase_results()

  integer :: ph
  character(len=:), allocatable :: group


  call results_closeGroup(results_addGroup('/current/phase/'))

  do ph = 1, size(material_name_phase)

    group = '/current/phase/'//trim(material_name_phase(ph))//'/'
    call results_closeGroup(results_addGroup(group))

    call mechanical_results(group,ph)
    call damage_results(group,ph)
    call thermal_results(group,ph)

  end do

end subroutine phase_results


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init()

  integer :: &
    ce, &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el, &                                                                                           !< counter in element loop
    en, ph
  class(tNode), pointer :: &
    num_crystallite, &
    phases


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

  if (num%subStepMinCryst   <= 0.0_pReal)      call IO_error(301,ext_msg='subStepMinCryst')
  if (num%subStepSizeCryst  <= 0.0_pReal)      call IO_error(301,ext_msg='subStepSizeCryst')
  if (num%stepIncreaseCryst <= 0.0_pReal)      call IO_error(301,ext_msg='stepIncreaseCryst')

  if (num%subStepSizeLp <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLp')
  if (num%subStepSizeLi <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLi')

  if (num%rtol_crystalliteState  <= 0.0_pReal) call IO_error(301,ext_msg='rtol_crystalliteState')
  if (num%rtol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='rtol_crystalliteStress')
  if (num%atol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='atol_crystalliteStress')

  if (num%iJacoLpresiduum < 1)                 call IO_error(301,ext_msg='iJacoLpresiduum')

  if (num%nState < 1)                          call IO_error(301,ext_msg='nState')
  if (num%nStress< 1)                          call IO_error(301,ext_msg='nStress')


  phases => config_material%get('phase')

  !$OMP PARALLEL DO PRIVATE(ce,ph,en)
  do el = 1, discretization_Nelems
    do ip = 1, discretization_nIPs
      ce = (el-1)*discretization_nIPs + ip
      do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
        en = material_phaseEntry(co,ce)
        ph = material_phaseID(co,ce)
        call crystallite_orientations(co,ip,el)
        call plastic_dependentState(ph,en)                                                          ! update dependent state variables to be consistent with basic states
     end do
    end do
  end do
  !$OMP END PARALLEL DO


end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations(co,ip,el)

  integer, intent(in) :: &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el                                                                                              !< counter in element loop

  integer :: ph, en


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  call phase_O(ph)%data(en)%fromMatrix(transpose(math_rotationalPart(mechanical_F_e(ph,en))))

  if (plasticState(material_phaseID(1,(el-1)*discretization_nIPs + ip))%nonlocal) &
    call plastic_nonlocal_updateCompatibility(phase_O,material_phaseID(1,(el-1)*discretization_nIPs + ip),ip,el)


end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(co,ce, tensor33)

  real(pReal), dimension(3,3), intent(in) :: tensor33
  integer, intent(in):: &
    co, &
    ce
  real(pReal), dimension(3,3) :: crystallite_push33ToRef

  real(pReal), dimension(3,3)             :: T
  integer :: ph, en

  ph = material_phaseID(co,ce)
  en = material_phaseEntry(co,ce)
  T = matmul(phase_O_0(ph)%data(en)%asMatrix(),transpose(math_inv33(phase_F(co,ce))))               ! ToDo: initial orientation correct?

  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


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
!> @brief Write restart data to file.
!--------------------------------------------------------------------------------------------------
subroutine phase_restartWrite(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ph


  groupHandle(1) = HDF5_addGroup(fileHandle,'phase')

  do ph = 1, size(material_name_phase)

    groupHandle(2) = HDF5_addGroup(groupHandle(1),material_name_phase(ph))

    call mechanical_restartWrite(groupHandle(2),ph)
    call thermal_restartWrite(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  end do

  call HDF5_closeGroup(groupHandle(1))

end subroutine phase_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read restart data from file.
!--------------------------------------------------------------------------------------------------
subroutine phase_restartRead(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ph


  groupHandle(1) = HDF5_openGroup(fileHandle,'phase')

  do ph = 1, size(material_name_phase)

    groupHandle(2) = HDF5_openGroup(groupHandle(1),material_name_phase(ph))

    call mechanical_restartRead(groupHandle(2),ph)
    call thermal_restartRead(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  end do

  call HDF5_closeGroup(groupHandle(1))

end subroutine phase_restartRead


end module phase
