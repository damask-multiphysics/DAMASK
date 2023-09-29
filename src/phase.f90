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
  use tables
  use IO
  use config
  use material
  use result
  use crystal
  use discretization
  use parallelization
  use HDF5
  use HDF5_utilities

  implicit none(type,external)
  private

  type :: tState
    integer :: &
      sizeState        = 0, &                                                                       !< size of state
      sizeDotState     = 0, &                                                                       !< size of dot state, i.e. state(1:sizeDot) follows time evolution by dotState rates
      offsetDeltaState = 0, &                                                                       !< index offset of delta state
      sizeDeltaState   = 0                                                                          !< size of delta state, i.e. state(offset+1:offset+sizeDelta) follows time evolution by deltaState increments
    real(pREAL), allocatable, dimension(:) :: &
      atol
    ! http://stackoverflow.com/questions/3948210
    real(pREAL), pointer,     dimension(:,:), contiguous :: &                                       !< is basically an allocatable+target, but in a type needs to be pointer
      state0, &
      state, &                                                                                      !< state
      dotState, &                                                                                   !< rate of state change
      deltaState                                                                                    !< increment of state change
    real(pREAL), pointer,     dimension(:,:)  :: &
      deltaState2
  end type

  type, extends(tState) :: tPlasticState
    logical :: nonlocal = .false.
  end type

  type :: tSourceState
    type(tState), dimension(:), allocatable :: p                                                    !< tState for each active source mechanism in a phase
  end type

  enum, bind(c); enumerator :: &
    UNDEFINED, &
    MECHANICAL_PLASTICITY_NONE, &
    MECHANICAL_PLASTICITY_ISOTROPIC, &
    MECHANICAL_PLASTICITY_PHENOPOWERLAW, &
    MECHANICAL_PLASTICITY_KINEHARDENING, &
    MECHANICAL_PLASTICITY_DISLOTWIN, &
    MECHANICAL_PLASTICITY_DISLOTUNGSTEN, &
    MECHANICAL_PLASTICITY_NONLOCAL, &
    MECHANICAL_EIGEN_THERMALEXPANSION, &
    DAMAGE_ISOBRITTLE, &
    DAMAGE_ANISOBRITTLE, &
    THERMAL_SOURCE_DISSIPATION, &
    THERMAL_SOURCE_EXTERNALHEAT
  end enum


  integer(kind(UNDEFINED)), dimension(:), allocatable :: &
    mechanical_plasticity_type, &                                                                   !< plasticity of each phase
    damage_type                                                                                     !< damage type of each phase
  integer(kind(UNDEFINED)),  dimension(:,:), allocatable :: &
    thermal_source_type, &
    mechanical_eigen_kinematics_type

  character(len=2), allocatable, dimension(:) :: phase_lattice
  real(pREAL),      allocatable, dimension(:) :: phase_cOverA
  real(pREAL),      allocatable, dimension(:) :: phase_rho

  type(tRotationContainer), dimension(:), allocatable :: &
    phase_O_0, &
    phase_O

  type :: tNumerics
    integer :: &
      iJacoLpresiduum, &                                                                            !< frequency of Jacobian update of residuum in Lp
      iJacoLiresiduum, &                                                                            !< frequency of Jacobian update of residuum in Li
      nState, &                                                                                     !< state loop limit
      nStress_Lp, &                                                                                 !< stress loop limit for Lp
      nStress_Li                                                                                    !< stress loop limit for Li
    real(pREAL) :: &
      stepMinCryst, &                                                                               !< minimum (relative) size of sub-step allowed during cutback
      stepSizeCryst, &                                                                              !< size of first substep when cutback
      stepSizeLp, &                                                                                 !< size of first substep when cutback in Lp calculation
      stepSizeLi, &                                                                                 !< size of first substep when cutback in Li calculation
      stepIncreaseCryst, &                                                                          !< increase of next substep size when previous substep converged
      rtol_crystalliteState, &
      rtol_Lp, &                                                                                    !< relative tolerance in stress loop for Lp
      atol_Lp, &                                                                                    !< absolute tolerance in stress loop for Lp
      rtol_Li, &                                                                                    !< relative tolerance in stress loop for Li
      atol_Li                                                                                       !< absolute tolerance in stress loop for Li
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

  type(tPlasticState), allocatable, dimension(:), public :: &
    plasticState
  type(tState),  allocatable, dimension(:), public :: &
    damageState


  interface

! == cleaned:begin =================================================================================
    module subroutine mechanical_init(phases,num_mech)
      type(tDict), pointer :: phases, num_mech
    end subroutine mechanical_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine thermal_init(phases)
      type(tDict), pointer :: phases
    end subroutine thermal_init


    module subroutine mechanical_result(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine mechanical_result

    module subroutine damage_result(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine damage_result

    module subroutine thermal_result(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine thermal_result

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
      real(pREAL), intent(in) :: Delta_t
      integer, intent(in) :: &
        co, &                                                                                       !< counter in constituent loop
        ce
      real(pREAL), dimension(3,3,3,3) :: dPdF
    end function phase_mechanical_dPdF

    module subroutine mechanical_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mechanical_restartWrite

    module subroutine thermal_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine thermal_restartWrite

    module subroutine damage_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine damage_restartWrite

    module subroutine mechanical_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mechanical_restartRead

    module subroutine thermal_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine thermal_restartRead

    module subroutine damage_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine damage_restartRead

    module function mechanical_S(ph,en) result(S)
      integer, intent(in) :: ph,en
      real(pREAL), dimension(3,3) :: S
    end function mechanical_S

    module function mechanical_L_p(ph,en) result(L_p)
      integer, intent(in) :: ph,en
      real(pREAL), dimension(3,3) :: L_p
    end function mechanical_L_p

    module function mechanical_F_e(ph,en) result(F_e)
      integer, intent(in) :: ph,en
      real(pREAL), dimension(3,3) :: F_e
    end function mechanical_F_e

    module function mechanical_F_i(ph,en) result(F_i)
      integer, intent(in) :: ph,en
      real(pREAL), dimension(3,3) :: F_i
    end function mechanical_F_i

    module function phase_F(co,ce) result(F)
      integer, intent(in) :: co, ce
      real(pREAL), dimension(3,3) :: F
    end function phase_F

    module function phase_P(co,ce) result(P)
      integer, intent(in) :: co, ce
      real(pREAL), dimension(3,3) :: P
    end function phase_P

    pure module function thermal_T(ph,en) result(T)
      integer, intent(in) :: ph,en
      real(pREAL) :: T
    end function thermal_T

    module function thermal_dot_T(ph,en) result(dot_T)
      integer, intent(in) :: ph,en
      real(pREAL) :: dot_T
    end function thermal_dot_T

    module function damage_phi(ph,en) result(phi)
      integer, intent(in) :: ph,en
      real(pREAL) :: phi
    end function damage_phi


    module subroutine phase_set_F(F,co,ce)
      real(pREAL), dimension(3,3), intent(in) :: F
      integer, intent(in) :: co, ce
    end subroutine phase_set_F

    module subroutine phase_thermal_setField(T,dot_T, co,ce)
      real(pREAL), intent(in) :: T, dot_T
      integer, intent(in) :: co, ce
    end subroutine phase_thermal_setField

    module subroutine phase_set_phi(phi,co,ce)
      real(pREAL), intent(in) :: phi
      integer, intent(in) :: co, ce
    end subroutine phase_set_phi


    module function phase_mu_phi(co,ce) result(mu)
      integer, intent(in) :: co, ce
      real(pREAL) :: mu
    end function phase_mu_phi

    module function phase_K_phi(co,ce) result(K)
      integer, intent(in) :: co, ce
      real(pREAL), dimension(3,3) :: K
    end function phase_K_phi


    module function phase_mu_T(co,ce) result(mu)
      integer, intent(in) :: co, ce
      real(pREAL) :: mu
    end function phase_mu_T

    module function phase_K_T(co,ce) result(K)
      integer, intent(in) :: co, ce
      real(pREAL), dimension(3,3) :: K
    end function phase_K_T

! == cleaned:end ===================================================================================

    module function phase_thermal_constitutive(Delta_t,ph,en) result(converged_)

      real(pREAL), intent(in) :: Delta_t
      integer, intent(in) :: ph, en
      logical :: converged_

    end function phase_thermal_constitutive

    module function phase_damage_constitutive(Delta_t,co,ce) result(converged_)
      real(pREAL), intent(in) :: Delta_t
      integer, intent(in) :: co, ce
      logical :: converged_
    end function phase_damage_constitutive

    module function phase_mechanical_constitutive(Delta_t,co,ce) result(converged_)
      real(pREAL), intent(in) :: Delta_t
      integer, intent(in) :: co, ce
      logical :: converged_
    end function phase_mechanical_constitutive

    !ToDo: Merge all the stiffness functions
    module function phase_homogenizedC66(ph,en) result(C)
      integer, intent(in) :: ph, en
      real(pREAL), dimension(6,6) :: C
    end function phase_homogenizedC66
    module function phase_damage_C66(C66,ph,en) result(C66_degraded)
      real(pREAL), dimension(6,6), intent(in)  :: C66
      integer,                     intent(in)  :: ph,en
      real(pREAL), dimension(6,6) :: C66_degraded
    end function phase_damage_C66

    module function phase_f_phi(phi,co,ce) result(f)
      integer, intent(in) :: ce,co
      real(pREAL), intent(in) :: &
        phi                                                                                         !< damage parameter
      real(pREAL) :: &
        f
    end function phase_f_phi

    module function phase_f_T(ph,en) result(f)
      integer, intent(in) :: ph, en
      real(pREAL) :: f
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


    module subroutine damage_anisobrittle_LiAndItsTangent(L_i, dL_i_dM_i, M_i, ph,en)
      integer, intent(in) :: ph, en
      real(pREAL),   intent(in),  dimension(3,3) :: &
        M_i
      real(pREAL),   intent(out), dimension(3,3) :: &
        L_i                                                                                         !< damage velocity gradient
      real(pREAL),   intent(out), dimension(3,3,3,3) :: &
        dL_i_dM_i                                                                                   !< derivative of L_i with respect to M_i
    end subroutine damage_anisobrittle_LiAndItsTangent

  end interface


#if __INTEL_COMPILER >= 1900
  public :: &
    prec, &
    math, &
    rotations, &
    IO, &
    config, &
    material, &
    result, &
    crystal, &
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
    phase_result, &
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
  type(tDict), pointer :: &
    phases, &
    phase, &
    num_phase, &
    num_mech
  character(len=:), allocatable :: refs


  print'(/,1x,a)', '<<<+-  phase init  -+>>>'; flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(phase_lattice(phases%length))
  allocate(phase_cOverA(phases%length),source=-1.0_pREAL)
  allocate(phase_rho(phases%length))
  allocate(phase_O_0(phases%length))

  do ph = 1,phases%length
    print'(/,1x,a,i0,a)', 'phase ',ph,': '//phases%key(ph)
    phase => phases%get_dict(ph)
    refs = config_listReferences(phase,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs
    phase_rho(ph) = phase%get_asReal('rho',defaultVal=0.0_pREAL)
    phase_lattice(ph) = phase%get_asStr('lattice')
    if (all(phase_lattice(ph) /= ['cF','cI','hP','tI'])) &
      call IO_error(130,ext_msg='phase_init: '//phase%get_asStr('lattice'))
    if (any(phase_lattice(ph) == ['hP','tI'])) &
      phase_cOverA(ph) = phase%get_asReal('c/a')
    allocate(phase_O_0(ph)%data(count(material_ID_phase==ph)))
  end do

  do ce = 1, size(material_ID_phase,2)
    ma = discretization_materialAt((ce-1)/discretization_nIPs+1)
    do co = 1,homogenization_Nconstituents(material_ID_homogenization(ce))
      ph = material_ID_phase(co,ce)
      phase_O_0(ph)%data(material_entry_phase(co,ce)) = material_O_0(ma)%data(co)
    end do
  end do

  allocate(phase_O(phases%length))
  do ph = 1,phases%length
    phase_O(ph)%data = phase_O_0(ph)%data
  end do

  num_phase => config_numerics%get_dict('phase',defaultVal=emptyDict)
  num_mech  => num_phase%get_dict('mechanical', defaultVal=emptyDict)

  call mechanical_init(phases,num_mech)
  call damage_init()
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

  allocate(state%atol             (sizeState),          source=0.0_pREAL)
  allocate(state%state0           (sizeState,NEntries), source=0.0_pREAL)
  allocate(state%state            (sizeState,NEntries), source=0.0_pREAL)

  allocate(state%dotState      (sizeDotState,NEntries), source=0.0_pREAL)

  allocate(state%deltaState  (sizeDeltaState,NEntries), source=0.0_pREAL)
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
subroutine phase_result()

  integer :: ph
  character(len=:), allocatable :: group


  call result_closeGroup(result_addGroup('/current/phase/'))

  do ph = 1, size(material_name_phase)

    group = '/current/phase/'//trim(material_name_phase(ph))//'/'
    call result_closeGroup(result_addGroup(group))

    call mechanical_result(group,ph)
    call damage_result(group,ph)
    call thermal_result(group,ph)

  end do

end subroutine phase_result


!--------------------------------------------------------------------------------------------------
!> @brief Allocate and initialize.
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init()

  integer :: &
    ce, &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el, &                                                                                           !< counter in element loop
    en, ph
  type(tDict), pointer :: &
    num_phase, &
    phases

  phases => config_material%get_dict('phase')

  !$OMP PARALLEL DO PRIVATE(ce,ph,en)
  do el = 1, discretization_Nelems
    do ip = 1, discretization_nIPs
      ce = (el-1)*discretization_nIPs + ip
      do co = 1,homogenization_Nconstituents(material_ID_homogenization(ce))
        en = material_entry_phase(co,ce)
        ph = material_ID_phase(co,ce)
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


  ph = material_ID_phase(co,(el-1)*discretization_nIPs + ip)
  en = material_entry_phase(co,(el-1)*discretization_nIPs + ip)

  call phase_O(ph)%data(en)%fromMatrix(transpose(math_rotationalPart(mechanical_F_e(ph,en))))

  if (plasticState(material_ID_phase(1,(el-1)*discretization_nIPs + ip))%nonlocal) &
    call plastic_nonlocal_updateCompatibility(phase_O,material_ID_phase(1,(el-1)*discretization_nIPs + ip),ip,el)


end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(co,ce, tensor33)

  real(pREAL), dimension(3,3), intent(in) :: tensor33
  integer, intent(in):: &
    co, &
    ce
  real(pREAL), dimension(3,3) :: crystallite_push33ToRef

  real(pREAL), dimension(3,3) :: T
  integer :: ph, en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)
  T = matmul(phase_O_0(ph)%data(en)%asMatrix(),transpose(math_inv33(phase_F(co,ce))))               ! ToDo: initial orientation correct?

  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief determines whether a point is converged
!--------------------------------------------------------------------------------------------------
logical pure function converged(residuum,state,atol)

  real(pREAL), intent(in), dimension(:) ::&
    residuum, state, atol
  real(pREAL) :: &
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
    call damage_restartWrite(groupHandle(2),ph)

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
    call damage_restartRead(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  end do

  call HDF5_closeGroup(groupHandle(1))

end subroutine phase_restartRead


end module phase
