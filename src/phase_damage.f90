!--------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all damage sources and kinematics constitutive models
!--------------------------------------------------------------------------------------------------
submodule(phase) damage

  type :: tDamageParameters
    real(pREAL) :: &
      mu = 0.0_pREAL, &                                                                             !< viscosity
      l_c = 0.0_pREAL                                                                               !< characteristic length
  end type tDamageParameters

  integer :: phase_damage_maxSizeDotState

  type :: tFieldQuantities
    real(pREAL), dimension(:), allocatable :: phi
  end type tFieldQuantities


  type(tFieldQuantities), dimension(:), allocatable :: current

  type(tDamageParameters), dimension(:), allocatable :: param

  interface

    module function anisobrittle_init() result(mySources)
      logical, dimension(:), allocatable :: mySources
    end function anisobrittle_init

    module function isobrittle_init() result(mySources)
      logical, dimension(:), allocatable :: mySources
    end function isobrittle_init


    module subroutine isobrittle_deltaState(C, Fe, ph, en)
      integer, intent(in) :: ph,en
      real(pREAL),  intent(in), dimension(3,3) :: &
        Fe
      real(pREAL),  intent(in), dimension(6,6) :: &
        C
    end subroutine isobrittle_deltaState


    module subroutine anisobrittle_dotState(M_i, ph, en)
      integer, intent(in) :: ph,en
      real(pREAL),  intent(in), dimension(3,3) :: &
        M_i
    end subroutine anisobrittle_dotState


    module subroutine anisobrittle_result(phase,group)
      integer,          intent(in) :: phase
      character(len=*), intent(in) :: group
    end subroutine anisobrittle_result

    module subroutine isobrittle_result(phase,group)
      integer,          intent(in) :: phase
      character(len=*), intent(in) :: group
    end subroutine isobrittle_result

 end interface

contains

!----------------------------------------------------------------------------------------------
!< @brief Initialize damage mechanisms.
!----------------------------------------------------------------------------------------------
module subroutine damage_init()

  integer :: &
    ph, &
    Nmembers
  type(tDict), pointer :: &
    phases, &
    phase, &
    source
  character(len=:), allocatable :: refs
  logical:: damage_active


  print'(/,1x,a)', '<<<+-  phase:damage init  -+>>>'


  phases => config_material%get_dict('phase')
  allocate(current(phases%length))
  allocate(damageState(phases%length))
  allocate(param(phases%length))

  damage_active = .false.
  do ph = 1,phases%length

    Nmembers = count(material_ID_phase == ph)

    allocate(current(ph)%phi(Nmembers),source=1.0_pREAL)

    phase => phases%get_dict(ph)
    source => phase%get_dict('damage',defaultVal=emptyDict)
    if (source%length > 0) then
      print'(/,1x,a,i0,a)', 'phase ',ph,': '//phases%key(ph)
      refs = config_listReferences(source,indent=3)
      if (len(refs) > 0) print'(/,1x,a)', refs
      damage_active = .true.
      param(ph)%mu = source%get_asReal('mu')
      param(ph)%l_c = source%get_asReal('l_c')
    end if

  end do

  allocate(damage_type(phases%length), source = UNDEFINED)

  if (damage_active) then
    where(isobrittle_init()  ) damage_type = DAMAGE_ISOBRITTLE
    where(anisobrittle_init()) damage_type = DAMAGE_ANISOBRITTLE
  end if

  phase_damage_maxSizeDotState = maxval(damageState%sizeDotState)

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P)
!--------------------------------------------------------------------------------------------------
module function phase_damage_constitutive(Delta_t,co,ce) result(converged_)

  real(pREAL), intent(in) :: Delta_t
  integer, intent(in) :: &
    co, &
    ce
  logical :: converged_

  integer :: &
    ph, en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)

  converged_ = .not. integrateDamageState(Delta_t,ph,en)

end function phase_damage_constitutive


!--------------------------------------------------------------------------------------------------
!> @brief returns the degraded/modified elasticity matrix
!--------------------------------------------------------------------------------------------------
module function phase_damage_C66(C66,ph,en) result(C66_degraded)

  real(pREAL), dimension(6,6), intent(in)  :: C66
  integer,                     intent(in)  :: ph,en
  real(pREAL), dimension(6,6) :: C66_degraded


  damageType: select case (damage_type(ph))
    case (DAMAGE_ISOBRITTLE) damageType
      C66_degraded = C66 * damage_phi(ph,en)**2
    case default damageType
      C66_degraded = C66
  end select damageType

end function phase_damage_C66


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
module subroutine damage_restore(ce)

  integer, intent(in) :: ce

  integer :: &
    co


  do co = 1,homogenization_Nconstituents(material_ID_homogenization(ce))
    if (damageState(material_ID_phase(co,ce))%sizeState > 0) &
    damageState(material_ID_phase(co,ce))%state( :,material_entry_phase(co,ce)) = &
      damageState(material_ID_phase(co,ce))%state0(:,material_entry_phase(co,ce))
  end do

end subroutine damage_restore


!----------------------------------------------------------------------------------------------
!< @brief returns local part of nonlocal damage driving force
!----------------------------------------------------------------------------------------------
module function phase_f_phi(phi,co,ce) result(f)

  integer, intent(in) :: ce,co
  real(pREAL), intent(in) :: &
    phi                                                                                             !< damage parameter
  real(pREAL) :: &
    f

  integer :: &
    ph, &
    en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)

  select case(damage_type(ph))
    case(DAMAGE_ISOBRITTLE,DAMAGE_ANISOBRITTLE)
      f = 1.0_pREAL &
        - 2.0_pREAL * phi*damageState(ph)%state(1,en)                                               ! ToDo: MD: seems to be phi**2
    case default
      f = 0.0_pREAL
  end select

end function phase_f_phi


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
function integrateDamageState(Delta_t,ph,en) result(broken)

  real(pREAL), intent(in) :: Delta_t
  integer, intent(in) :: &
    ph, &
    en
  logical :: broken

  integer :: &
    NiterationState, &                                                                              !< number of iterations in state loop
    size_so
  real(pREAL) :: &
    zeta
  real(pREAL), dimension(phase_damage_maxSizeDotState) :: &
    r                                                                                               ! state residuum
  real(pREAL), dimension(phase_damage_maxSizeDotState,2) :: source_dotState
  logical :: &
    converged_


  if (damageState(ph)%sizeState == 0) then
    broken = .false.
    return
  end if

  converged_ = .true.
  broken = phase_damage_collectDotState(ph,en)
  if (broken) return

  size_so = damageState(ph)%sizeDotState
  damageState(ph)%state(1:size_so,en) = damageState(ph)%state0  (1:size_so,en) &
                                      + damageState(ph)%dotState(1:size_so,en) * Delta_t
  source_dotState(1:size_so,2) = 0.0_pREAL

  iteration: do NiterationState = 1, num%nState

    if (nIterationState > 1) source_dotState(1:size_so,2) = source_dotState(1:size_so,1)
    source_dotState(1:size_so,1) = damageState(ph)%dotState(:,en)

    broken = phase_damage_collectDotState(ph,en)
    if (broken) exit iteration


      zeta = damper(damageState(ph)%dotState(:,en),source_dotState(1:size_so,1),source_dotState(1:size_so,2))
      damageState(ph)%dotState(:,en) = damageState(ph)%dotState(:,en) * zeta &
                                     + source_dotState(1:size_so,1)* (1.0_pREAL - zeta)
      r(1:size_so) = damageState(ph)%state   (1:size_so,en)  &
                   - damageState(ph)%State0  (1:size_so,en)  &
                   - damageState(ph)%dotState(1:size_so,en) * Delta_t
      damageState(ph)%state(1:size_so,en) = damageState(ph)%state(1:size_so,en) - r(1:size_so)
      converged_ = converged_  .and. converged(r(1:size_so), &
                                               damageState(ph)%state(1:size_so,en), &
                                               damageState(ph)%atol(1:size_so))


    if (converged_) then
      broken = phase_damage_deltaState(mechanical_F_e(ph,en),ph,en)
      exit iteration
    end if

  end do iteration

  broken = broken .or. .not. converged_


  contains
  !--------------------------------------------------------------------------------------------------
  !> @brief Calculate the damping for correction of state and dot state.
  !--------------------------------------------------------------------------------------------------
  real(pREAL) pure function damper(omega_0,omega_1,omega_2)

  real(pREAL), dimension(:), intent(in) :: &
    omega_0, omega_1, omega_2

  real(pREAL) :: dot_prod12, dot_prod22

  dot_prod12 = dot_product(omega_0-omega_1, omega_1-omega_2)
  dot_prod22 = dot_product(omega_1-omega_2, omega_1-omega_2)

  if (min(dot_product(omega_0,omega_1),dot_prod12) < 0.0_pREAL .and. dot_prod22 > 0.0_pREAL) then
    damper = 0.75_pREAL + 0.25_pREAL * tanh(2.0_pREAL + 4.0_pREAL * dot_prod12 / dot_prod22)
  else
    damper = 1.0_pREAL
  end if

  end function damper

end function integrateDamageState


module subroutine damage_restartWrite(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer, intent(in) :: ph


  select case(damage_type(ph))
    case(DAMAGE_ISOBRITTLE,DAMAGE_ANISOBRITTLE)
      call HDF5_write(damageState(ph)%state,groupHandle,'omega_damage')
  end select

end subroutine damage_restartWrite


module subroutine damage_restartRead(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer, intent(in) :: ph


  select case(damage_type(ph))
    case(DAMAGE_ISOBRITTLE,DAMAGE_ANISOBRITTLE)
  call HDF5_read(damageState(ph)%state0,groupHandle,'omega_damage')
  end select


end subroutine damage_restartRead


!----------------------------------------------------------------------------------------------
!< @brief writes damage sources results to HDF5 output file
!----------------------------------------------------------------------------------------------
module subroutine damage_result(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph


  if (damage_type(ph) /= UNDEFINED) &
    call result_closeGroup(result_addGroup(group//'damage'))

  sourceType: select case (damage_type(ph))

    case (DAMAGE_ISOBRITTLE) sourceType
      call isobrittle_result(ph,group//'damage/')

    case (DAMAGE_ANISOBRITTLE) sourceType
      call anisobrittle_result(ph,group//'damage/')

  end select sourceType

end subroutine damage_result


!--------------------------------------------------------------------------------------------------
!> @brief Constitutive equation for calculating the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
function phase_damage_collectDotState(ph,en) result(broken)

  integer, intent(in) :: &
    ph, &
    en                                                                                         !< counter in source loop
  logical :: broken


  broken = .false.

  if (damageState(ph)%sizeState > 0) then

    sourceType: select case (damage_type(ph))

      case (DAMAGE_ANISOBRITTLE) sourceType
        call anisobrittle_dotState(mechanical_S(ph,en), ph,en) ! ToDo: use M_d

    end select sourceType

    broken = broken .or. any(IEEE_is_NaN(damageState(ph)%dotState(:,en)))

  end if

end function phase_damage_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief Damage viscosity.
!--------------------------------------------------------------------------------------------------
module function phase_mu_phi(co,ce) result(mu)

  integer, intent(in) :: co, ce
  real(pREAL) :: mu


  mu = param(material_ID_phase(co,ce))%mu

end function phase_mu_phi


!--------------------------------------------------------------------------------------------------
!> @brief Damage conductivity/diffusivity in reference configuration.
!--------------------------------------------------------------------------------------------------
module function phase_K_phi(co,ce) result(K)

  integer, intent(in) :: co, ce
  real(pREAL), dimension(3,3) :: K


  K = crystallite_push33ToRef(co,ce,param(material_ID_phase(co,ce))%l_c**2*math_I3)

end function phase_K_phi


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
function phase_damage_deltaState(Fe, ph, en) result(broken)

  integer, intent(in) :: &
    ph, &
    en
  real(pREAL),   intent(in), dimension(3,3) :: &
    Fe                                                                                              !< elastic deformation gradient

  integer :: &
    myOffset, &
    mySize
  logical :: &
    broken


  broken = .false.

  if (damageState(ph)%sizeState == 0) return

   sourceType: select case (damage_type(ph))

    case (DAMAGE_ISOBRITTLE) sourceType
      call isobrittle_deltaState(phase_homogenizedC66(ph,en), Fe, ph,en)
      broken = any(IEEE_is_NaN(damageState(ph)%deltaState(:,en)))
      if (.not. broken) then
        myOffset = damageState(ph)%offsetDeltaState
        mySize   = damageState(ph)%sizeDeltaState
        damageState(ph)%state(myOffset + 1: myOffset + mySize,en) = &
        damageState(ph)%state(myOffset + 1: myOffset + mySize,en) + damageState(ph)%deltaState(1:mySize,en)
      end if

  end select sourceType


end function phase_damage_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Check if a source mechanism is active or not.
!--------------------------------------------------------------------------------------------------
function source_active(source_label)  result(active_source)

  character(len=*), intent(in)         :: source_label                                              !< name of source mechanism
  logical, dimension(:), allocatable  :: active_source

  type(tDict), pointer :: &
    phases, &
    phase, &
    src
  integer :: ph


  phases => config_material%get_dict('phase')
  allocate(active_source(phases%length))
  do ph = 1, phases%length
    phase => phases%get_dict(ph)
    src => phase%get_dict('damage',defaultVal=emptyDict)
    active_source(ph) = src%get_asStr('type',defaultVal = 'x') == source_label
  end do


end function source_active


!----------------------------------------------------------------------------------------------
!< @brief Set damage parameter
!----------------------------------------------------------------------------------------------
module subroutine phase_set_phi(phi,co,ce)

  real(pREAL), intent(in) :: phi
  integer, intent(in) :: ce, co


  current(material_ID_phase(co,ce))%phi(material_entry_phase(co,ce)) = phi

end subroutine phase_set_phi


module function damage_phi(ph,en) result(phi)

  integer, intent(in) :: ph, en
  real(pREAL) :: phi


  phi = current(ph)%phi(en)

end function damage_phi



!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
!--------------------------------------------------------------------------------------------------
module subroutine damage_forward()

  integer :: ph


  do ph = 1, size(damageState)
    if (damageState(ph)%sizeState > 0) &
      damageState(ph)%state0 = damageState(ph)%state
  end do

end subroutine damage_forward


end submodule damage
