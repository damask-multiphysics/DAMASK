!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all thermal sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(constitutive) constitutive_thermal

  enum, bind(c); enumerator :: &
    THERMAL_UNDEFINED_ID ,&
    THERMAL_DISSIPATION_ID, &
    THERMAL_EXTERNALHEAT_ID
  end enum

  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: T
  end type tDataContainer
  integer(kind(THERMAL_UNDEFINED_ID)),     dimension(:,:), allocatable :: &
    thermal_source

  type(tDataContainer), dimension(:), allocatable :: current

  integer :: thermal_source_maxSizeDotState
  interface

  module function source_thermal_dissipation_init(source_length) result(mySources)
    integer, intent(in) :: source_length
    logical, dimension(:,:), allocatable :: mySources
  end function source_thermal_dissipation_init

  module function source_thermal_externalheat_init(source_length) result(mySources)
    integer, intent(in) :: source_length
    logical, dimension(:,:), allocatable :: mySources
  end function source_thermal_externalheat_init

  module function kinematics_thermal_expansion_init(kinematics_length) result(myKinematics)
    integer, intent(in) :: kinematics_length
    logical, dimension(:,:), allocatable :: myKinematics
  end function kinematics_thermal_expansion_init


    module subroutine source_thermal_externalheat_dotState(ph, me)
      integer, intent(in) :: &
        ph, &
        me
    end subroutine source_thermal_externalheat_dotState


  module subroutine thermal_dissipation_getRate(TDot, Tstar,Lp,phase)
    integer, intent(in) :: &
      phase                                                                                         !< phase ID of element
    real(pReal),  intent(in), dimension(3,3) :: &
      Tstar                                                                                         !< 2nd Piola Kirchhoff stress tensor for a given element
    real(pReal),  intent(in), dimension(3,3) :: &
      Lp                                                                                            !< plastic velocuty gradient for a given element
    real(pReal),  intent(out) :: &
      TDot
  end subroutine thermal_dissipation_getRate

  module subroutine thermal_externalheat_getRate(TDot, ph,me)
    integer, intent(in) :: &
      ph, &
      me
    real(pReal),  intent(out) :: &
      TDot
  end subroutine thermal_externalheat_getRate

 end interface

contains

!----------------------------------------------------------------------------------------------
!< @brief initializes thermal sources and kinematics mechanism
!----------------------------------------------------------------------------------------------
module subroutine thermal_init(phases)

  class(tNode), pointer :: &
    phases

  class(tNode), pointer :: &
    phase, thermal, sources

  integer :: &
    ph, so, &
    Nconstituents


  print'(/,a)', ' <<<+-  constitutive_thermal init  -+>>>'

  allocate(current(phases%length))

  allocate(thermalState (phases%length))
  allocate(thermal_Nsources(phases%length),source = 0)

  do ph = 1, phases%length

    Nconstituents = count(material_phaseAt == ph) * discretization_nIPs

    allocate(current(ph)%T(Nconstituents),source=300.0_pReal)
    phase => phases%get(ph)
    if(phase%contains('thermal')) then
      thermal => phase%get('thermal')
      sources => thermal%get('source',defaultVal=emptyList)

      thermal_Nsources(ph) = sources%length
    endif
    allocate(thermalstate(ph)%p(thermal_Nsources(ph)))
  enddo

  allocate(thermal_source(maxval(thermal_Nsources),phases%length), source = THERMAL_UNDEFINED_ID)

  if(maxval(thermal_Nsources) /= 0) then
    where(source_thermal_dissipation_init (maxval(thermal_Nsources))) thermal_source = THERMAL_DISSIPATION_ID
    where(source_thermal_externalheat_init(maxval(thermal_Nsources))) thermal_source = THERMAL_EXTERNALHEAT_ID
  endif

  thermal_source_maxSizeDotState = 0
  PhaseLoop2:do ph = 1,phases%length

    do so = 1,thermal_Nsources(ph)
      thermalState(ph)%p(so)%partitionedState0 = thermalState(ph)%p(so)%state0
      thermalState(ph)%p(so)%state             = thermalState(ph)%p(so)%partitionedState0
    enddo

    thermal_source_maxSizeDotState   = max(thermal_source_maxSizeDotState, &
                                                maxval(thermalState(ph)%p%sizeDotState))
  enddo PhaseLoop2

!--------------------------------------------------------------------------------------------------
!initialize kinematic mechanisms
  if(maxval(phase_Nkinematics) /= 0) where(kinematics_thermal_expansion_init(maxval(phase_Nkinematics))) &
                                           phase_kinematics = KINEMATICS_thermal_expansion_ID

end subroutine thermal_init


!----------------------------------------------------------------------------------------------
!< @brief calculates thermal dissipation rate
!----------------------------------------------------------------------------------------------
module subroutine constitutive_thermal_getRate(TDot, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(out) :: &
    TDot

  real(pReal) :: &
    my_Tdot
  integer :: &
    ph, &
    homog, &
    instance, &
    me, &
    so, &
    co

   homog  = material_homogenizationAt(el)
   instance = thermal_typeInstance(homog)

  TDot = 0.0_pReal
  do co = 1, homogenization_Nconstituents(homog)
     ph = material_phaseAt(co,el)
     me = material_phasememberAt(co,ip,el)
     do so = 1, thermal_Nsources(ph)
       select case(thermal_source(so,ph))
         case (THERMAL_DISSIPATION_ID)
          call thermal_dissipation_getRate(my_Tdot, mech_S(ph,me),mech_L_p(ph,me),ph)

         case (THERMAL_EXTERNALHEAT_ID)
          call thermal_externalheat_getRate(my_Tdot, ph,me)

         case default
          my_Tdot = 0.0_pReal
       end select
       Tdot = Tdot + my_Tdot
     enddo
   enddo

end subroutine constitutive_thermal_getRate


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_thermal_collectDotState(ph,me) result(broken)

  integer, intent(in) :: ph, me
  logical :: broken

  integer :: i


  broken = .false.

  SourceLoop: do i = 1, thermal_Nsources(ph)

    if (thermal_source(i,ph) == THERMAL_EXTERNALHEAT_ID) &
      call source_thermal_externalheat_dotState(ph,me)

    broken = broken .or. any(IEEE_is_NaN(thermalState(ph)%p(i)%dotState(:,me)))

  enddo SourceLoop

end function constitutive_thermal_collectDotState


module function thermal_stress(Delta_t,ph,me) result(converged_)

  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: ph, me
  logical :: converged_

  integer :: so


  do so = 1, thermal_Nsources(ph)
    thermalState(ph)%p(so)%state(:,me) = thermalState(ph)%p(so)%subState0(:,me)
  enddo
  converged_ = .not. integrateThermalState(Delta_t,ph,me)

end function thermal_stress


!--------------------------------------------------------------------------------------------------
!> @brief integrate state with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
function integrateThermalState(Delta_t, ph,me) result(broken)

  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: ph, me
  logical :: &
    broken

  integer :: &
    so, &
    sizeDotState

  broken = constitutive_thermal_collectDotState(ph,me)
  if(broken) return

  do so = 1, thermal_Nsources(ph)
    sizeDotState = thermalState(ph)%p(so)%sizeDotState
    thermalState(ph)%p(so)%state(1:sizeDotState,me) = thermalState(ph)%p(so)%subState0(1:sizeDotState,me) &
                                                    + thermalState(ph)%p(so)%dotState(1:sizeDotState,me) * Delta_t
  enddo

end function integrateThermalState


module subroutine constitutive_thermal_initializeRestorationPoints(ph,me)

  integer, intent(in) :: ph, me

  integer :: so


  do so = 1, size(thermalState(ph)%p)
    thermalState(ph)%p(so)%partitionedState0(:,me) = thermalState(ph)%p(so)%state0(:,me)
  enddo

end subroutine constitutive_thermal_initializeRestorationPoints



module subroutine thermal_forward()

  integer :: ph, so


  do ph = 1, size(thermalState)
    do so = 1, size(thermalState(ph)%p)
      thermalState(ph)%p(so)%state0 = thermalState(ph)%p(so)%state
    enddo
  enddo

end subroutine thermal_forward


!----------------------------------------------------------------------------------------------
!< @brief Get temperature (for use by non-thermal physics)
!----------------------------------------------------------------------------------------------
module function thermal_T(ph,me) result(T)

  integer, intent(in) :: ph, me
  real(pReal) :: T


  T = current(ph)%T(me)

end function thermal_T


!----------------------------------------------------------------------------------------------
!< @brief Set temperature
!----------------------------------------------------------------------------------------------
module subroutine constitutive_thermal_setT(T,co,ce)

  real(pReal), intent(in) :: T
  integer, intent(in) :: ce, co


  current(material_phaseAt2(co,ce))%T(material_phaseMemberAt2(co,ce)) = T

end subroutine constitutive_thermal_setT



!--------------------------------------------------------------------------------------------------
!> @brief checks if a source mechanism is active or not
!--------------------------------------------------------------------------------------------------
function thermal_active(source_label,src_length)  result(active_source)

  character(len=*), intent(in)         :: source_label                                              !< name of source mechanism
  integer,          intent(in)         :: src_length                                                !< max. number of sources in system
  logical, dimension(:,:), allocatable :: active_source

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, thermal, &
    src
  integer :: p,s

  phases => config_material%get('phase')
  allocate(active_source(src_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get(p)
    if (phase%contains('thermal')) then
      thermal =>  phase%get('thermal',defaultVal=emptyList)
      sources =>  thermal%get('source',defaultVal=emptyList)
      do s = 1, sources%length
        src => sources%get(s)
        if(src%get_asString('type') == source_label) active_source(s,p) = .true.
      enddo
    endif
  enddo


end function thermal_active


end submodule constitutive_thermal
