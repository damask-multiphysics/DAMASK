!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all thermal sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(phase) thermal

  enum, bind(c); enumerator :: &
    THERMAL_UNDEFINED_ID ,&
    THERMAL_DISSIPATION_ID, &
    THERMAL_EXTERNALHEAT_ID
  end enum

  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: T, dot_T
  end type tDataContainer
  integer(kind(THERMAL_UNDEFINED_ID)),  dimension(:,:), allocatable :: &
    thermal_source

  type(tDataContainer), dimension(:), allocatable :: current

  integer :: thermal_source_maxSizeDotState


  interface

    module function dissipation_init(source_length) result(mySources)
      integer, intent(in) :: source_length
      logical, dimension(:,:), allocatable :: mySources
    end function dissipation_init

    module function externalheat_init(source_length) result(mySources)
      integer, intent(in) :: source_length
      logical, dimension(:,:), allocatable :: mySources
    end function externalheat_init

    module function kinematics_thermal_expansion_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_thermal_expansion_init


    module subroutine externalheat_dotState(ph, me)
      integer, intent(in) :: &
        ph, &
        me
    end subroutine externalheat_dotState

    module subroutine dissipation_getRate(TDot, ph,me)
      integer, intent(in) :: &
        ph, &
        me
      real(pReal),  intent(out) :: &
        TDot
    end subroutine dissipation_getRate

    module subroutine externalheat_getRate(TDot, ph,me)
      integer, intent(in) :: &
        ph, &
        me
      real(pReal),  intent(out) :: &
        TDot
    end subroutine externalheat_getRate

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

    Nconstituents = count(material_phaseAt2 == ph)

    allocate(current(ph)%T(Nconstituents),source=300.0_pReal)
    allocate(current(ph)%dot_T(Nconstituents),source=0.0_pReal)
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
    where(dissipation_init (maxval(thermal_Nsources))) thermal_source = THERMAL_DISSIPATION_ID
    where(externalheat_init(maxval(thermal_Nsources))) thermal_source = THERMAL_EXTERNALHEAT_ID
  endif

  thermal_source_maxSizeDotState = 0
  PhaseLoop2:do ph = 1,phases%length

    do so = 1,thermal_Nsources(ph)
      thermalState(ph)%p(so)%state  = thermalState(ph)%p(so)%state0
    enddo

    thermal_source_maxSizeDotState  = max(thermal_source_maxSizeDotState, &
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
module subroutine constitutive_thermal_getRate(TDot, ph,me)

  integer, intent(in) :: ph, me
  real(pReal), intent(out) :: &
    TDot

  real(pReal) :: &
    my_Tdot
  integer :: &
    so


  TDot = 0.0_pReal

  do so = 1, thermal_Nsources(ph)
   select case(thermal_source(so,ph))
     case (THERMAL_DISSIPATION_ID)
      call dissipation_getRate(my_Tdot, ph,me)

     case (THERMAL_EXTERNALHEAT_ID)
      call externalheat_getRate(my_Tdot, ph,me)

     case default
      my_Tdot = 0.0_pReal
   end select
   Tdot = Tdot + my_Tdot
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
      call externalheat_dotState(ph,me)

    broken = broken .or. any(IEEE_is_NaN(thermalState(ph)%p(i)%dotState(:,me)))

  enddo SourceLoop

end function constitutive_thermal_collectDotState


module function thermal_stress(Delta_t,ph,me) result(converged_)

  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: ph, me
  logical :: converged_

  integer :: so


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
    thermalState(ph)%p(so)%state(1:sizeDotState,me) = thermalState(ph)%p(so)%state0(1:sizeDotState,me) &
                                                    + thermalState(ph)%p(so)%dotState(1:sizeDotState,me) * Delta_t
  enddo

end function integrateThermalState


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
module subroutine constitutive_thermal_setField(T,dot_T, co,ce)

  real(pReal), intent(in) :: T, dot_T
  integer, intent(in) :: ce, co


  current(material_phaseAt2(co,ce))%T(material_phaseMemberAt2(co,ce)) = T
  current(material_phaseAt2(co,ce))%dot_T(material_phaseMemberAt2(co,ce)) = dot_T

end subroutine constitutive_thermal_setField



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


end submodule thermal
