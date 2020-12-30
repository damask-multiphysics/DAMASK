!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all thermal sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(constitutive) constitutive_thermal

  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: T
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: current

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


  module subroutine source_thermal_dissipation_getRateAndItsTangent(TDot, dTDot_dT, Tstar, Lp, phase)
    integer, intent(in) :: &
      phase                                                                                         !< phase ID of element
    real(pReal),  intent(in), dimension(3,3) :: &
      Tstar                                                                                         !< 2nd Piola Kirchhoff stress tensor for a given element
    real(pReal),  intent(in), dimension(3,3) :: &
      Lp                                                                                            !< plastic velocuty gradient for a given element
    real(pReal),  intent(out) :: &
      TDot, &
      dTDot_dT
  end subroutine source_thermal_dissipation_getRateAndItsTangent

  module subroutine source_thermal_externalheat_getRateAndItsTangent(TDot, dTDot_dT, phase, of)
    integer, intent(in) :: &
      phase, &
      of
    real(pReal),  intent(out) :: &
      TDot, &
      dTDot_dT
  end subroutine source_thermal_externalheat_getRateAndItsTangent

 end interface

contains

!----------------------------------------------------------------------------------------------
!< @brief initializes thermal sources and kinematics mechanism
!----------------------------------------------------------------------------------------------
module subroutine thermal_init(phases)

  class(tNode), pointer :: &
    phases

  integer :: &
    ph, &
    Nconstituents


  print'(/,a)', ' <<<+-  constitutive_mech init  -+>>>'

  allocate(current(phases%length))


  do ph = 1, phases%length

    Nconstituents = count(material_phaseAt == ph) * discretization_nIPs

    allocate(current(ph)%T(Nconstituents))

  enddo

! initialize source mechanisms
  if(maxval(phase_Nsources) /= 0) then
    where(source_thermal_dissipation_init (maxval(phase_Nsources))) phase_source = SOURCE_thermal_dissipation_ID
    where(source_thermal_externalheat_init(maxval(phase_Nsources))) phase_source = SOURCE_thermal_externalheat_ID
  endif

!--------------------------------------------------------------------------------------------------
!initialize kinematic mechanisms
  if(maxval(phase_Nkinematics) /= 0) where(kinematics_thermal_expansion_init(maxval(phase_Nkinematics))) &
                                           phase_kinematics = KINEMATICS_thermal_expansion_ID

end subroutine thermal_init


!----------------------------------------------------------------------------------------------
!< @brief calculates thermal dissipation rate
!----------------------------------------------------------------------------------------------
module subroutine constitutive_thermal_getRateAndItsTangents(TDot, dTDot_dT, T, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(in) :: &
    T                                                                                             !< plastic velocity gradient
  real(pReal), intent(inout) :: &
    TDot, &
    dTDot_dT

  real(pReal) :: &
    my_Tdot, &
    my_dTdot_dT
  integer :: &
    ph, &
    homog, &
    instance, &
    me, &
    so, &
    co

   homog  = material_homogenizationAt(el)
   instance = thermal_typeInstance(homog)

  do co = 1, homogenization_Nconstituents(homog)
     ph = material_phaseAt(co,el)
     me = material_phasememberAt(co,ip,el)
     do so = 1, phase_Nsources(ph)
       select case(phase_source(so,ph))
         case (SOURCE_thermal_dissipation_ID)
          call source_thermal_dissipation_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                               mech_S(ph,me),mech_L_p(ph,me), ph)

         case (SOURCE_thermal_externalheat_ID)
          call source_thermal_externalheat_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                                ph, me)

         case default
          my_Tdot = 0.0_pReal
          my_dTdot_dT = 0.0_pReal
       end select
       Tdot = Tdot + my_Tdot
       dTdot_dT = dTdot_dT + my_dTdot_dT
     enddo
   enddo

end subroutine constitutive_thermal_getRateAndItsTangents


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
module subroutine constitutive_thermal_setT(T,co,ip,el)

  real(pReal), intent(in) :: T
  integer, intent(in) :: co, ip, el


  current(material_phaseAt(co,el))%T(material_phaseMemberAt(co,ip,el)) = T

end subroutine constitutive_thermal_setT


end submodule constitutive_thermal
