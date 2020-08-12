!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all thermal sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(constitutive) constitutive_thermal

  interface

  module subroutine source_thermal_dissipation_init
  end subroutine source_thermal_dissipation_init

  module subroutine source_thermal_externalheat_init
  end subroutine source_thermal_externalheat_init

  module subroutine kinematics_thermal_expansion_init
  end subroutine kinematics_thermal_expansion_init


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
module subroutine thermal_init

! initialize source mechanisms
  if (any(phase_source == SOURCE_thermal_dissipation_ID))     call source_thermal_dissipation_init
  if (any(phase_source == SOURCE_thermal_externalheat_ID))    call source_thermal_externalheat_init

!--------------------------------------------------------------------------------------------------
!initialize kinematic mechanisms
  if (any(phase_kinematics == KINEMATICS_thermal_expansion_ID)) call kinematics_thermal_expansion_init

end subroutine thermal_init


!----------------------------------------------------------------------------------------------
!< @brief calculates thermal dissipation rate
!----------------------------------------------------------------------------------------------
module subroutine constitutive_thermal_getRateAndItsTangents(TDot, dTDot_dT, T, S, Lp, ip, el)
  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(in) :: &
    T
  real(pReal),  intent(in), dimension(:,:,:,:,:) :: &
    S, &                                                                                            !< current 2nd Piola Kirchhoff stress
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), intent(inout) :: &
    TDot, &
    dTDot_dT

  real(pReal) :: &
    my_Tdot, &
    my_dTdot_dT
  integer :: &
    phase, &
    homog, &
    instance, &
    grain, &
    source, &
    constituent

   homog  = material_homogenizationAt(el)
   instance = thermal_typeInstance(homog)

  do grain = 1, homogenization_Ngrains(homog)
     phase = material_phaseAt(grain,el)
     constituent = material_phasememberAt(grain,ip,el)
     do source = 1, phase_Nsources(phase)
       select case(phase_source(source,phase))
         case (SOURCE_thermal_dissipation_ID)
          call source_thermal_dissipation_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                               S(1:3,1:3,grain,ip,el), &
                                                                Lp(1:3,1:3,grain,ip,el), &
                                                                    phase)

         case (SOURCE_thermal_externalheat_ID)
          call source_thermal_externalheat_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                                phase, constituent)

         case default
          my_Tdot = 0.0_pReal
          my_dTdot_dT = 0.0_pReal
       end select
       Tdot = Tdot + my_Tdot
       dTdot_dT = dTdot_dT + my_dTdot_dT
     enddo
   enddo

end subroutine constitutive_thermal_getRateAndItsTangents


end submodule constitutive_thermal
