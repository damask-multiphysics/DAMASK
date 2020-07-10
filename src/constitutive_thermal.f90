submodule(constitutive) constitutive_thermal

  implicit none

  interface

  module subroutine source_thermal_dissipation_init
  end subroutine source_thermal_dissipation_init
 
  module subroutine source_thermal_externalheat_init
  end subroutine source_thermal_externalheat_init

  module subroutine kinematics_thermal_expansion_init
  end subroutine kinematics_thermal_expansion_init


  module subroutine source_thermal_dissipation_getRateAndItsTangent(TDot, dTDot_dT, Tstar, Lp, phase)

    integer, intent(in) :: &
      phase
    real(pReal),  intent(in), dimension(3,3) :: &
      Tstar
    real(pReal),  intent(in), dimension(3,3) :: &
      Lp

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

module subroutine thermal_init

! initialize source mechanisms
  if (any(phase_source == SOURCE_thermal_dissipation_ID))     call source_thermal_dissipation_init
  if (any(phase_source == SOURCE_thermal_externalheat_ID))    call source_thermal_externalheat_init
 
!--------------------------------------------------------------------------------------------------
!initialize kinematic mechanisms
  if (any(phase_kinematics == KINEMATICS_thermal_expansion_ID)) call kinematics_thermal_expansion_init

end subroutine thermal_init


module procedure constitutive_thermal_getRateAndItsTangents

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
                                                               Tstar(1:3,1:3,grain,ip,el), &
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
 
end procedure constitutive_thermal_getRateAndItsTangents

end submodule
