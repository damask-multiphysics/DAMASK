submodule(constitutive) constitutive_thermal

  implicit none

  interface

  module subroutine source_thermal_dissipation_init
  end subroutine source_thermal_dissipation_init
 
  module subroutine source_thermal_externalheat_init
  end subroutine source_thermal_externalheat_init

  module subroutine kinematics_thermal_expansion_init
  end subroutine kinematics_thermal_expansion_init

  module subroutine source_thermal_externalheat_dotState(phase, of)
    integer, intent(in) :: &
      phase, &
      of
  end subroutine source_thermal_externalheat_dotState


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

!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
module function thermal_dotState(S, FArray, Fi, FpArray, subdt, ipc, ip, el,phase,of) result(broken_thermal)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                              !< element
    phase, &
    of
  real(pReal),  intent(in) :: &
    subdt                                                                                           !< timestep
  real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,discretization_nIP,discretization_nElem) :: &
    FArray, &                                                                                       !< elastic deformation gradient
    FpArray                                                                                         !< plastic deformation gradient
  real(pReal),  intent(in), dimension(3,3) :: &
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),  intent(in), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
  logical :: broken_thermal
  integer :: i

  SourceLoop: do i = 1, phase_Nsources(phase)

    sourceType: select case (phase_source(i,phase))

      case (SOURCE_thermal_externalheat_ID) sourceType
        call source_thermal_externalheat_dotState(phase,of)

    end select sourceType

    broken_thermal = any(IEEE_is_NaN(sourceState(phase)%p(i)%dotState(:,of)))

  enddo sourceLoop

end function thermal_dotState


module subroutine thermal_source_getRateAndItsTangents(Tdot, dTdot_dT, T, Tstar, Lp, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(in) :: &
    T
  real(pReal),  intent(in), dimension(:,:,:,:,:) :: &
    Tstar, &
    Lp
  real(pReal), intent(inout) :: &
    Tdot, &
    dTdot_dT

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
 
end subroutine thermal_source_getRateAndItsTangents

end submodule
