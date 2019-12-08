!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isothermal temperature field
!--------------------------------------------------------------------------------------------------
module thermal_isothermal
  use prec
  use config
  use material

  implicit none
  private
  
  public :: &
    thermal_isothermal_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine thermal_isothermal_init
 
  integer :: &
    homog, &
    NofMyHomog

  write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_isothermal_label//' init  -+>>>'

  initializeInstances: do homog = 1, material_Nhomogenization
    
    if (thermal_type(homog) /= THERMAL_isothermal_ID) cycle
    NofMyHomog = count(material_homogenizationAt == homog)
    thermalState(homog)%sizeState = 0
    allocate(thermalState(homog)%state0   (0,NofMyHomog), source=0.0_pReal)
    allocate(thermalState(homog)%subState0(0,NofMyHomog), source=0.0_pReal)
    allocate(thermalState(homog)%state    (0,NofMyHomog), source=0.0_pReal)
      
    deallocate(temperature    (homog)%p)
    allocate  (temperature    (homog)%p(1), source=thermal_initialT(homog))
    deallocate(temperatureRate(homog)%p)
    allocate  (temperatureRate(homog)%p(1), source=0.0_pReal)

  enddo initializeInstances

end subroutine thermal_isothermal_init

end module thermal_isothermal
