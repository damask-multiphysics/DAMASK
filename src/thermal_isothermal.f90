!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isothermal temperature field
!--------------------------------------------------------------------------------------------------
module thermal_isothermal
  use config
  use material

  implicit none
  public

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine thermal_isothermal_init

  integer :: h,NofMyHomog

  print'(/,a)',   ' <<<+-  thermal_isothermal init  -+>>>'; flush(6)

  do h = 1, material_Nhomogenization
    if (thermal_type(h) /= THERMAL_isothermal_ID) cycle

    NofMyHomog = count(material_homogenizationAt == h)
    thermalState(h)%sizeState = 0
    allocate(thermalState(h)%state0   (0,NofMyHomog))
    allocate(thermalState(h)%subState0(0,NofMyHomog))
    allocate(thermalState(h)%state    (0,NofMyHomog))

    deallocate(temperature    (h)%p)
    allocate  (temperature    (h)%p(1), source=thermal_initialT(h))
    deallocate(temperatureRate(h)%p)
    allocate  (temperatureRate(h)%p(1))

  enddo

end subroutine thermal_isothermal_init

end module thermal_isothermal
