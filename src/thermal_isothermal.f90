!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isothermal temperature field
!--------------------------------------------------------------------------------------------------
module thermal_isothermal
  use prec
  use config
  use material

  implicit none
  public

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine thermal_isothermal_init

  integer :: h,Nmaterialpoints

  print'(/,a)',   ' <<<+-  thermal_isothermal init  -+>>>'; flush(6)

  do h = 1, size(material_name_homogenization)
    if (thermal_type(h) /= THERMAL_isothermal_ID) cycle

    Nmaterialpoints = count(material_homogenizationAt == h)

    allocate(temperature    (h)%p(Nmaterialpoints),source=thermal_initialT(h))
    allocate(temperatureRate(h)%p(Nmaterialpoints),source = 0.0_pReal)

  enddo

end subroutine thermal_isothermal_init

end module thermal_isothermal
