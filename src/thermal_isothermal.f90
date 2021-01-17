!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isothermal temperature field
!--------------------------------------------------------------------------------------------------
module thermal_isothermal
  use prec
  use config
  use material
  use discretization

  implicit none
  public

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine thermal_isothermal_init(T)

  real(pReal), dimension(:), intent(inout) ::  T

  integer :: Ninstances,Nmaterialpoints,ho,ip,el,ce

  print'(/,a)',   ' <<<+-  thermal_isothermal init  -+>>>'; flush(6)

  do ho = 1, size(thermal_type)
    if (thermal_type(ho) /= THERMAL_isothermal_ID) cycle

    Nmaterialpoints = count(material_homogenizationAt == ho)

    allocate(temperature    (ho)%p(Nmaterialpoints),source=thermal_initialT(ho))
    allocate(temperatureRate(ho)%p(Nmaterialpoints),source = 0.0_pReal)

  enddo

  ce = 0
  do el = 1, discretization_Nelems
    do ip = 1, discretization_nIPs
      ce = ce + 1
      ho = material_homogenizationAt(el)
      if (thermal_type(ho) == THERMAL_isothermal_ID) T(ce) = thermal_initialT(ho)
    enddo
  enddo

end subroutine thermal_isothermal_init

end module thermal_isothermal
