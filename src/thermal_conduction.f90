!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for temperature evolution from heat conduction
!--------------------------------------------------------------------------------------------------
module thermal_conduction
  use prec
  use material
  use config
  use lattice
  use results
  use constitutive
  use YAML_types
  use discretization

  implicit none
  private

  public :: &
    thermal_conduction_init, &
    thermal_conduction_getSource, &
    thermal_conduction_putTemperatureAndItsRate

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_init()

  integer :: Nmaterialpoints,ho
  class(tNode), pointer :: &
    material_homogenization


  print'(/,a)', ' <<<+-  thermal_conduction init  -+>>>'; flush(6)

  material_homogenization => config_material%get('homogenization')
  do ho = 1, size(material_name_homogenization)
    if (thermal_type(ho) /= THERMAL_conduction_ID) cycle

    Nmaterialpoints=count(material_homogenizationAt==ho)

    allocate  (temperature    (ho)%p(Nmaterialpoints), source=thermal_initialT(ho))
    allocate  (temperatureRate(ho)%p(Nmaterialpoints), source=0.0_pReal)

  enddo

end subroutine thermal_conduction_init


!--------------------------------------------------------------------------------------------------
!> @brief return heat generation rate
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_getSource(Tdot, ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(out) :: &
    Tdot

 integer :: &
    homog

  homog = material_homogenizationAt(el)
  call constitutive_thermal_getRate(TDot, ip,el)

  Tdot = Tdot/real(homogenization_Nconstituents(homog),pReal)

end subroutine thermal_conduction_getSource


!--------------------------------------------------------------------------------------------------
!> @brief updates thermal state with solution from heat conduction PDE
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_putTemperatureAndItsRate(T,Tdot,ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in) :: &
    T, &
    Tdot
  integer :: &
    homog, &
    offset

  homog  = material_homogenizationAt(el)
  offset = material_homogenizationMemberAt(ip,el)
  temperature    (homog)%p(offset) = T
  temperatureRate(homog)%p(offset) = Tdot

end subroutine thermal_conduction_putTemperatureAndItsRate


end module thermal_conduction
