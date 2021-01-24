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

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters),             dimension(:), allocatable :: &
    param

  public :: &
    thermal_conduction_init, &
    thermal_conduction_getSource, &
    thermal_conduction_putTemperatureAndItsRate, &
    thermal_conduction_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_init()

  integer :: Ninstances,Nmaterialpoints,ho,ip,el,ce
  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogThermal


  print'(/,a)', ' <<<+-  thermal_conduction init  -+>>>'; flush(6)

  Ninstances = count(thermal_type == THERMAL_conduction_ID)
  allocate(param(Ninstances))

  material_homogenization => config_material%get('homogenization')
  do ho = 1, size(material_name_homogenization)
    if (thermal_type(ho) /= THERMAL_conduction_ID) cycle
    homog => material_homogenization%get(ho)
    homogThermal => homog%get('thermal')
    associate(prm => param(thermal_typeInstance(ho)))

#if defined (__GFORTRAN__)
    prm%output = output_asStrings(homogThermal)
#else
    prm%output = homogThermal%get_asStrings('output',defaultVal=emptyStringArray)
#endif

    Nmaterialpoints=count(material_homogenizationAt==ho)

    allocate  (temperature    (ho)%p(Nmaterialpoints), source=thermal_initialT(ho))
    allocate  (temperatureRate(ho)%p(Nmaterialpoints), source=0.0_pReal)

    end associate
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


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_results(homog,group)

  integer,          intent(in) :: homog
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(damage_typeInstance(homog)))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('T')
        call results_writeDataset(group,temperature(homog)%p,'T',&
                                  'temperature','K')
    end select
  enddo outputsLoop
  end associate

end subroutine thermal_conduction_results

end module thermal_conduction
