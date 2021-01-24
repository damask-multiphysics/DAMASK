!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) homogenization_thermal


  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: T, dot_T
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: current

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters),             dimension(:), allocatable :: &
    param


contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_init()

  class(tNode), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationThermal
  integer :: ho


  print'(/,a)',   ' <<<+-  homogenization_thermal init  -+>>>'


  allocate(homogenization_T(discretization_nIPs*discretization_Nelems))
  allocate(homogenization_dot_T(discretization_nIPs*discretization_Nelems))

  configHomogenizations => config_material%get('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    allocate(current(ho)%T(count(material_homogenizationAt2==ho)), source=thermal_initialT(ho))
    allocate(current(ho)%dot_T(count(material_homogenizationAt2==ho)), source=0.0_pReal)
    configHomogenization => configHomogenizations%get(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('thermal')) then
        configHomogenizationThermal => configHomogenization%get('thermal')
#if defined (__GFORTRAN__)
        prm%output = output_asStrings(configHomogenizationThermal)
#else
        prm%output = configHomogenizationThermal%get_asStrings('output',defaultVal=emptyStringArray)
#endif
      else
        prm%output = emptyStringArray
      endif
    end associate
  enddo

end subroutine thermal_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_partition(T,ce)

  real(pReal), intent(in) :: T
  integer,     intent(in) :: ce

  integer :: co

  do co = 1, homogenization_Nconstituents(material_homogenizationAt2(ce))
    call constitutive_thermal_setT(T,co,ce)
  enddo

end subroutine thermal_partition


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize temperature rates
!--------------------------------------------------------------------------------------------------
module subroutine thermal_homogenize(ip,el)

  integer, intent(in) :: ip,el

  call constitutive_thermal_getRate(homogenization_dot_T((el-1)*discretization_nIPs+ip), ip,el)

end subroutine thermal_homogenize


end submodule homogenization_thermal
