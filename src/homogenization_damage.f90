!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) homogenization_damage

  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: phi
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
module subroutine damage_init()

  class(tNode), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationDamage
  integer :: ho


  print'(/,a)',   ' <<<+-  homogenization_damage init  -+>>>'


  configHomogenizations => config_material%get('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    allocate(current(ho)%phi(count(material_homogenizationAt2==ho)), source=1.0_pReal)
    configHomogenization => configHomogenizations%get(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('damage')) then
        configHomogenizationDamage => configHomogenization%get('damage')
#if defined (__GFORTRAN__)
        prm%output = output_asStrings(configHomogenizationDamage)
#else
        prm%output = configHomogenizationDamage%get_asStrings('output',defaultVal=emptyStringArray)
#endif
      else
        prm%output = emptyStringArray
      endif
    end associate
  enddo

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine damage_partition(ce)

  real(pReal) :: phi
  integer,     intent(in) :: ce

  integer :: co


  phi     = current(material_homogenizationAt2(ce))%phi(material_homogenizationMemberAt2(ce))
  do co = 1, homogenization_Nconstituents(material_homogenizationAt2(ce))
    call constitutive_damage_set_phi(phi,co,ce)
  enddo

end subroutine damage_partition


end submodule homogenization_damage
