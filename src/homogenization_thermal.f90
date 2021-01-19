!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) homogenization_thermal


contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_init()


  print'(/,a)',   ' <<<+-  homogenization_thermal init  -+>>>'

  allocate(homogenization_T(discretization_nIPs*discretization_Nelems))
  allocate(homogenization_dot_T(discretization_nIPs*discretization_Nelems))

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
