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


end subroutine thermal_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition T onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_partition(T,ip,el)

  real(pReal), intent(in) :: T
  integer,     intent(in) :: &
    ip, &                                                                                           !< integration point
    el                                                                                              !< element number


  call constitutive_thermal_setT(T,1,ip,el)

end subroutine thermal_partition


end submodule homogenization_thermal
