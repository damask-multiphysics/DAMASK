!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) homogenization_damage


contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine damage_init()


  print'(/,a)',   ' <<<+-  homogenization_damage init  -+>>>'

  allocate(homogenization_phi(discretization_nIPs*discretization_Nelems))
  allocate(homogenization_dot_phi(discretization_nIPs*discretization_Nelems))

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine damage_partition(phi,ce)

  real(pReal), intent(in) :: phi
  integer,     intent(in) :: ce

  integer :: co

  do co = 1, homogenization_Nconstituents(material_homogenizationAt2(ce))
    call constitutive_damage_set_phi(phi,co,ce)
  enddo

end subroutine damage_partition


end submodule homogenization_damage
