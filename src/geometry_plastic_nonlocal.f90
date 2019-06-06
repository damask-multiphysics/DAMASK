!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Geometric information about the IP cells needed for the nonlocal
! plasticity model
!--------------------------------------------------------------------------------------------------
module geometry_plastic_nonlocal
  use prec

  implicit none
  private

  real(pReal), dimension(:,:),     allocatable, public, protected :: &
    geometry_plastic_nonlocal_IPvolume0                                                             !< volume associated with IP (initially!)
 
  real(pReal), dimension(:,:,:),   allocatable, public, protected :: &
    geometry_plastic_nonlocal_IParea0                                                               !< area of interface to neighboring IP (initially!)
  
  real(pReal), dimension(:,:,:,:), allocatable, public, protected :: &
    geometry_plastic_nonlocal_IPareaNormal0                                                         !< area normal of interface to neighboring IP (initially!)
  
  integer,     dimension(:,:,:,:), allocatable, public, protected :: &
    geometry_plastic_nonlocal_IPneighborhood                                                        !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]


  public :: &
    geometry_plastic_nonlocal_set_IPneighborhood, &
    geometry_plastic_nonlocal_set_IPvolume
    
  contains
  
subroutine geometry_plastic_nonlocal_set_IPneighborhood(IPneighborhood)

  integer, dimension(:,:,:,:), intent(in) :: IPneighborhood

  geometry_plastic_nonlocal_IPneighborhood = IPneighborhood

end subroutine geometry_plastic_nonlocal_set_IPneighborhood


subroutine geometry_plastic_nonlocal_set_IPvolume(IPvolume)

  real(pReal), dimension(:,:), intent(in) :: IPvolume

  geometry_plastic_nonlocal_IPvolume0 = IPvolume

end subroutine geometry_plastic_nonlocal_set_IPvolume


end module geometry_plastic_nonlocal
