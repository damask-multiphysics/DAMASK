!--------------------------------------------------------------------------------------------------
!> @brief spatial discretization
!--------------------------------------------------------------------------------------------------
module discretization

  use, intrinsic :: iso_c_binding
  use prec
  use results

  implicit none
  private
  
  integer, public, protected :: &
    discretization_nIP, &
    discretization_nElem


  real(pReal), dimension(:,:), allocatable :: &
    discretization_IPcoords0, &
    discretization_NodeCoords0, &
    discretization_IPcoords, &
    discretization_NodeCoords

  public :: &
    discretization_init, &
    discretization_results, &
    discretization_setIPcoords

contains
  

subroutine discretization_init(nElem,IPcoords0,NodeCoords0)

  integer,                     intent(in) :: &
    nElem
  real(pReal), dimension(:,:), intent(in) :: &
    IPcoords0, &
    NodeCoords0

  write(6,'(/,a)')   ' <<<+-  discretization init  -+>>>'
   
  discretization_nElem = nElem 
  discretization_nIP   = size(IPcoords0,2)

  discretization_IPcoords0   = IPcoords0
  discretization_IPcoords    = IPcoords0
  discretization_NodeCoords0 = NodeCoords0
  discretization_NodeCoords  = NodeCoords0
  
end subroutine discretization_init


subroutine discretization_results

  real(pReal), dimension(:,:), allocatable :: u
  
  u =  discretization_NodeCoords -discretization_NodeCoords0
  call results_writeDataset('current',U,'U','nodal displacements','m')
  
  u = discretization_IPcoords -discretization_IPcoords0
  call results_writeDataset('current',u,'u','IP displacements','m')

end subroutine discretization_results


subroutine discretization_setIPcoords(IPcoords)

  real(pReal), dimension(:,:), intent(in) :: IPcoords
  
  discretization_IPcoords = IPcoords

end subroutine discretization_setIPcoords


end module discretization
