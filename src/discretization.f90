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
    discretization_nElem, &
    discretization_nIP

  real(pReal), dimension(:,:), allocatable :: &
    discretization_Centers_disp, &
    discretization_Nodes_disp

  public :: &
    discretization_init, &
    discretization_results

contains
  

subroutine discretization_init(nElem,nIP,nNodes)

  integer, intent(in) :: &
    nElem, &
    nIP, &
    nNodes

  write(6,'(/,a)')   ' <<<+-  discretization init  -+>>>'
   
  discretization_nElem = nElem
  discretization_nIP   = nIP

  allocate(discretization_Centers_disp(3,nIP),source    = 0.0_pReal)
  allocate(discretization_Nodes_disp(  3,nNodes),source = 0.0_pReal)

end subroutine discretization_init


subroutine discretization_results

  call results_writeDataset('current',discretization_Centers_disp,'U','disp','m')
  call results_writeDataset('current',discretization_Nodes_disp,'u','disp','m')

end subroutine discretization_results

end module discretization
