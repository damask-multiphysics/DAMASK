!--------------------------------------------------------------------------------------------------
!> @brief spatial discretization
!> @details serves as an abstraction layer between the different solvers and DAMASK
!--------------------------------------------------------------------------------------------------
module discretization

  use prec
  use results
#if defined(PETSc) || defined(DAMASK_HDF5)
  use HDF5_utilities
#endif

  implicit none
  private
  
  integer,     public, protected :: &
    discretization_nIP, &
    discretization_nElem
    
  integer,     public, protected, dimension(:),   allocatable :: &
    discretization_homogenizationAt, &
    discretization_microstructureAt   

  real(pReal), public, protected, dimension(:,:), allocatable :: & 
    discretization_IPcoords0, &
    discretization_NodeCoords0, &
    discretization_IPcoords, &
    discretization_NodeCoords

  public :: &
    discretization_init, &
    discretization_results, &
    discretization_setIPcoords

contains
  
!--------------------------------------------------------------------------------------------------
!> @brief stores the relevant information in globally accesible variables
!--------------------------------------------------------------------------------------------------
subroutine discretization_init(homogenizationAt,microstructureAt,IPcoords0,NodeCoords0)

  integer,     dimension(:),   intent(in) :: &
    homogenizationAt, &
    microstructureAt
  real(pReal), dimension(:,:), intent(in) :: &
    IPcoords0, &
    NodeCoords0

  write(6,'(/,a)')   ' <<<+-  discretization init  -+>>>'

  discretization_nElem = size(microstructureAt,1)
  discretization_nIP   = size(IPcoords0,2)/discretization_nElem

  discretization_homogenizationAt = homogenizationAt
  discretization_microstructureAt = microstructureAt  

  discretization_IPcoords0   = IPcoords0
  discretization_IPcoords    = IPcoords0

  discretization_NodeCoords0 = NodeCoords0
  discretization_NodeCoords  = NodeCoords0
  
end subroutine discretization_init


!--------------------------------------------------------------------------------------------------
!> @brief write the displacements
!--------------------------------------------------------------------------------------------------
subroutine discretization_results
#if defined(PETSc) || defined(DAMASK_HDF5)
  real(pReal), dimension(:,:), allocatable :: u
  
  call HDF5_closeGroup(results_addGroup(trim('current/geometry')))
  
  u =  discretization_NodeCoords - discretization_NodeCoords0
  call results_writeDataset('current/geometry',u,'u_n','nodal displacements','m')
  
  u = discretization_IPcoords - discretization_IPcoords0
  call results_writeDataset('current/geometry',u,'u_c','cell center displacements','m')
#endif
end subroutine discretization_results


!--------------------------------------------------------------------------------------------------
!> @brief stores current IP coordinates
!--------------------------------------------------------------------------------------------------
subroutine discretization_setIPcoords(IPcoords)

  real(pReal), dimension(:,:), intent(in) :: IPcoords
  
  discretization_IPcoords = IPcoords

end subroutine discretization_setIPcoords

end module discretization