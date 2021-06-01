!--------------------------------------------------------------------------------------------------
!> @brief spatial discretization
!> @details serves as an abstraction layer between the different solvers and DAMASK
!--------------------------------------------------------------------------------------------------
module discretization

  use prec
  use results

  implicit none
  private

  integer,     public, protected :: &
    discretization_nIPs, &
    discretization_Nelems

  integer,     public, protected, dimension(:),   allocatable :: &
    discretization_materialAt                                                                       !ToDo: discretization_ID_material

  real(pReal), public, protected, dimension(:,:), allocatable :: &
    discretization_IPcoords0, &
    discretization_IPcoords, &
    discretization_NodeCoords0, &
    discretization_NodeCoords

  integer :: &
    discretization_sharedNodesBegin

  public :: &
    discretization_init, &
    discretization_results, &
    discretization_setIPcoords, &
    discretization_setNodeCoords

contains

!--------------------------------------------------------------------------------------------------
!> @brief stores the relevant information in globally accesible variables
!--------------------------------------------------------------------------------------------------
subroutine discretization_init(materialAt,&
                               IPcoords0,NodeCoords0,&
                               sharedNodesBegin)

  integer,     dimension(:),   intent(in) :: &
    materialAt
  real(pReal), dimension(:,:), intent(in) :: &
    IPcoords0, &
    NodeCoords0
  integer, optional,           intent(in) :: &
    sharedNodesBegin                                                                                !< index of first node shared among different processes (MPI)

  print'(/,a)', ' <<<+-  discretization init  -+>>>'; flush(6)

  discretization_Nelems = size(materialAt,1)
  discretization_nIPs   = size(IPcoords0,2)/discretization_Nelems

  discretization_materialAt = materialAt

  discretization_IPcoords0   = IPcoords0
  discretization_IPcoords    = IPcoords0

  discretization_NodeCoords0 = NodeCoords0
  discretization_NodeCoords  = NodeCoords0

  if(present(sharedNodesBegin)) then
    discretization_sharedNodesBegin = sharedNodesBegin
  else
    discretization_sharedNodesBegin = size(discretization_NodeCoords0,2)
  endif

end subroutine discretization_init


!--------------------------------------------------------------------------------------------------
!> @brief write the displacements
!--------------------------------------------------------------------------------------------------
subroutine discretization_results

  real(pReal), dimension(:,:), allocatable :: u

  call results_closeGroup(results_addGroup('current/geometry'))

  u = discretization_NodeCoords (1:3,:discretization_sharedNodesBegin) &
    - discretization_NodeCoords0(1:3,:discretization_sharedNodesBegin)
  call results_writeDataset(u,'current/geometry','u_n','displacements of the nodes','m')

  u = discretization_IPcoords &
    - discretization_IPcoords0
  call results_writeDataset(u,'current/geometry','u_p','displacements of the materialpoints (cell centers)','m')

end subroutine discretization_results


!--------------------------------------------------------------------------------------------------
!> @brief stores current IP coordinates
!--------------------------------------------------------------------------------------------------
subroutine discretization_setIPcoords(IPcoords)

  real(pReal), dimension(:,:), intent(in) :: IPcoords

  discretization_IPcoords = IPcoords

end subroutine discretization_setIPcoords


!--------------------------------------------------------------------------------------------------
!> @brief stores current IP coordinates
!--------------------------------------------------------------------------------------------------
subroutine discretization_setNodeCoords(NodeCoords)

  real(pReal), dimension(:,:), intent(in) :: NodeCoords

  discretization_NodeCoords = NodeCoords

end subroutine discretization_setNodeCoords


end module discretization
