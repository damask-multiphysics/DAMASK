!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module discretization_mesh
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdmda.h>
  use PETScdmplex
  use PETScdmda
  use PETScis

  use DAMASK_interface
  use parallelization
  use IO
  use config
  use discretization
  use results
  use FEM_quadrature
  use YAML_types
  use prec

  implicit none
  private

  integer, public, protected :: &
    mesh_Nboundaries, &
    mesh_NcpElemsGlobal

  integer, public, protected :: &
    mesh_NcpElems                                                                                   !< total number of CP elements in mesh

!!!! BEGIN DEPRECATED !!!!!
  integer, public, protected :: &
    mesh_maxNips                                                                                    !< max number of IPs in any CP element
!!!! BEGIN DEPRECATED !!!!!

  real(pReal), dimension(:,:), allocatable :: &
    mesh_ipVolume, &                                                                                !< volume associated with IP (initially!)
    mesh_node0                                                                                      !< node x,y,z coordinates (initially!)

  real(pReal), pointer, dimension(:) :: &
    mesh_node0_temp
  
  real(pReal), dimension(:,:,:), allocatable :: &
    mesh_ipCoordinates                                                                              !< IP x,y,z coordinates (after deformation!)

  DM, public :: geomMesh

  PetscInt, dimension(:), allocatable, public, protected :: &
    mesh_boundaries

  public :: &
    discretization_mesh_init, &
    mesh_FEM_build_ipVolumes, &
    mesh_FEM_build_ipCoordinates

contains


!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine discretization_mesh_init(restart)

  logical, intent(in) :: restart

  integer, allocatable, dimension(:) :: chunkPos
  integer :: dimPlex, &
    mesh_Nnodes, &                                                                                  !< total number of nodes in mesh
    j, l, &
    debug_element, debug_ip
  PetscSF :: sf
  DM :: globalMesh
  PetscInt :: nFaceSets
  PetscInt, pointer, dimension(:) :: pFaceSets
  character(len=pStringLen), dimension(:), allocatable :: fileContent
  IS :: faceSetIS
  PetscErrorCode :: ierr
  integer, dimension(:), allocatable :: &
    materialAt
  class(tNode), pointer :: &
    num_mesh
  integer :: integrationOrder                                                                       !< order of quadrature rule required
  type(tvec)  :: coords_node0

  print'(/,a)',   ' <<<+-  discretization_mesh init  -+>>>'

!--------------------------------------------------------------------------------
! read numerics parameter
  num_mesh => config_numerics%get('mesh',defaultVal=emptyDict)
  integrationOrder = num_mesh%get_asInt('integrationorder',defaultVal = 2)

!---------------------------------------------------------------------------------
! read debug parameters
  debug_element = config_debug%get_asInt('element',defaultVal=1)
  debug_ip      = config_debug%get_asInt('integrationpoint',defaultVal=1)

  call DMPlexCreateFromFile(PETSC_COMM_WORLD,interface_geomFile,PETSC_TRUE,globalMesh,ierr)
  CHKERRQ(ierr)
  call DMGetDimension(globalMesh,dimPlex,ierr)
  CHKERRQ(ierr)
  call DMGetStratumSize(globalMesh,'depth',dimPlex,mesh_NcpElemsGlobal,ierr)
  CHKERRQ(ierr)
  call DMView(globalMesh, PETSC_VIEWER_STDOUT_WORLD,ierr)
  CHKERRQ(ierr)

  ! get number of IDs in face sets (for boundary conditions?)
  call DMGetLabelSize(globalMesh,'Face Sets',mesh_Nboundaries,ierr)
  CHKERRQ(ierr)
  call MPI_Bcast(mesh_Nboundaries,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  call MPI_Bcast(mesh_NcpElemsGlobal,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  call MPI_Bcast(dimPlex,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)

  if (worldrank == 0) then
    call DMClone(globalMesh,geomMesh,ierr)
  else
    call DMPlexDistribute(globalMesh,0,sf,geomMesh,ierr)
  endif
  CHKERRQ(ierr)

  allocate(mesh_boundaries(mesh_Nboundaries), source = 0)
  call DMGetLabelSize(globalMesh,'Face Sets',nFaceSets,ierr)
  CHKERRQ(ierr)
  call DMGetLabelIdIS(globalMesh,'Face Sets',faceSetIS,ierr)
  CHKERRQ(ierr)
  if (nFaceSets > 0) then
    call ISGetIndicesF90(faceSetIS,pFaceSets,ierr)
    CHKERRQ(ierr)
    mesh_boundaries(1:nFaceSets) = pFaceSets
    CHKERRQ(ierr)
    call ISRestoreIndicesF90(faceSetIS,pFaceSets,ierr)
  endif
  call MPI_Bcast(mesh_boundaries,mesh_Nboundaries,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)

  call DMDestroy(globalMesh,ierr); CHKERRQ(ierr)

  call DMGetStratumSize(geomMesh,'depth',dimPlex,mesh_NcpElems,ierr)
  CHKERRQ(ierr)
  call DMGetStratumSize(geomMesh,'depth',0,mesh_Nnodes,ierr)
  CHKERRQ(ierr)

! Get initial nodal coordinates
  call DMGetCoordinates(geomMesh,coords_node0,ierr)
  CHKERRQ(ierr)
  allocate(mesh_node0_temp(dimPlex*mesh_Nnodes))
  call VecGetArrayF90(coords_node0, mesh_node0_temp,ierr)
  CHKERRQ(ierr)
  
  mesh_maxNips = FEM_nQuadrature(dimPlex,integrationOrder)

  call mesh_FEM_build_ipCoordinates(dimPlex,FEM_quadrature_points(dimPlex,integrationOrder)%p)
  call mesh_FEM_build_ipVolumes(dimPlex)

  allocate(materialAt(mesh_NcpElems))
  do j = 1, mesh_NcpElems
    call DMGetLabelValue(geomMesh,'Cell Sets',j-1,materialAt(j),ierr)
    CHKERRQ(ierr)
  end do
  materialAt = materialAt + 1

  if (debug_element < 1 .or. debug_element > mesh_NcpElems) call IO_error(602,ext_msg='element')
  if (debug_ip < 1 .or. debug_ip > mesh_maxNips)            call IO_error(602,ext_msg='IP')

  allocate(mesh_node0(3,mesh_Nnodes),source=0.0_pReal)
  mesh_node0 = reshape(mesh_node0_temp,[dimPlex,mesh_Nnodes])

  call discretization_init(materialAt,&
                           reshape(mesh_ipCoordinates,[3,mesh_maxNips*mesh_NcpElems]), &
                           mesh_node0)

  call writeGeometry(reshape(mesh_ipCoordinates,[3,mesh_maxNips*mesh_NcpElems]),mesh_node0)  

end subroutine discretization_mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume. Allocates global array 'mesh_ipVolume'
!--------------------------------------------------------------------------------------------------
subroutine mesh_FEM_build_ipVolumes(dimPlex)

  PetscInt,intent(in):: dimPlex
  PetscReal          :: vol
  PetscReal, pointer,dimension(:) :: pCent, pNorm
  PetscInt           :: cellStart, cellEnd, cell
  PetscErrorCode     :: ierr

  allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems),source=0.0_pReal)

  call DMPlexGetHeightStratum(geomMesh,0,cellStart,cellEnd,ierr); CHKERRQ(ierr)
  allocate(pCent(dimPlex))
  allocate(pNorm(dimPlex))
  do cell = cellStart, cellEnd-1
    call  DMPlexComputeCellGeometryFVM(geomMesh,cell,vol,pCent,pNorm,ierr)
    CHKERRQ(ierr)
    mesh_ipVolume(:,cell+1) = vol/real(mesh_maxNips,pReal)
  enddo

end subroutine mesh_FEM_build_ipVolumes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP Coordinates. Allocates global array 'mesh_ipCoordinates'
!--------------------------------------------------------------------------------------------------
subroutine mesh_FEM_build_ipCoordinates(dimPlex,qPoints)

  PetscInt,      intent(in) :: dimPlex
  PetscReal,     intent(in) :: qPoints(mesh_maxNips*dimPlex)

  PetscReal,        pointer,dimension(:) :: pV0, pCellJ, pInvcellJ
  PetscReal                 :: detJ
  PetscInt                  :: cellStart, cellEnd, cell, qPt, dirI, dirJ, qOffset
  PetscErrorCode            :: ierr


  allocate(mesh_ipCoordinates(3,mesh_maxNips,mesh_NcpElems),source=0.0_pReal)

  allocate(pV0(dimPlex))
  allocatE(pCellJ(dimPlex**2))
  allocatE(pinvCellJ(dimPlex**2))
  call DMPlexGetHeightStratum(geomMesh,0,cellStart,cellEnd,ierr); CHKERRQ(ierr)
  do cell = cellStart, cellEnd-1                                                                     !< loop over all elements
    call DMPlexComputeCellGeometryAffineFEM(geomMesh,cell,pV0,pCellJ,pInvcellJ,detJ,ierr)
    CHKERRQ(ierr)
    qOffset = 0
    do qPt = 1, mesh_maxNips
      do dirI = 1, dimPlex
        mesh_ipCoordinates(dirI,qPt,cell+1) = pV0(dirI)
        do dirJ = 1, dimPlex
          mesh_ipCoordinates(dirI,qPt,cell+1) = mesh_ipCoordinates(dirI,qPt,cell+1) + &
                                                pCellJ((dirI-1)*dimPlex+dirJ)*(qPoints(qOffset+dirJ) + 1.0)
        enddo
      enddo
      qOffset = qOffset + dimPlex
    enddo
  enddo

end subroutine mesh_FEM_build_ipCoordinates

!--------------------------------------------------------------------------------------------------
!> @brief Write all information needed for the DADF5 geometry
!--------------------------------------------------------------------------------------------------
subroutine writeGeometry(coordinates_points,coordinates_nodes)

  real(pReal), dimension(:,:), intent(in) :: &
  coordinates_nodes, &
  coordinates_points
  
  call results_openJobFile
  call results_closeGroup(results_addGroup('geometry'))
  
  call results_writeDataset(coordinates_nodes,'geometry','x_n', &
        'initial coordinates of the nodes','m')
  
  call results_writeDataset(coordinates_points,'geometry','x_p', &
        'initial coordinates of the materialpoints (cell centers)','m')
  
  call results_closeJobFile
  
  end subroutine writeGeometry

end module discretization_mesh
