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
 use IO
 use debug
 use discretization
 use numerics
 use FEsolving
 use FEM_quadrature
 use prec
 
 implicit none
 private
 
 integer, public, protected :: &
   mesh_Nboundaries, &
   mesh_NcpElemsGlobal

 integer :: &
   mesh_NcpElems                                                                                     !< total number of CP elements in mesh

!!!! BEGIN DEPRECATED !!!!!
 integer, public, protected :: &
   mesh_maxNips                                                                                     !< max number of IPs in any CP element
!!!! BEGIN DEPRECATED !!!!!

 integer, dimension(:,:), allocatable :: &
   mesh_element !DEPRECATED

 real(pReal), dimension(:,:), allocatable :: &
   mesh_ipVolume, &                                                                                 !< volume associated with IP (initially!)
   mesh_node0                                                                                       !< node x,y,z coordinates (initially!)
 
 real(pReal), dimension(:,:,:), allocatable :: &
   mesh_ipCoordinates                                                                               !< IP x,y,z coordinates (after deformation!)

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

  integer, intent(in) :: restart

  integer, dimension(1), parameter:: FE_geomtype = [1]                                              !< geometry type of particular element type
  integer, dimension(1) :: FE_Nips                                                                  !< number of IPs in a specific type of element
  
  integer, allocatable, dimension(:) :: chunkPos
  integer :: dimPlex, &
    mesh_Nnodes, &                                                                                  !< total number of nodes in mesh
    j, l
  integer, parameter :: &
    mesh_ElemType=1                                                                                 !< Element type of the mesh (only support homogeneous meshes)
  PetscSF :: sf
  DM :: globalMesh
  PetscInt :: nFaceSets
  PetscInt, pointer, dimension(:) :: pFaceSets
  character(len=pStringLen), dimension(:), allocatable :: fileContent
  IS :: faceSetIS 
  PetscErrorCode :: ierr

  
  write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'

  ! read in file
  call DMPlexCreateFromFile(PETSC_COMM_WORLD,geometryFile,PETSC_TRUE,globalMesh,ierr)
  CHKERRQ(ierr)
  ! get spatial dimension (2 or 3?)
  call DMGetDimension(globalMesh,dimPlex,ierr)
  CHKERRQ(ierr)
  write(6,*) 'dimension',dimPlex;flush(6)
  call DMGetStratumSize(globalMesh,'depth',dimPlex,mesh_NcpElemsGlobal,ierr)
  CHKERRQ(ierr)
  ! get number of IDs in face sets (for boundary conditions?)
  call DMGetLabelSize(globalMesh,'Face Sets',mesh_Nboundaries,ierr)
  CHKERRQ(ierr)
  write(6,*) 'number of "Face Sets"',mesh_Nboundaries;flush(6)
  call MPI_Bcast(mesh_Nboundaries,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  call MPI_Bcast(mesh_NcpElemsGlobal,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
  call MPI_Bcast(dimPlex,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)

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

  if (worldrank == 0) then
    fileContent = IO_read_ASCII(geometryFile)
    l = 0
    do
      l = l + 1
      if (IO_isBlank(fileContent(l))) cycle         ! need also to ignore C and C++ style comments?
      if (trim(fileContent(l)) == '$Elements') then
        j = 0
        l = l + 1
        do
          l = l + 1
          if (trim(fileContent(l)) == '$EndElements') exit
          chunkPos = IO_stringPos(fileContent(l))
          if (chunkPos(1) == 3+IO_intValue(fileContent(l),chunkPos,3)+dimPlex+1) then
            call DMSetLabelValue(globalMesh,'material',j,IO_intValue(fileContent(l),chunkPos,4),ierr)
            CHKERRQ(ierr)
            j = j + 1
          endif
        enddo
        exit
      endif
    enddo
    call DMClone(globalMesh,geomMesh,ierr)
    CHKERRQ(ierr)
  else 
    call DMPlexDistribute(globalMesh,0,sf,geomMesh,ierr)
    CHKERRQ(ierr)
  endif  

  call DMDestroy(globalMesh,ierr); CHKERRQ(ierr)
  
  call DMGetStratumSize(geomMesh,'depth',dimPlex,mesh_NcpElems,ierr)
  CHKERRQ(ierr)
  call DMGetStratumSize(geomMesh,'depth',0,mesh_Nnodes,ierr)
  CHKERRQ(ierr)

  FE_Nips(FE_geomtype(1)) = FEM_nQuadrature(dimPlex,integrationOrder)
  mesh_maxNips = FE_Nips(1)
  
  write(6,*) 'mesh_maxNips',mesh_maxNips
  call mesh_FEM_build_ipCoordinates(dimPlex,FEM_quadrature_points(dimPlex,integrationOrder)%p)
  call mesh_FEM_build_ipVolumes(dimPlex)
  
  allocate (mesh_element (4,mesh_NcpElems)); mesh_element = 0
  do j = 1, mesh_NcpElems
    mesh_element( 1,j) = -1                                                                     ! DEPRECATED
    mesh_element( 2,j) = mesh_elemType                                                               ! elem type
    mesh_element( 3,j) = 1                                                                      ! homogenization
    call DMGetLabelValue(geomMesh,'material',j-1,mesh_element(4,j),ierr)
    CHKERRQ(ierr)
  end do 

  if (debug_e < 1 .or. debug_e > mesh_NcpElems) &
    call IO_error(602,ext_msg='element')                                                        ! selected element does not exist
  if (debug_i < 1 .or. debug_i > FE_Nips(FE_geomtype(mesh_element(2,debug_e)))) &
    call IO_error(602,ext_msg='IP')                                                             ! selected element does not have requested IP
  
  FEsolving_execElem = [ 1,mesh_NcpElems ]                                                      ! parallel loop bounds set to comprise all DAMASK elements
  FEsolving_execIP   = [1,FE_Nips(FE_geomtype(mesh_element(2,1)))]
  
  allocate(mesh_node0(3,mesh_Nnodes),source=0.0_pReal)

  call discretization_init(mesh_element(3,:),mesh_element(4,:),&
                           reshape(mesh_ipCoordinates,[3,mesh_maxNips*mesh_NcpElems]), &
                           mesh_node0)
 
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

end module discretization_mesh
