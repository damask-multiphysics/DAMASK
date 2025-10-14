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
  use PETScDMplex
  use PETScDMDA
  use PETScIS
#if (PETSC_VERSION_MAJOR==3 && (PETSC_VERSION_MINOR>=18 && PETSC_VERSION_MINOR<23))
  use PETScDT
#endif
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  use MPI_f08
#endif

  use CLI
  use prec
  use parallelization
  use config
  use IO
  use discretization
  use result
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<18)
  use FEM_quadrature
#endif
  use types
  use prec

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  implicit none(type,external)
#else
  implicit none
#endif
  private

  PetscInt, public, protected :: &
    mesh_Nboundaries                                                                                !< Number of defined BC (total)

  PetscInt, dimension(2), public, protected :: &
    mesh_BCTypeSetSize                                                                              !< Number of each BC type (1: Vertex, 2: Faces)

  PetscInt, dimension(:), allocatable, public, protected :: &
    mesh_boundariesIS, &                                                                            !< Index Set (tag values) of BC in mesh file
    mesh_boundariesIdx                                                                              !< BC Type index (BCType_Vertex, BCType_Face)

  PetscInt, public, protected :: &
    mesh_NcpElemsGlobal

  PetscInt, public, protected :: &
    mesh_NcpElems                                                                                   !< total number of CP elements in mesh

!!!! BEGIN DEPRECATED !!!!!
  PetscInt, public, protected :: &
    mesh_maxNips                                                                                    !< max number of IPs in any CP element
!!!! END DEPRECATED !!!!!

  PetscInt, parameter, public :: &
    mesh_BCTypeFace   = 1_pPETSCINT, &
    mesh_BCTypeVertex = 2_pPETSCINT

  character(len=*), dimension(2), public, parameter :: &
    mesh_BCTypeLabel = ['Face Sets  ', 'Vertex Sets']                                               !< PETSc default BC labels

  DM, public :: geomMesh

#if (PETSC_VERSION_MINOR<16)
  external :: &
    DMDestroy
#endif
#if (PETSC_VERSION_MINOR<23)
  external :: &
#if (PETSC_VERSION_MINOR<22)
    DMAddField, &
#endif
    PetscFEGetDimension, &
    PetscFEDestroy, &
    PetscFEGetDualSpace, &
    PetscDualSpaceGetFunctional
#endif

  public :: &
    discretization_mesh_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize the mesh.
!--------------------------------------------------------------------------------------------------
subroutine discretization_mesh_init()

  PetscInt :: dimPlex, &
    mesh_Nnodes, &                                                                                  !< total number of nodes in mesh
    j
  PetscSF :: sf
  DM :: globalMesh
  PetscInt :: NelemsGlobal, Nelems
  PetscInt :: label
  IS :: setIS
  PetscInt, pointer, dimension(:) :: pSets
#if (PETSC_VERSION_MINOR>=18)
  PetscQuadrature :: quadrature
#endif
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscInt, dimension(:), allocatable :: &
    materialAt
  type(tDict), pointer :: &
    num_solver, &
    num_mesh
  PetscInt :: p_i, p_s
  integer:: dim
#if (PETSC_VERSION_MINOR>=18)
  real(pREAL), dimension(:),     pointer     :: qPointsP
#endif
#if (PETSC_VERSION_MINOR<23)
  real(pREAL), dimension(:),     pointer     :: PETSC_NULL_REAL_POINTER => NULL()
#endif
  real(pREAL), dimension(:,:),   allocatable :: v_0, &                                              ! volume associated with IP (initially!)
                                                x_n                                                 ! node x,y,z coordinates (initially!)
  real(pREAL), dimension(:,:,:), allocatable :: x_p                                                 ! IP x,y,z coordinates (after deformation!)
  PetscInt,    dimension(:,:),   allocatable :: T_e                                                 ! element connectivity (node numbers in each cell)


  print'(/,1x,a)',   '<<<+-  discretization_mesh init  -+>>>'

!--------------------------------------------------------------------------------
! read numerics parameter
  num_solver => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_mesh   => num_solver%get_dict('mesh',defaultVal=emptyDict)
  p_i = int(num_mesh%get_asInt('p_i',defaultVal=2),pPETSCINT)
  p_s = int(num_mesh%get_asInt('p_s',defaultVal=2),pPETSCINT)

  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,                        &
                                ' -dm_plex_filename ' // CLI_geomFile // ' &
                                & -dm_plex_interpolate 1                   &
                                & -dm_plex_gmsh_mark_vertices ',           &
                                err_PETSc)
  call DMCreate(PETSC_COMM_WORLD, globalMesh, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetType(globalMesh, DMPLEX, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetFromOptions(globalMesh, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDimension(globalMesh,dimPlex,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetStratumSize(globalMesh,'depth',dimPlex,NelemsGlobal,err_PETSc)
  CHKERRQ(err_PETSc)
  mesh_NcpElemsGlobal = NelemsGlobal

  dim = int(dimPlex)
  call MPI_Bcast(dim,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  dimPlex = int(dim,pPETSCINT)

  if (worldsize == 1) then
    call DMClone(globalMesh,geomMesh,err_PETSc)
  else
    call DMPlexDistribute(globalMesh,0_pPETSCINT,sf,geomMesh,err_PETSc)
  end if
  CHKERRQ(err_PETSc)

  mesh_Nboundaries = 0_pPETSCINT
  allocate(mesh_boundariesIS(mesh_Nboundaries))
  mesh_BCTypeSetSize = 0_pPETSCINT
  do label = 1_pPETSCINT, size(mesh_BCTypeLabel)
    call DMGetLabelSize(globalMesh,mesh_BCTypeLabel(label),mesh_BCTypeSetSize(label),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGetLabelIdIS(globalMesh,mesh_BCTypeLabel(label),setIS,err_PETSc)
    CHKERRQ(err_PETSc)
    if (mesh_BCTypeSetSize(label) > 0_pPETSCINT) then
      call ISGetIndices(setIS,pSets,err_PETSc)
      CHKERRQ(err_PETSc)
      mesh_boundariesIS = [mesh_boundariesIS, pSets]
      call ISRestoreIndices(setIS,pSets,err_PETSc)
      CHKERRQ(err_PETSc)
      mesh_Nboundaries = mesh_Nboundaries+mesh_BCTypeSetSize(label)
    end if
  end do
  allocate(mesh_boundariesIdx(mesh_Nboundaries), source = mesh_BCTypeFace)
  mesh_boundariesIdx(mesh_BCTypeSetSize(mesh_BCTypeFace)+1_pPETSCINT:mesh_Nboundaries) = mesh_BCTypeVertex

  call MPI_Bcast(mesh_BCTypeSetSize,int(size(mesh_BCTypeSetSize),kind=MPI_INTEGER_KIND),MPI_INTEGER,0_MPI_INTEGER_KIND, &
                 MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(mesh_Nboundaries,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(mesh_boundariesIdx,int(mesh_Nboundaries,kind=MPI_INTEGER_KIND),MPI_INTEGER,0_MPI_INTEGER_KIND, &
                 MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(mesh_boundariesIS,int(mesh_Nboundaries),MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  call DMDestroy(globalMesh,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMGetStratumSize(geomMesh,'depth',dimPlex,Nelems,err_PETSc)
  CHKERRQ(err_PETSc)
  mesh_NcpElems = Nelems
  call DMGetStratumSize(geomMesh,'depth',0_pPETSCINT,mesh_Nnodes,err_PETSc)
  CHKERRQ(err_PETSc)

! Get initial nodal coordinates
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=18)
  call PetscDTSimplexQuadrature(dimplex, p_i, PETSCDTSIMPLEXQUAD_DEFAULT, quadrature, err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MINOR>=22)
  call PetscQuadratureGetData(quadrature,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                              mesh_maxNips,qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#else
  call PetscQuadratureGetData(quadrature,PETSC_NULL_INTEGER(1),PETSC_NULL_INTEGER(1), &
                              mesh_maxNips,qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  x_p = build_coordinates_IP(dimPlex,qPointsP)
  v_0 = build_volume_IP(dimPlex)
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  call build_nodes_and_connectivity(x_n,T_e, &
                                    p_s, geomMesh)
#else
  call build_nodes_and_connectivity(x_n,p_s)
#endif

#if (PETSC_VERSION_MINOR>=22)
  call PetscQuadratureRestoreData(quadrature,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                                  PETSC_NULL_INTEGER,qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#else
  call PetscQuadratureRestoreData(quadrature,PETSC_NULL_INTEGER(1),PETSC_NULL_INTEGER(1), &
                                  PETSC_NULL_INTEGER(1),qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  call PetscQuadratureDestroy(quadrature, err_PETSc)
  CHKERRQ(err_PETSc)
#else
  mesh_maxNips = FEM_nQuadrature(dimPlex,p_i)

  x_p = build_coordinates_IP(dimPlex,FEM_quadrature_points(dimPlex,p_i)%p)
  v_0 = build_volume_IP(dimPlex)
#endif

  allocate(materialAt(mesh_NcpElems))
  do j = 1, mesh_NcpElems
    call DMGetLabelValue(geomMesh,'Cell Sets',j-1,materialAt(j),err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  call discretization_init(int(materialAt),&
                           reshape(x_p,[3,int(mesh_maxNips*mesh_NcpElems)]), &
                           x_n)

  call writeGeometry(reshape(x_p,[3,int(mesh_maxNips*mesh_NcpElems)]), &
                     x_n,T_e)

end subroutine discretization_mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP volume.
!--------------------------------------------------------------------------------------------------
function build_volume_IP(dimPlex) result(v_0)

  real(pREAL), dimension(:,:), allocatable :: v_0
  PetscInt,    intent(in) :: dimPlex

  PetscReal      :: vol
  PetscInt       :: cellStart, cellEnd, cell
  PetscErrorCode :: err_PETSc
  PetscReal, pointer,dimension(:) :: pCent, pNorm


  allocate(v_0(mesh_maxNips,mesh_NcpElems),source=0.0_pREAL)

  call DMPlexGetHeightStratum(geomMesh,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
  allocate(pCent(dimPlex))
  allocate(pNorm(dimPlex))
  do cell = cellStart, cellEnd-1
    call DMPlexComputeCellGeometryFVM(geomMesh,cell,vol,pCent,pNorm,err_PETSc)
    CHKERRQ(err_PETSc)
    v_0(:,cell+1) = vol/real(mesh_maxNips,pREAL)
  end do

end function build_volume_IP


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP Coordinates.
!--------------------------------------------------------------------------------------------------
function build_coordinates_IP(dimPlex,qPoints) result(x_p)

  PetscInt,                                     intent(in) :: dimPlex
  PetscReal,   dimension(mesh_maxNips*dimPlex), intent(in) :: qPoints
  real(pREAL), dimension(:,:,:), allocatable :: x_p


  PetscReal, pointer, dimension(:) :: pV0, pCellJ, pInvcellJ
  PetscReal      :: detJ
  PetscInt       :: cellStart, cellEnd, cell, qPt, dirI, dirJ, qOffset
  PetscErrorCode :: err_PETSc


  allocate(x_p(3,mesh_maxNips,mesh_NcpElems),source=0.0_pREAL)

  allocate(pV0(dimPlex))
  allocatE(pCellJ(dimPlex**2))
  allocatE(pinvCellJ(dimPlex**2))
  call DMPlexGetHeightStratum(geomMesh,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
  do cell = cellStart, cellEnd-1_pPETSCINT                                                          ! loop over all elements
    call DMPlexComputeCellGeometryAffineFEM(geomMesh,cell,pV0,pCellJ,pInvcellJ,detJ,err_PETSc)
    CHKERRQ(err_PETSc)
    qOffset = 0
    do qPt = 1_pPETSCINT, mesh_maxNips
      do dirI = 1_pPETSCINT, dimPlex
        x_p(dirI,qPt,cell+1) = pV0(dirI)
        do dirJ = 1_pPETSCINT, dimPlex
          x_p(dirI,qPt,cell+1) = x_p(dirI,qPt,cell+1) &
                               + pCellJ((dirI-1)*dimPlex+dirJ)*(qPoints(qOffset+dirJ) + 1.0_pREAL)
        end do
      end do
      qOffset = qOffset + dimPlex
    end do
  end do

end function build_coordinates_IP


!--------------------------------------------------------------------------------------------------
!> @brief Get mesh node coordinates and element connectivity matrix.
!--------------------------------------------------------------------------------------------------
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
subroutine build_nodes_and_connectivity(x_n, connectivity_matrix, &
                                       p_s, geomMesh)
#else
subroutine build_nodes_and_connectivity(x_n, p_s)
#endif

  real(pREAL), dimension(:,:), allocatable, intent(out) :: &
    x_n                                                                                             !< mesh nodes coordinates (including high-order approximation)
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  PetscInt,    dimension(:,:), allocatable, intent(out) :: &
    connectivity_matrix                                                                             !< element connectivity (node numbers in each cell)
#endif
  PetscInt, intent(in) :: &
    p_s                                                                                             !< order of approximation space
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  DM,       intent(in) :: &
    geomMesh
#endif

  DM       :: coordDM                                                                               ! approximation space DM
  Vec      :: coordVec                                                                              ! local nodes coordinates
  PetscInt :: coordDim,    &                                                                        ! DOF-per-node
              nLocalNodes, &                                                                        ! local number of nodes
              nCellNodes,  &                                                                        ! number of nodes in a cell
              feDim,       &                                                                        ! DOF per cell (nNodes x DOF per node)
              feBasis,     &
              cellStart, cellEnd, cell
  PetscDS  :: coordDS
  PetscFE  :: coordFE
  PetscSection    :: globalSection, localSection                                                    ! section (to retrieve DOF)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>23)
  DMPolytopeType  :: cell_type                                                                      ! tri, quad, tet, hex
#endif
  PetscDualSpace  :: coordDualSpace
  PetscQuadrature :: refQuadrature
  PetscErrorCode  :: err_PETSc, ierr

  real(pREAL), dimension(:), pointer     :: coords, &                                               ! local nodes coordinates
                                            nodeCoords                                              ! single node coordinates
  real(pREAL), dimension(:), allocatable :: refCoords                                               ! node coordinates in reference element [-1,+1]^d
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=23)
  real(pREAL), dimension(:), allocatable :: mappedCoords                                            ! real (mesh) node coordinates
#else
  real(pREAL), dimension(:), allocatable, target :: mappedCoords                                    ! real (mesh) node coordinates
  real(pREAL), dimension(:), pointer     :: pMappedCoords, &
                                            PETSC_NULL_REAL_POINTER => NULL()
#endif
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>23)
  PetscInt,    dimension(:), pointer     :: indices                                                 ! cell closure DOF indices
  integer,     dimension(:), allocatable :: node_map                                                ! PETSc to VTK node order mapping
#endif


  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                '-coord_petscspace_degree ' // IO_intAsStr(int(p_s)) // ' &
                                &-coord_petscdualspace_lagrange_node_type equispaced &
                                &-coord_petscdualspace_lagrange_node_endpoints 1 ',  &
                                err_PETSc)
  CHKERRQ(err_PETSc)
  call DMClone(geomMesh, coordDM, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetFromOptions(coordDM, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDimension(coordDM, coordDim, err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFECreateDefault(PETSC_COMM_SELF, coordDim, coordDim, PETSC_TRUE, 'coord_', p_s, &
                            coordFE, err_PETSc)
  CHKERRQ(err_PETSc)

  call PetscFEGetDimension(coordFE, feDim, err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MINOR>18 && PETSC_VERSION_MINOR<23)
  call DMAddField(coordDM, PETSC_NULL_DMLABEL, coordFE, err_PETSc)
#else
  call DMAddField(coordDM, PETSC_NULL_DMLABEL, PetscObjectCast(coordFE), err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  call DMCreateDS(coordDM, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDS(coordDM, coordDS, err_PETSc)
  CHKERRQ(err_PETSc)

  call DMGetGlobalSection(coordDM, globalSection, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLocalSection(coordDM, localSection, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateLocalVector(coordDM, coordVec, err_PETSc)
  CHKERRQ(err_PETSc)

  call PetscFEGetDualSpace(coordFE, coordDualSpace, err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFEDestroy(coordFE, err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexGetHeightStratum(coordDM, 0_pPETSCINT, cellStart, cellEnd, err_PETSc)
  CHKERRQ(err_PETSc)

  allocate(nodeCoords(coordDim))
  allocate(refCoords(feDim))
  allocate(mappedCoords(feDim))
  nCellNodes = feDim / coordDim

  do feBasis = 0_pPETSCINT, feDim-1_pPETSCINT, coordDim                                             ! coordinates in the reference cell in [-1,+1]^d
    call PetscDualSpaceGetFunctional(coordDualSpace, feBasis, refQuadrature, err_PETSc)
    CHKERRQ(err_PETSc)
    call PetscQuadratureGetData(refQuadrature, coordDim, coordDim, PETSC_NULL_INTEGER, &
                                nodeCoords, PETSC_NULL_REAL_POINTER, err_PETSc)
    CHKERRQ(err_PETSc)
    refCoords(feBasis + 1:feBasis + coordDim) = nodeCoords
  end do

  do cell = cellStart, cellEnd - 1_pPETSCINT                                                        ! map reference to real (mesh) coordinates
    call DMPlexReferenceToCoordinates(coordDM, cell, nCellNodes, refCoords, &
                                      mappedCoords, err_PETSc)
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MINOR>=23)
    PetscCall(DMPlexVecSetClosure(coordDM, localSection, coordVec, cell, mappedCoords,INSERT_VALUES, ierr))
#else
    pMappedCoords => mappedCoords
    call DMPlexVecSetClosure(coordDM, localSection, coordVec, cell, pMappedCoords, &
                             INSERT_VALUES, err_PETSc)
    CHKERRQ(err_PETSc)
#endif
  end do

  call VecGetArrayRead(coordVec,coords,err_PETSc)
  nLocalNodes = size(coords) / coordDim
  allocate(x_n(3, nLocalNodes), source = 0.0_pREAL)
  x_n(1:coordDim, 1:nLocalNodes) = reshape(coords, [coordDim, nLocalNodes])
  call VecRestoreArrayRead(coordVec, coords, err_PETSc)

#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  PetscCall(DMPlexGetHeightStratum(coordDM, 0_pPETSCINT, cellStart, cellEnd, ierr))
  PetscCall(DMPlexGetCellType(coordDM, cellStart, cell_type, ierr))
  node_map = PETSc_to_VTK_node_order(cell_type, p_s)

  allocate(connectivity_matrix(nCellNodes, cellEnd - cellStart), source = -1_pPETSCINT)
  do cell = cellStart, cellEnd - 1_pPETSCINT
    call DMPlexGetClosureIndices(coordDM, localSection, globalSection, cell, PETSC_TRUE, &
                                 PETSC_NULL_INTEGER, indices, PETSC_NULL_INTEGER_ARRAY, &
                                 PETSC_NULL_REAL_POINTER, ierr)
    CHKERRQ(ierr)
    connectivity_matrix(1:nCellNodes, cell + 1_pPETSCINT) = indices(coordDim * node_map - 1_pPETSCINT) &
                                                          / coordDim
    call DMPlexRestoreClosureIndices(coordDM, localSection, globalSection, cell, PETSC_TRUE, &
                                     PETSC_NULL_INTEGER, indices, PETSC_NULL_INTEGER_ARRAY, &
                                     PETSC_NULL_REAL_POINTER, ierr)
    CHKERRQ(ierr)
  end do
#endif

  call DMDestroy(coordDM, err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine build_nodes_and_connectivity


!--------------------------------------------------------------------------------------------------
!> @brief Map node order from PETSc to VTK.
!--------------------------------------------------------------------------------------------------
function PETSc_to_VTK_node_order(cell_type, order) result(mapping)

  integer, allocatable, dimension(:) :: &
    mapping                                                                                         !< node mapping
  DMPolytopeType, intent(in) :: &
    cell_type                                                                                       !< tri, quad, tet, hex
  integer,        intent(in) :: &
    order                                                                                           !< approximation space order


  if (cell_type == DM_POLYTOPE_TRIANGLE) then
    select case (order)
      case (1)
        mapping = [ 1,  2,  3]
      case (2)
        mapping = [ 4,  5,  6,  1,  2,  3]
      case (3)
        mapping = [ 8,  9, 10,  2,  3,  4,  5,  6,  7,  1]
      case (4)
        mapping = [13, 14, 15,  4,  5,  6,  7,  8,  9, 10, 11, 12,  1,  2,  3]
      case (5)
        mapping = [19, 20, 21,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1,  3,  6,  2,  5,  4]
    end select
  else if (cell_type == DM_POLYTOPE_TETRAHEDRON) then
    select case (order)
      case (1)
        mapping = [ 2,  4,  1,  3]
      case (2)
        mapping = [ 8, 10,  7,  9,  5,  4,  1,  2,  6,  3]
      case (3)
        mapping = [18, 20, 17, 19, 14, 13, 12, 11,  5,  6,  7,  8, 16, 15, 10,  9,  4,  3,  1,  2]
      case (4)
        mapping = [33, 35, 32, 34, 28, 27, 26, 25, 24, 23, 14, 15, 16, 17, 18, 19, 31, 30, 29, 22, &
                   21, 20, 12, 13, 11,  8,  9, 10,  3,  4, 2,  7,  5,  6,  1]
      case (5)
        mapping = [54, 56, 53, 55, 48, 47, 46, 45, 44, 43, 42, 41, 29, 30, 31, 32, 33, 34, 35, 36, &
                   52, 51, 50, 49, 40, 39, 38, 37, 25, 28, 23, 27, 26, 24, 17, 19, 22, 18, 21, 20, &
                    7, 10,  5,  9,  8,  6, 16, 11, 13, 14, 12, 15,  3,  4,  1,  2]
    end select
  end if

end function PETSc_to_VTK_node_order


!--------------------------------------------------------------------------------------------------
!> @brief Write all information needed for the DADF5 geometry.
!--------------------------------------------------------------------------------------------------
subroutine writeGeometry(coordinates_points,x_n,connectivity_matrix)

  real(pREAL), dimension(:,:), intent(in) :: &
    x_n, &                                                                                          !< mesh nodes coordinates (including high-order approximation)
    coordinates_points                                                                              !< ip coordinates
  PetscInt,    dimension(:,:), intent(in), optional :: &
    connectivity_matrix                                                                             !< element connectivity (node numbers in each cell)


  call result_openJobFile()
  call result_closeGroup(result_addGroup('geometry'))

  call result_writeDataset(x_n,'geometry','x_n', &
                           'initial coordinates of the nodes','m')

  call result_writeDataset(coordinates_points,'geometry','x_p', &
                           'initial coordinates of the materialpoints (cell centers)','m')

  if (present(connectivity_matrix)) &
    call result_writeDataset(int(connectivity_matrix),'geometry','T_e', &
                             'connectivity matrix','1')

  call result_closeJobFile()

end subroutine writeGeometry

end module discretization_mesh
