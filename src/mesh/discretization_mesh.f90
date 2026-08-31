! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Javier Velo, KU Leuven
!--------------------------------------------------------------------------------------------------
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdmda.h>
module discretization_mesh
  use PETScDMplex
  use PETScDMDA
  use PETScIS
#if PETSC_VERSION_MINOR<23
  use PETScDT
#endif
#ifndef PETSC_EXPOSES_MPI
  use MPI_f08
#endif

  use CLI
  use prec
  use parallelization
  use config
  use IO
  use discretization
  use result
  use types
  use prec

#ifndef PETSC_EXPOSES_MPIF90
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
    mesh_boundariesIdx                                                                              !< PETSC_BC_TYPE_X

  PetscInt, public, protected :: &
    mesh_nElems

!!!! BEGIN DEPRECATED !!!!!
  PetscInt, public, protected :: &
    mesh_maxNips                                                                                    !< max number of IPs in any CP element
!!!! END DEPRECATED !!!!!

  PetscInt, parameter, public :: &
    mesh_BCTypeFace   = 1_pPETSCINT, &
    mesh_BCTypeVertex = 2_pPETSCINT

  real(pREAL), dimension(:,:), allocatable, public, protected :: &
    x_n                                                                                             !< node x,z,y coordinates

  character(len=pSTRLEN), dimension(:), allocatable, public, protected :: &
    mesh_BCLabels

  enum, bind(c); enumerator :: &
    PETSC_BC_TYPE_CELL = 1, &
    PETSC_BC_TYPE_FACE, &
    PETSC_BC_TYPE_EDGE, &
    PETSC_BC_TYPE_VERTEX
  end enum

  character(len=*), &
  dimension(size([PETSC_BC_TYPE_CELL,PETSC_BC_TYPE_FACE,PETSC_BC_TYPE_EDGE,PETSC_BC_TYPE_VERTEX])), &
  parameter, public :: &
    PETSC_GENERIC_LABELS = ['Cell Sets  ', 'Face Sets  ', 'Edge Sets  ', 'Vertex Sets']             ! PETSc generic labels

  DM, public :: geom

#if PETSC_VERSION_MINOR<23
  external :: &
#if PETSC_VERSION_MINOR<22
    DMAddField, &
#endif
    PetscDualSpaceGetFunctional, &
    PetscFEDestroy, &
    PetscFEGetDimension, &
    PetscFEGetDualSpace, &
    PetscFESetQuadrature
#endif

  public :: &
    discretization_mesh_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize the mesh discretization.
!--------------------------------------------------------------------------------------------------
subroutine discretization_mesh_init()

  DM        :: globalMesh
  PetscInt  :: dimPlex, &                                                                           ! mesh dimension
               p_s, p_i, &                                                                          ! shape function/integration order
               cell_sets_size, &                                                                    ! size of 'Cell Sets' label
               n_mesh_labels, &                                                                     ! total number of labels in mesh file
               cell_start, cell_end, point_start, &
               j
  IS        :: label_values_IS                                                                      ! BC label values IS
  PetscBool :: has_label, &                                                                         ! label exists in the mesh
               is_simplex                                                                           ! simplex mesh
#if PETSC_VERSION_MINOR>=24
  IS        :: cell_types_IS                                                                        ! 'celltype' label IS
  PetscInt  :: n_polytopes                                                                          ! number of different polytopes in the mesh
  PetscInt,    dimension(:),     pointer     :: cell_types                                          ! cell_types_IS values
  PetscInt,    dimension(:,:),   allocatable :: T_e                                                 ! element connectivity (node numbers in each cell)
#else
  real(pREAL), dimension(:),     pointer     :: qPointsP                                            ! quadrature points coordinates
  real(pREAL), dimension(:),     pointer     :: PETSC_NULL_REAL_POINTER => NULL()
#endif
  DMLabel :: label
  PetscSF :: SF
  PetscDS :: global_DS
  PetscFE :: global_FE
  PetscQuadrature :: quadrature
  PetscErrorCode  :: err_PETSc

  type(tDict), pointer :: &
    num_solver, &
    num_mesh
  integer                                    :: dim, n, m, k
  integer(MPI_INTEGER_KIND)                  :: err_MPI
  integer,     dimension(:),     allocatable :: BC_set_idx                                          ! index for PETSc set labels (no 'edges' in 2D)
  PetscInt,    dimension(:),     pointer     :: label_values                                        ! BC label values (from IS)
  PetscInt,    dimension(:),     allocatable :: materialAt, &                                       ! material ID per cell
                                                label_tmp
  real(pREAL), dimension(:,:),   allocatable :: v_0                                                 ! volume associated with IP (initially!)
  real(pREAL), dimension(:,:,:), allocatable :: x_p                                                 ! IP x,y,z coordinates

  character(pSTRLEN)            :: BC_label                                                         ! label (string, defined in mesh file)
  character(len=:), allocatable :: PETSc_options                                                    ! options to set up DM (from numerics file)


  print'(/,1x,a)',   '<<<+-  discretization_mesh init  -+>>>'; flush(IO_STDOUT)

!--------------------------------------------------------------------------------------------------
! read numerics parameter
  num_solver => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_mesh   => num_solver%get_dict('mesh',defaultVal=emptyDict)
  p_i = int(num_mesh%get_asInt('p_i',defaultVal=2),pPETSCINT)
  p_s = int(num_mesh%get_asInt('p_s',defaultVal=2),pPETSCINT)

  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,                      &
                                ' -dm_plex_filename ' // CLI_geomFile // &
                                ' -dm_plex_interpolate 1                 &
                                & -dm_plex_gmsh_use_generic              &
                                & -dm_plex_gmsh_use_regions              &
                                & -dm_plex_gmsh_multiple_tags            &
                                & -dm_plex_gmsh_mark_vertices ',         &
                                err_PETSc)
  call DMCreate(PETSC_COMM_WORLD,globalMesh,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetType(globalMesh,DMPLEX,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetFromOptions(globalMesh,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDimension(globalMesh,dimPlex,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexIsSimplex(globalMesh,is_simplex,err_PETSc)
  CHKERRQ(err_PETSc)
  if (.not. is_simplex) then
#if PETSC_VERSION_MINOR<24
    call IO_error(800_pI16, 'mesh is not a simplex')
#endif
    p_i = p_i + 1_pPETSCINT                                                                         ! adjust for quad/hex (non-simplex)
  end if
  call DMGetStratumSize(globalMesh,'depth',dimPlex,mesh_nElems,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! check invalid mesh: empty 'Cell Sets' or mixed elements
  call DMGetLabelSize(globalMesh,'Cell Sets',cell_sets_size,err_PETSc)                              ! cannot be empty, needed to assign material ID
  if (cell_sets_size == 0_pPETSCINT) &
    call IO_error(800_pI16,'missing definition of ',&
                            trim(merge('surface','volume ',dimPlex == 2_pPETSCINT)),&
                            ' group(s) for assigning material IDs')
#if PETSC_VERSION_MINOR>=24
  call DMGetLabelIdIS(globalMesh,'celltype',cell_types_IS,err_PETSc)
  call ISGetSize(cell_types_IS,n_polytopes,err_PETSc)
  if (n_polytopes /= dimPlex + 1_pPETSCINT) then                                                    ! at most one polytope type per dimension (0..dimPlex)
    call ISGetIndices(cell_types_IS,cell_types,err_PETSc)
    if (any(cell_types == DM_POLYTOPE_SEG_PRISM_TENSOR%v .or. &
            cell_types >  DM_POLYTOPE_HEXAHEDRON%v)) then
      call IO_error(800_pI16,'mesh contains elements other than tri/quad/tet/hex')
    else if (count(cell_types == DM_POLYTOPE_TRIANGLE%v .and. &
                   cell_types == DM_POLYTOPE_QUADRILATERAL%v,dim = 1) > 1) then
      call IO_error(800_pI16,'mixed element types (triangle and quadrilateral)')
    else if (count(cell_types == DM_POLYTOPE_TETRAHEDRON%v .and. &
                   cell_types == DM_POLYTOPE_HEXAHEDRON%v,dim = 1) > 1) then
      call IO_error(800_pI16,'mixed element types (tetrahedron and hexahedron)')
    end if
  end if
  call ISDestroy(cell_types_IS,err_PETSc)
#endif

!--------------------------------------------------------------------------------------------------
! read mesh tags
  if (dimPlex == 2_pPETSCINT) then
    allocate(BC_set_idx,source = [PETSC_BC_TYPE_FACE, PETSC_BC_TYPE_VERTEX])
  else
    allocate(BC_set_idx,source = [PETSC_BC_TYPE_FACE, PETSC_BC_TYPE_EDGE, PETSC_BC_TYPE_VERTEX])
  end if

  mesh_Nboundaries = 0_pPETSCINT
  do n = 1, size(BC_set_idx)
    call DMHasLabel(globalMesh,PETSC_GENERIC_LABELS(BC_set_idx(n)),has_label,err_PETSc)
    if (has_label) then
      call DMGetLabel(globalMesh,PETSC_GENERIC_LABELS(BC_set_idx(n)),label,err_PETSc)
      CHKERRQ(err_PETSc)
      call DMLabelGetNonEmptyStratumValuesIS(label, label_values_IS, err_PETSc)
      CHKERRQ(err_PETSc)
      call ISGetIndices(label_values_IS,label_values,err_PETSc)
      CHKERRQ(err_PETSc)
      if (.not. allocated(mesh_boundariesIS)) then
        allocate(mesh_boundariesIS,source = label_values)
        allocate(mesh_boundariesIdx(size(label_values)),source = int(BC_set_idx(n),pPETSCINT))
        mesh_Nboundaries = mesh_Nboundaries + int(size(label_values),pPETSCINT)
      else
        allocate(label_tmp,mold = label_values)
        k = 0
        do m = 1, size(label_values)
          if (any(label_values(m) == mesh_boundariesIS)) cycle
          k = k + 1
          label_tmp(k) = label_values(m)
        end do
        mesh_boundariesIS = [mesh_boundariesIS,label_tmp(1:k)]
        mesh_boundariesIdx = [mesh_boundariesIdx,[(int(BC_set_idx(n),pPETSCINT),m = 1,k)]]
        mesh_Nboundaries = mesh_Nboundaries + int(k,pPETSCINT)
        deallocate(label_tmp)
      end if
      call ISRestoreIndices(label_values_IS,label_values,err_PETSc)
      CHKERRQ(err_PETSc)
    end if
  end do
  deallocate(BC_set_idx)
  if (mesh_Nboundaries == 0_pPETSCINT) &
    call IO_error(800_pI16,'no groups available for boundary condition assignment')

!--------------------------------------------------------------------------------------------------
! read mesh labels
  call DMGetNumLabels(globalMesh,n_mesh_labels,err_PETSc)
  CHKERRQ(err_PETSc)
  if (n_mesh_labels > 2_pPETSCINT) then                                                             ! there are user-defined labels (for BC/material ID)
    allocate(character(len=pSTRLEN) :: mesh_BCLabels(mesh_Nboundaries))
    mesh_BCLabels = ''

    call DMPlexGetHeightStratum(globalMesh,0_pPETSCINT,cell_start,cell_end,err_PETSc)
    do j = 2_pPETSCINT, n_mesh_labels - 1_pPETSCINT                                                 ! skip 'celltype' and 'depth' labels; 0-indexing in PETSc
      call DMGetLabelName(globalMesh,j,BC_label,err_PETSc)
      CHKERRQ(err_PETSc)
      if (any(BC_label == PETSC_GENERIC_LABELS)) cycle
      call DMGetLabel(globalMesh,BC_label,label,err_PETSc)
      call DMLabelGetBounds(label,point_start,PETSC_NULL_INTEGER,err_PETSc)
      if (point_start < cell_end) cycle
      call DMLabelGetNonEmptyStratumValuesIS(label,label_values_IS,err_PETSc)
      CHKERRQ(err_PETSc)
      call ISGetIndices(label_values_IS,label_values,err_PETSc)
      CHKERRQ(err_PETSc)
      n = findloc(mesh_boundariesIS,label_values(1),dim = 1)
      mesh_BCLabels(n) = BC_label
      call ISRestoreIndices(label_values_IS,label_values,err_PETSc)
      CHKERRQ(err_PETSc)
    end do
  else
    mesh_BCLabels = emptyStrArray
  end if

  dim = int(dimPlex)
  call MPI_Bcast(dim,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  dimPlex = int(dim,pPETSCINT)

  if (worldsize == 1) then
    call DMClone(globalMesh,geom,err_PETSc)
  else
    call DMPlexDistribute(globalMesh,0_pPETSCINT,SF,geom,err_PETSc)
  end if
  CHKERRQ(err_PETSc)

  call MPI_Bcast(mesh_Nboundaries,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(mesh_boundariesIS,int(mesh_Nboundaries,MPI_INTEGER_KIND),MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(mesh_boundariesIdx,int(mesh_Nboundaries,MPI_INTEGER_KIND),MPI_INTEGER,0_MPI_INTEGER_KIND, &
                 MPI_COMM_WORLD,err_MPI)
  if (worldrank /= 0) allocate(character(len=pSTRLEN) :: mesh_BCLabels(mesh_Nboundaries))
  call MPI_Bcast(mesh_BCLabels,int(pSTRLEN*size(mesh_BCLabels),MPI_INTEGER_KIND),MPI_CHARACTER, &
                 0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

!--------------------------------------------------------------------------------------------------
! set up quadrature, fields (mechanical, thermal...) and DS
#if PETSC_VERSION_MINOR>=24
  if (is_simplex) then
    call PetscDTSimplexQuadrature(dimPlex,p_i,PETSCDTSIMPLEXQUAD_DEFAULT,quadrature,err_PETSc)
  else
    call PetscDTGaussTensorQuadrature(dimPlex,dimPlex,p_i,-1.0_pREAL,1.0_pREAL, &
                                      quadrature,err_PETSc)
  end if
#elif PETSC_VERSION_MINOR==23
  call PetscDTSimplexQuadrature(dimPlex,p_i,PETSCDTSIMPLEXQUAD_DEFAULT,quadrature,err_PETSc)
#else
  call PetscDTSimplexQuadrature(dimPlex,p_i,-1,quadrature,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  PETSc_options = ' -petscspace_degree ' // IO_intAsStr(int(p_s)) // &
                  ' -petscdualspace_lagrange_node_type equispaced    &
                  & -petscdualspace_lagrange_node_endpoints 1     '
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMcreateFEDefault(geom,dimPlex,'',p_i,global_FE,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFESetQuadrature(global_FE,quadrature,err_PETSc)
  CHKERRQ(err_PETSc)
#if PETSC_VERSION_MINOR>22
  call DMAddField(geom,PETSC_NULL_DMLABEL,PetscObjectCast(global_FE),err_PETSc)
#else
  call DMAddField(geom,PETSC_NULL_DMLABEL,global_FE,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  call DMCreateDS(geom,err_PETSc)                                                                   ! create after adding all fields
  CHKERRQ(err_PETSc)
  call DMGetDS(geom, global_DS, err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscDSSetForceQuad(global_DS,PETSC_FALSE,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDestroy(globalMesh,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! set up geometry values (node coordinates, IP, volume, connectivity, material ID)
#if PETSC_VERSION_MINOR>=24
  call PetscQuadratureGetData(quadrature,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                              mesh_maxNips,PETSC_NULL_REAL_POINTER, &
                              PETSC_NULL_REAL_POINTER,err_PETSc)
#elif PETSC_VERSION_MINOR>21
  call PetscQuadratureGetData(quadrature,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                              mesh_maxNips,qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#else
  call PetscQuadratureGetData(quadrature,PETSC_NULL_INTEGER(1),PETSC_NULL_INTEGER(1), &
                              mesh_maxNips,qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

#if PETSC_VERSION_MINOR>=24
  x_p = build_coordinates_IP(dimPlex,quadrature)
#else
  x_p = build_coordinates_IP(dimPlex,qPointsP)
#endif
  v_0 = build_volume_IP(dimPlex)

#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  call build_nodes_and_connectivity(x_n,T_e,p_s)
#else
  call build_nodes_and_connectivity(x_n,p_s)
#endif

#if (PETSC_VERSION_MINOR==22 || PETSC_VERSION_MINOR==23)
  call PetscQuadratureRestoreData(quadrature,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                                  PETSC_NULL_INTEGER,qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#elif (PETSC_VERSION_MINOR<22)
  call PetscQuadratureRestoreData(quadrature,PETSC_NULL_INTEGER(1),PETSC_NULL_INTEGER(1), &
                                  PETSC_NULL_INTEGER(1),qPointsP,PETSC_NULL_REAL_POINTER,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  call PetscQuadratureDestroy(quadrature,err_PETSc)
  CHKERRQ(err_PETSc)

  allocate(materialAt(mesh_nElems))
  do j = 1_pPETSCINT, mesh_nElems
    call DMGetLabelValue(geom,'Cell Sets',j-1_pPETSCINT,materialAt(j),err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  call discretization_init(int(materialAt),reshape(x_p,[3,int(mesh_maxNips*mesh_nElems)]),x_n)

#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  call writeGeometry(reshape(x_p,[3,int(mesh_maxNips*mesh_nElems)]),x_n,T_e)
#else
  call writeGeometry(reshape(x_p,[3,int(mesh_maxNips*mesh_nElems)]),x_n)
#endif

end subroutine discretization_mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP volume.
!--------------------------------------------------------------------------------------------------
function build_volume_IP(dimPlex) result(v_0)

  real(pREAL), dimension(:,:), allocatable :: v_0
  PetscInt,    intent(in) :: dimPlex

  PetscReal      :: vol
  PetscInt       :: cell_start, cell_end, cell
  PetscErrorCode :: err_PETSc
  PetscReal, pointer,dimension(:) :: pCent, pNorm


  allocate(v_0(mesh_maxNips,mesh_nElems),source=0.0_pREAL)

  call DMPlexGetHeightStratum(geom,0_pPETSCINT,cell_start,cell_end,err_PETSc)
  CHKERRQ(err_PETSc)
  allocate(pCent(dimPlex))
  allocate(pNorm(dimPlex))
  do cell = cell_start, cell_end - 1_pPETSCINT
    call DMPlexComputeCellGeometryFVM(geom,cell,vol,pCent,pNorm,err_PETSc)
    CHKERRQ(err_PETSc)
    v_0(:,cell+1) = vol/real(mesh_maxNips,pREAL)
  end do

end function build_volume_IP


!--------------------------------------------------------------------------------------------------
!> @brief Calculate IP Coordinates.
!--------------------------------------------------------------------------------------------------
#if (PETSC_VERSION_MINOR>=24)
function build_coordinates_IP(dimPlex,quadrature) result(x_p)
#else
function build_coordinates_IP(dimPlex,qPoints) result(x_p)
#endif

  PetscInt,                                     intent(in) :: dimPlex
#if PETSC_VERSION_MINOR>=24
  PetscQuadrature,                              intent(in) :: quadrature
#else
  PetscReal,   dimension(mesh_maxNips*dimPlex), intent(in) :: qPoints
#endif
  real(pREAL), dimension(:,:,:), allocatable :: x_p


#if PETSC_VERSION_MINOR>=24
  PetscReal, pointer, dimension(:) :: pV0, pCellJ, pInvcellJ, pDetJ
#else
  PetscReal, pointer, dimension(:) :: pV0, pCellJ, pInvcellJ
  PetscReal      :: detJ
  PetscInt       :: qPt, dirI, dirJ, qOffset
#endif
  PetscInt       :: cell_start, cell_end, cell
  PetscErrorCode :: err_PETSc


  call DMPlexGetHeightStratum(geom,0_pPETSCINT,cell_start,cell_end,err_PETSc)
  CHKERRQ(err_PETSc)

  allocate(x_p(3,mesh_maxNips,mesh_nElems),source=0.0_pREAL)
#if PETSC_VERSION_MINOR>=24
  allocate(pV0(mesh_maxNips*dimPlex))
  allocate(pCellJ(mesh_maxNips*dimPlex**2))
  allocate(pInvCellJ(mesh_maxNips*dimPlex**2))
  allocate(pDetJ(mesh_maxNips))

  do cell = cell_start, cell_end - 1_pPETSCINT
    call DMPlexComputeCellGeometryFEM(geom,cell,quadrature,pV0,pCellJ,pInvCellJ,pDetJ,err_PETSc)
    CHKERRQ(err_PETSc)
    x_p(1:dimPlex,1:mesh_maxNips,cell+1_pPETSCINT) = reshape(pV0,[dimPlex,mesh_maxNips])
  end do
#else
  allocate(pV0(dimPlex))
  allocate(pCellJ(dimPlex**2))
  allocate(pinvCellJ(dimPlex**2))

  do cell = cell_start, cell_end - 1_pPETSCINT                                                      ! loop over all elements
    call DMPlexComputeCellGeometryAffineFEM(geom,cell,pV0,pCellJ,pInvcellJ,detJ,err_PETSc)
    CHKERRQ(err_PETSc)
    qOffset = 0_pPETSCINT
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
#endif

end function build_coordinates_IP


!--------------------------------------------------------------------------------------------------
!> @brief Get mesh node coordinates and build element connectivity matrix.
!--------------------------------------------------------------------------------------------------
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
subroutine build_nodes_and_connectivity(x_n, T_e, p_s)
#else
subroutine build_nodes_and_connectivity(x_n, p_s)
#endif

  real(pREAL), dimension(:,:), allocatable, intent(out) :: x_n                                      !< mesh nodes coordinates (including high-order approximation)
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  PetscInt,    dimension(:,:), allocatable, intent(out) :: T_e                                      !< element connectivity (node numbers in each cell)
#endif
  PetscInt,                                 intent(in)  :: p_s                                      !< order of approximation space

  Vec      :: coords_vec                                                                            ! local nodes coordinates
  PetscInt :: dimPlex, &                                                                            ! DoF per node
              n_nodes, &                                                                            ! (local) number of nodes
              n_cell_nodes, &                                                                       ! number of nodes per cell
              FE_dim, &                                                                             ! DoF per cell (n_nodes x DoF per node)
              cell_start, cell_end, &                                                               ! (local) first & last+1 cell number
              cell, basis                                                                           ! loop
  PetscDS  :: DS
  PetscFE  :: FE
#if PETSC_VERSION_MINOR>22
  PetscObject     :: FE_obj
#endif
  PetscSection    :: global_section, local_section                                                  ! section (to retrieve DoF)
#if PETSC_VERSION_MINOR>=24
  DMPolytopeType  :: cell_type                                                                      ! tri, quad, tet, hex
#endif
  PetscDualSpace  :: dual_space                                                                     ! FE dual space
  PetscQuadrature :: quadrature                                                                     ! functional from dual space
  PetscErrorCode  :: err_PETSc
  real(pREAL), dimension(:), pointer     :: coords, &                                               ! local nodes coordinates
                                            node_coords                                             ! single node coordinates
  real(pREAL), dimension(:), allocatable :: ref_coords                                              ! node coordinates in reference element [-1,+1]^d
#if PETSC_VERSION_MINOR>22
  real(pREAL), dimension(:), allocatable :: mapped_coords                                           ! real (mesh) node coordinates
#else
  real(pREAL), dimension(:), allocatable, target :: mapped_coords                                   ! real (mesh) node coordinates
  real(pREAL), dimension(:), pointer     :: mapped_coords_, &                                       ! pointer to mapped_coords
                                            PETSC_NULL_REAL_POINTER => NULL()
#endif
#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  PetscInt,    dimension(:), pointer     :: indices                                                 ! cell closure DoF indices
  integer,     dimension(:), allocatable :: node_map                                                ! PETSc to VTK node order mapping
#endif


  call DMGetDimension(geom,dimPlex,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDS(geom,DS,err_PETSc)
  CHKERRQ(err_PETSc)
#if PETSC_VERSION_MINOR>22
  call PetscDSGetDiscretization(DS,0_pPETSCINT,FE_obj,err_PETSc)
  CHKERRQ(err_PETSc)
  PetscObjectSpecificCast(FE,FE_obj)
#else
  call PetscDSGetDiscretization(DS,0_pPETSCINT,FE,err_PETSc)
#endif
  call PetscFEGetDimension(FE,FE_dim,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFEGetDualSpace(FE,dual_space,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFEDestroy(FE,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMGetGlobalSection(geom,global_section,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLocalSection(geom,local_section,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateLocalVector(geom,coords_vec,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexGetHeightStratum(geom,0_pPETSCINT,cell_start,cell_end,err_PETSc)
  CHKERRQ(err_PETSc)

  allocate(node_coords(dimPlex))
  allocate(ref_coords(FE_dim))
  do basis = 0_pPETSCINT, FE_dim - 1_pPETSCINT, dimPlex                                             ! coordinates in the reference cell in [-1,+1]^d
    call PetscDualSpaceGetFunctional(dual_space,basis,quadrature,err_PETSc)
    CHKERRQ(err_PETSc)
#if PETSC_VERSION_MINOR>21
    call PetscQuadratureGetData(quadrature,dimPlex,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
#else
    call PetscQuadratureGetData(quadrature,dimPlex,PETSC_NULL_INTEGER(1), &
                                PETSC_NULL_INTEGER(1),&
#endif
                                node_coords,PETSC_NULL_REAL_POINTER,err_PETSc)
    CHKERRQ(err_PETSc)
    ref_coords(basis+1_pPETSCINT:basis+dimPlex) = node_coords
  end do

  allocate(mapped_coords(FE_dim))
  n_cell_nodes = FE_dim / dimPlex
  do cell = cell_start, cell_end - 1_pPETSCINT                                                      ! map reference to real (mesh) coordinates
    call DMPlexReferenceToCoordinates(geom,cell,n_cell_nodes,ref_coords,&
                                      mapped_coords,err_PETSc)
    CHKERRQ(err_PETSc)
#if PETSC_VERSION_MINOR>22
    call DMPlexVecSetClosure(geom,local_section,coords_vec,cell,mapped_coords, &
#else
    mapped_coords_ => mapped_coords
    call DMPlexVecSetClosure(geom,local_section,coords_vec,cell,mapped_coords_, &
#endif
                             INSERT_VALUES,err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  call VecGetArrayRead(coords_vec,coords,err_PETSc)
  n_nodes = size(coords)/dimPlex
  allocate(x_n(3,n_nodes),source = 0.0_pREAL)
  x_n(1:dimPlex,1:n_nodes) = reshape(coords,[dimPlex,n_nodes])
  call VecRestoreArrayRead(coords_vec,coords,err_PETSc)

#if (PETSC_VERSION_MINOR>24 || (PETSC_VERSION_MINOR==24 && PETSC_VERSION_SUBMINOR>=1))
  call DMPlexGetHeightStratum(geom,0_pPETSCINT,cell_start,cell_end,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexGetCellType(geom,cell_start,cell_type,err_PETSc)
  CHKERRQ(err_PETSc)
  node_map = PETSc_to_VTK_node_order(cell_type,int(p_s))

  allocate(T_e(n_cell_nodes, cell_end - cell_start), source = -1_pPETSCINT)
  do cell = cell_start,cell_end - 1_pPETSCINT
    call DMPlexGetClosureIndices(geom,local_section,global_section,cell,PETSC_TRUE, &
                                 PETSC_NULL_INTEGER,indices,PETSC_NULL_INTEGER_ARRAY, &
                                 PETSC_NULL_REAL_POINTER,err_PETSc)
    CHKERRQ(err_PETSc)
    T_e(1_pPETSCINT:n_cell_nodes,cell+1_pPETSCINT) = indices(dimPlex*node_map-1_pPETSCINT)/dimPlex
    call DMPlexRestoreClosureIndices(geom,local_section,global_section,cell,PETSC_TRUE, &
                                     PETSC_NULL_INTEGER,indices,PETSC_NULL_INTEGER_ARRAY, &
                                     PETSC_NULL_REAL_POINTER,err_PETSc)
    CHKERRQ(err_PETSc)
  end do
#endif

end subroutine build_nodes_and_connectivity


!--------------------------------------------------------------------------------------------------
!> @brief Map node order from PETSc to VTK.
!--------------------------------------------------------------------------------------------------
function PETSc_to_VTK_node_order(cell_type,order) result(mapping)

  integer, allocatable, dimension(:) :: &
    mapping                                                                                         !< node mapping
  DMPolytopeType, intent(in) :: &
    cell_type                                                                                       !< tri, quad, tet, hex
  integer,        intent(in) :: &
    order                                                                                           !< approximation space order


  if (cell_type == DM_POLYTOPE_TRIANGLE) then
    select case (order)
      case (1)
        mapping = [  1,   2,   3]
      case (2)
        mapping = [  4,   5,   6,   1,   2,   3]
      case (3)
        mapping = [  8,   9,  10,   2,   3,   4,   5,   6,   7,   1]
      case (4)
        mapping = [ 13,  14,  15,   4,   5,   6,   7,   8,   9,  10,  11,  12,   1,   2,   3]
      case (5)
        mapping = [ 19,  20,  21,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  1,  &
                     3,   6,   2,   5,   4]
    end select
  else if (cell_type == DM_POLYTOPE_QUADRILATERAL) then
    select case (order)
      case (1)
        mapping = [  1,   2,   3,   4]
      case (2)
        mapping = [  6,   7,   8,   9,   2,   3,   4,   5,   1]
      case (3)
        mapping = [  13, 14,  15,  16,   5,   6,   7,   8,  10,  9,  12,   11,  1,   2,   3,   4]
      case (4)
        mapping = [ 22,  23,  24,  25,  10,  11,  12,  13,  14,  15,  18,  17,  16,  21,  20,  19, &
                     1,   2,   3,   4,   5,   6,   7,   8,   9]
      case (5)
        mapping = [ 33,  34,  35,  36,  17,  18,  19,  20,  21,  22,  23,  24,  28,  27,  26,  25, &
                    32,  31,  30,  29,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12, &
                    13,  14,  15,  16]
    end select
  else if (cell_type == DM_POLYTOPE_TETRAHEDRON) then
    select case (order)
      case (1)
        mapping = [  2,   4,   1,   3]
      case (2)
        mapping = [  8,  10,   7,   9,   5,   4,   1,   2,   6,   3]
      case (3)
        mapping = [ 18,  20,  17,  19,  14,  13,  12,  11,   5,   6,   7,   8,  16,  15,  10,   9, &
                     4,   3,   1,   2]
      case (4)
        mapping = [ 33,  35,  32,  34,  28,  27,  26,  25,  24,  23,  14,  15,  16,  17,  18,  19, &
                    31,  30,  29,  22,  21,  20,  12,  13,  11,   8,   9,  10,   3,   4,   2,   7, &
                     5,   6,   1]
      case (5)
        mapping = [ 54,  56,  53,  55,  48,  47,  46,  45,  44,  43,  42,  41,  29,  30,  31,  32, &
                    33,  34,  35,  36,  52,  51,  50,  49,  40,  39,  38,  37,  25,  28,  23,  27, &
                    26,  24,  17,  19,  22,  18,  21,  20,   7,  10,   5,   9,   8,   6,  16,  11, &
                    13,  14,  12,  15,   3,   4,   1,   2]
    end select
  else if (cell_type == DM_POLYTOPE_HEXAHEDRON) then
    select case (order)
      case (1)
        mapping = [  1,   4,   3,   2,   5,   6,   7,   8]
      case (2)
        mapping = [ 20,  23,  22,  21,  24,  25,  26,  27,  11,  10,   9,   8,  12,  13,  14,  15, &
                    17,  16,  19,  18,   7,   6,   4,   5,   2,   3,   1]
      case (3)
        mapping = [ 57,  60,  59,  58,  61,  62,  63,  64,  40,  39,  38,  37,  35,  36,  33,  34, &
                    41,  42,  43,  44,  46,  45,  48,  47,  52,  51,  49,  50,  56,  55,  53,  54, &
                    29,  31,  30,  32,  25,  26,  27,  28,  17,  18,  19,  20,  22,  21,  24,  23, &
                     9,  11,  10,  12,  13,  14,  15,  16,   1,   2,   3,   4,   5,  6,   7,   8]
      case (4)
        mapping = [118, 121, 120, 119, 122, 123, 124, 125,  93,  92,  91,  90,  89,  88,  85,  86, &
                    87,  82,  83,  84,  94,  95,  96,  97,  98,  99, 102, 101, 100, 105, 104, 103, &
                   111, 110, 109, 106, 107, 108, 117, 116, 115, 112, 113, 114,  73,  76,  79,  74, &
                    77,  80,  75,  78,  81,  64,  65,  66,  67,  68,  69,  70,  71,  72,  46,  47, &
                    48,  49,  50,  51,  52,  53,  54,  57,  56,  55,  60,  59,  58,  63,  62,  61, &
                    28,  31,  34,  29,  32,  35,  30,  33,  36,  37,  38,  39,  40,  41,  42,  43, &
                    44,  45,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14, &
                    15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27]
      case (5)
        mapping = [209, 212, 211, 210, 213, 214, 215, 216, 176, 175, 174, 173, 172, 171, 170, 169, &
                   165, 166, 167, 168, 161, 162, 163, 164, 177, 178, 179, 180, 181, 182, 183, 184, &
                   188, 187, 186, 185, 192, 191, 190, 189, 200, 199, 198, 197, 193, 194, 195, 196, &
                   208, 207, 206, 205, 201, 202, 203, 204, 145, 149, 153, 157, 146, 150, 154, 158, &
                   147, 151, 155, 159, 148, 152, 156, 160, 129, 130, 131, 132, 133, 134, 135, 136, &
                   137, 138, 139, 140, 141, 142, 143, 144,  97,  98,  99, 100, 101, 102, 103, 104, &
                   105, 106, 107, 108, 109, 110, 111, 112, 116, 115, 114, 113, 120, 119, 118, 117, &
                   124, 123, 122, 121, 128, 127, 126, 125,  65,  69,  73,  77,  66,  70,  74,  78, &
                    67,  71,  75,  79,  68,  72,  76,  80,  81,  82,  83,  84,  85,  86,  87,  88, &
                    89,  90,  91,  92,  93,  94,  95,  96,   1,   2,   3,   4,   5,   6,   7,   8, &
                     9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24, &
                    25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                    41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56, &
                    57,  58,  59,  60,  61,  62,  63,  64]
    end select
  end if

end function PETSc_to_VTK_node_order


!--------------------------------------------------------------------------------------------------
!> @brief Write all information needed for the DADF5 geometry.
!--------------------------------------------------------------------------------------------------
subroutine writeGeometry(x_p,x_n,T_e)

  real(pREAL), dimension(:,:), intent(in) :: &
    x_n, &                                                                                          !< mesh nodes coordinates (including high-order approximation)
    x_p                                                                                             !< ip coordinates
  PetscInt,    dimension(:,:), intent(in), optional :: &
    T_e                                                                                             !< element connectivity (node numbers in each cell)


  call result_openJobFile()
  call result_closeGroup(result_addGroup('geometry'))

  call result_writeDataset(x_n,'geometry','x_n','initial coordinates of the nodes','m')

  call result_writeDataset(x_p,'geometry','x_p', &
                           'initial coordinates of the materialpoints (cell centers)','m')

  if (present(T_e)) &
    call result_writeDataset(int(T_e),'geometry','T_e','connectivity matrix','1')

  call result_closeJobFile()

end subroutine writeGeometry

end module discretization_mesh
