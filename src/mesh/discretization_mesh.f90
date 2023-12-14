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
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  use CLI
  use parallelization
  use IO
  use config
  use discretization
  use result
  use FEM_quadrature
  use types
  use prec

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
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

  DM, public :: geomMesh

  PetscInt, dimension(:), allocatable, public, protected :: &
    mesh_boundaries

  real(pREAL), dimension(:,:), allocatable :: &
    mesh_ipVolume, &                                                                                !< volume associated with IP (initially!)
    mesh_node0                                                                                      !< node x,y,z coordinates (initially!)

  real(pREAL), dimension(:,:,:), allocatable :: &
    mesh_ipCoordinates                                                                              !< IP x,y,z coordinates (after deformation!)

#if defined(PETSC_USE_64BIT_INDICES) || PETSC_VERSION_MINOR < 17
  external :: &
    DMDestroy
#endif
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

  PetscInt :: dimPlex, &
    mesh_Nnodes, &                                                                                  !< total number of nodes in mesh
    j
  PetscSF :: sf
  DM :: globalMesh
  PetscInt :: nFaceSets, Nboundaries, NelemsGlobal, Nelems
  PetscInt, pointer, dimension(:) :: pFaceSets
  IS :: faceSetIS
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscInt, dimension(:), allocatable :: &
    materialAt
  type(tDict), pointer :: &
    num_solver, &
    num_mesh
  integer :: p_i, dim                                                                               !< integration order (quadrature rule)
  type(tvec) :: coords_node0
  real(pREAL), pointer, dimension(:) :: &
    mesh_node0_temp

  print'(/,1x,a)',   '<<<+-  discretization_mesh init  -+>>>'

!--------------------------------------------------------------------------------
! read numerics parameter
  num_solver => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_mesh   => num_solver%get_dict('mesh',defaultVal=emptyDict)
  p_i = num_mesh%get_asInt('p_i',defaultVal=2)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>16)
  call DMPlexCreateFromFile(PETSC_COMM_WORLD,CLI_geomFile,'n/a',PETSC_TRUE,globalMesh,err_PETSc)
#else
  call DMPlexCreateFromFile(PETSC_COMM_WORLD,CLI_geomFile,PETSC_TRUE,globalMesh,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
  call DMGetDimension(globalMesh,dimPlex,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetStratumSize(globalMesh,'depth',dimPlex,NelemsGlobal,err_PETSc)
  CHKERRQ(err_PETSc)
  mesh_NcpElemsGlobal = int(NelemsGlobal)

  ! get number of IDs in face sets (for boundary conditions?)
  call DMGetLabelSize(globalMesh,'Face Sets',Nboundaries,err_PETSc)
  CHKERRQ(err_PETSc)
  mesh_Nboundaries = int(Nboundaries)
  call MPI_Bcast(mesh_Nboundaries,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(mesh_NcpElemsGlobal,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  dim = int(dimPlex)
  call MPI_Bcast(dim,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  dimPlex = int(dim,pPETSCINT)
  call parallelization_chkerr(err_MPI)

  if (worldsize == 1) then
    call DMClone(globalMesh,geomMesh,err_PETSc)
  else
    call DMPlexDistribute(globalMesh,0_pPETSCINT,sf,geomMesh,err_PETSc)
  end if
  CHKERRQ(err_PETSc)

  allocate(mesh_boundaries(mesh_Nboundaries), source = 0_pPETSCINT)
  call DMGetLabelSize(globalMesh,'Face Sets',nFaceSets,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLabelIdIS(globalMesh,'Face Sets',faceSetIS,err_PETSc)
  CHKERRQ(err_PETSc)
  if (nFaceSets > 0) then
    call ISGetIndicesF90(faceSetIS,pFaceSets,err_PETSc)
    CHKERRQ(err_PETSc)
    mesh_boundaries(1:nFaceSets) = pFaceSets
    CHKERRQ(err_PETSc)
    call ISRestoreIndicesF90(faceSetIS,pFaceSets,err_PETSc)
  end if
  call MPI_Bcast(mesh_boundaries,mesh_Nboundaries,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  call DMDestroy(globalMesh,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMGetStratumSize(geomMesh,'depth',dimPlex,Nelems,err_PETSc)
  CHKERRQ(err_PETSc)
  mesh_NcpElems = int(Nelems)
  call DMGetStratumSize(geomMesh,'depth',0_pPETSCINT,mesh_Nnodes,err_PETSc)
  CHKERRQ(err_PETSc)

! Get initial nodal coordinates
  call DMGetCoordinatesLocal(geomMesh,coords_node0,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecGetArrayF90(coords_node0, mesh_node0_temp,err_PETSc)
  CHKERRQ(err_PETSc)

  mesh_maxNips = FEM_nQuadrature(dimPlex,p_i)

  call mesh_FEM_build_ipCoordinates(dimPlex,FEM_quadrature_points(dimPlex,p_i)%p)
  call mesh_FEM_build_ipVolumes(dimPlex)

  allocate(materialAt(mesh_NcpElems))
  do j = 1, mesh_NcpElems
    call DMGetLabelValue(geomMesh,'Cell Sets',j-1,materialAt(j),err_PETSc)
    CHKERRQ(err_PETSc)
  end do
  materialAt = materialAt + 1_pPETSCINT

  allocate(mesh_node0(3,mesh_Nnodes),source=0.0_pREAL)
  mesh_node0(1:dimPlex,:) = reshape(mesh_node0_temp,[dimPlex,mesh_Nnodes])


  call discretization_init(int(materialAt),&
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
  PetscErrorCode     :: err_PETSc

  allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems),source=0.0_pREAL)

  call DMPlexGetHeightStratum(geomMesh,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
  allocate(pCent(dimPlex))
  allocate(pNorm(dimPlex))
  do cell = cellStart, cellEnd-1
    call  DMPlexComputeCellGeometryFVM(geomMesh,cell,vol,pCent,pNorm,err_PETSc)
    CHKERRQ(err_PETSc)
    mesh_ipVolume(:,cell+1) = vol/real(mesh_maxNips,pREAL)
  end do

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
  PetscErrorCode            :: err_PETSc


  allocate(mesh_ipCoordinates(3,mesh_maxNips,mesh_NcpElems),source=0.0_pREAL)

  allocate(pV0(dimPlex))
  allocatE(pCellJ(dimPlex**2))
  allocatE(pinvCellJ(dimPlex**2))
  call DMPlexGetHeightStratum(geomMesh,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
  do cell = cellStart, cellEnd-1                                                                     !< loop over all elements
    call DMPlexComputeCellGeometryAffineFEM(geomMesh,cell,pV0,pCellJ,pInvcellJ,detJ,err_PETSc)
    CHKERRQ(err_PETSc)
    qOffset = 0
    do qPt = 1_pPETSCINT, mesh_maxNips
      do dirI = 1_pPETSCINT, dimPlex
        mesh_ipCoordinates(dirI,qPt,cell+1) = pV0(dirI)
        do dirJ = 1_pPETSCINT, dimPlex
          mesh_ipCoordinates(dirI,qPt,cell+1) = mesh_ipCoordinates(dirI,qPt,cell+1) + &
                                                pCellJ((dirI-1)*dimPlex+dirJ)*(qPoints(qOffset+dirJ) + 1.0_pREAL)
        end do
      end do
      qOffset = qOffset + dimPlex
    end do
  end do

end subroutine mesh_FEM_build_ipCoordinates

!--------------------------------------------------------------------------------------------------
!> @brief Write all information needed for the DADF5 geometry
!--------------------------------------------------------------------------------------------------
subroutine writeGeometry(coordinates_points,coordinates_nodes)

  real(pREAL), dimension(:,:), intent(in) :: &
  coordinates_nodes, &
  coordinates_points

  call result_openJobFile()
  call result_closeGroup(result_addGroup('geometry'))

  call result_writeDataset(coordinates_nodes,'geometry','x_n', &
                           'initial coordinates of the nodes','m')

  call result_writeDataset(coordinates_points,'geometry','x_p', &
                           'initial coordinates of the materialpoints (cell centers)','m')

  call result_closeJobFile()

  end subroutine writeGeometry

end module discretization_mesh
