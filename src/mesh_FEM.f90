!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module mesh     
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
 use FEM_Zoo
 use prec
 use mesh_base
 
 implicit none
 private
 
 integer, public, protected :: &
   mesh_Nboundaries, &
   mesh_NcpElems, &                                                                                 !< total number of CP elements in mesh
   mesh_NcpElemsGlobal, &
   mesh_Nnodes                                                                                      !< total number of nodes in mesh

!!!! BEGIN DEPRECATED !!!!!
 integer, public, protected :: &
   mesh_maxNips                                                                                     !< max number of IPs in any CP element
!!!! BEGIN DEPRECATED !!!!!

 integer, dimension(:,:), allocatable :: &
   mesh_element !DEPRECATED

 real(pReal), dimension(:,:), allocatable  :: &
   mesh_node                                                                                        !< node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)
 
 real(pReal), dimension(:,:), allocatable :: &
   mesh_ipVolume, &                                                                                 !< volume associated with IP (initially!)
   mesh_node0                                                                                       !< node x,y,z coordinates (initially!)
 
 real(pReal), dimension(:,:,:), allocatable, public :: &
   mesh_ipCoordinates                                                                               !< IP x,y,z coordinates (after deformation!)

 DM, public :: geomMesh
 
 PetscInt, dimension(:), allocatable, public, protected :: &
   mesh_boundaries

                       
  type, public, extends(tMesh) :: tMesh_FEM

   
   contains
   procedure, pass(self) :: tMesh_FEM_init
   generic, public :: init => tMesh_FEM_init
 end type tMesh_FEM
 
 type(tMesh_FEM), public, protected :: theMesh
 

 public :: &
   mesh_init, &
   mesh_FEM_build_ipVolumes, &
   mesh_FEM_build_ipCoordinates, &
   mesh_cellCenterCoordinates

contains

subroutine tMesh_FEM_init(self,dimen,order,nodes)
 
 integer, intent(in) :: dimen
 integer, intent(in) :: order
 real(pReal), intent(in), dimension(:,:) :: nodes
 class(tMesh_FEM) :: self
 
 if (dimen == 2) then
   if (order == 1) call  self%tMesh%init('mesh',1,nodes)
   if (order == 2) call  self%tMesh%init('mesh',2,nodes)
 elseif(dimen == 3) then
   if (order == 1) call self%tMesh%init('mesh',6,nodes)
   if (order == 2) call self%tMesh%init('mesh',8,nodes)
 endif

 end subroutine tMesh_FEM_init



!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init

 integer, dimension(1), parameter:: FE_geomtype = [1]                      !< geometry type of particular element type

 integer, dimension(1) :: FE_Nips                         !< number of IPs in a specific type of element

 
 integer, parameter :: FILEUNIT = 222
 integer :: j
 integer, allocatable, dimension(:) :: chunkPos
 integer :: dimPlex
 integer, parameter :: &
   mesh_ElemType=1                                                                             !< Element type of the mesh (only support homogeneous meshes)
 character(len=512) :: &
   line
 logical :: flag
 PetscSF :: sf
 DM :: globalMesh
 PetscInt :: face, nFaceSets
 PetscInt, pointer :: pFaceSets(:)
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
 if (nFaceSets > 0) call ISGetIndicesF90(faceSetIS,pFaceSets,ierr)
 do face = 1, nFaceSets
   mesh_boundaries(face) = pFaceSets(face)
 enddo  
 if (nFaceSets > 0) call ISRestoreIndicesF90(faceSetIS,pFaceSets,ierr)
 call MPI_Bcast(mesh_boundaries,mesh_Nboundaries,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)

 ! this read in function should ignore C and C++ style comments
 ! it is used for BC only?
 if (worldrank == 0) then
   j = 0
   flag = .false.
   call IO_open_file(FILEUNIT,trim(geometryFile))
   do
     read(FILEUNIT,'(a512)') line
     if (trim(line) == IO_EOF) exit                                                                    ! skip empty lines
     if (trim(line) == '$Elements') then
       read(FILEUNIT,'(a512)') line ! number of elements (ignore)
       read(FILEUNIT,'(a512)') line 
       flag = .true.  
     endif
     if (trim(line) == '$EndElements') exit
     if (flag) then
       chunkPos = IO_stringPos(line)
       if (chunkPos(1) == 3+IO_intValue(line,chunkPos,3)+dimPlex+1) then
         call DMSetLabelValue(globalMesh,'material',j,IO_intValue(line,chunkPos,4),ierr)
         CHKERRQ(ierr)
         j = j + 1  
       endif                                                                                           ! count all identifiers to allocate memory and do sanity check
     endif
   enddo
   close (FILEUNIT)
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

 FE_Nips(FE_geomtype(1)) = FEM_Zoo_nQuadrature(dimPlex,integrationOrder)
 mesh_maxNips = FE_Nips(1)
 
 write(6,*) 'mesh_maxNips',mesh_maxNips
 call mesh_FEM_build_ipCoordinates(dimPlex,FEM_Zoo_QuadraturePoints(dimPlex,integrationOrder)%p)
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
 if (allocated(FEsolving_execIP)) deallocate(FEsolving_execIP)
 allocate(FEsolving_execIP(2,mesh_NcpElems)); FEsolving_execIP = 1                        ! parallel loop bounds set to comprise from first IP...
 forall (j = 1:mesh_NcpElems) FEsolving_execIP(2,j) = FE_Nips(FE_geomtype(mesh_element(2,j)))  ! ...up to own IP count for each element
 
 allocate(mesh_node0(3,mesh_Nnodes),source=0.0_pReal)
 call theMesh%init(dimplex,integrationOrder,mesh_node0)
 call theMesh%setNelems(mesh_NcpElems)

   call discretization_init(mesh_element(3,:),mesh_element(4,:),&
                           reshape(mesh_ipCoordinates,[3,mesh_maxNips*mesh_NcpElems]), &
                           mesh_node0)
 
end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculates cell center coordinates.
!--------------------------------------------------------------------------------------------------
pure function mesh_cellCenterCoordinates(ip,el)
 
 integer, intent(in) :: el, &                                                                  !< element number
                        ip                                                                     !< integration point number
 real(pReal), dimension(3) :: mesh_cellCenterCoordinates                                             !< x,y,z coordinates of the cell center of the requested IP cell

end function mesh_cellCenterCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume. Allocates global array 'mesh_ipVolume'
!> @details The IP volume is calculated differently depending on the cell type.
!> 2D cells assume an element depth of one in order to calculate the volume.
!> For the hexahedral cell we subdivide the cell into subvolumes of pyramidal
!> shape with a cell face as basis and the central ip at the tip. This subvolume is
!> calculated as an average of four tetrahedals with three corners on the cell face 
!> and one corner at the central ip.
!--------------------------------------------------------------------------------------------------
subroutine mesh_FEM_build_ipVolumes(dimPlex)
 
  PetscInt           :: dimPlex
  PetscReal          :: vol
  PetscReal,  target :: cent(dimPlex), norm(dimPlex)
  PetscReal, pointer :: pCent(:), pNorm(:)
  PetscInt           :: cellStart, cellEnd, cell
  PetscErrorCode     :: ierr
 
  if (.not. allocated(mesh_ipVolume)) then
    allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems))
    mesh_ipVolume = 0.0_pReal 
  endif
 
  call DMPlexGetHeightStratum(geomMesh,0,cellStart,cellEnd,ierr); CHKERRQ(ierr)
  pCent => cent
  pNorm => norm
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
  
  PetscReal,         target :: v0(dimPlex), cellJ(dimPlex*dimPlex), invcellJ(dimPlex*dimPlex)
  PetscReal,        pointer :: pV0(:), pCellJ(:), pInvcellJ(:)
  PetscReal                 :: detJ
  PetscInt                  :: cellStart, cellEnd, cell, qPt, dirI, dirJ, qOffset
  PetscErrorCode            :: ierr
 
 
  allocate(mesh_ipCoordinates(3,mesh_maxNips,mesh_NcpElems),source=0.0_pReal)
 
  pV0 => v0
  pCellJ => cellJ
  pInvcellJ => invcellJ
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

end module mesh
