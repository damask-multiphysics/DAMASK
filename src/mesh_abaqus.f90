!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc, Abaqus and the spectral solver
!--------------------------------------------------------------------------------------------------
module mesh
 use, intrinsic :: iso_c_binding
 use prec, only: pReal, pInt
 use mesh_base

 implicit none
 private
 integer(pInt), public, protected :: &
   mesh_NcpElems, &                                                                                 !< total number of CP elements in local mesh
   mesh_elemType, &                                                                                 !< Element type of the mesh (only support homogeneous meshes)
   mesh_Nnodes, &                                                                                   !< total number of nodes in mesh
   mesh_Ncellnodes, &                                                                               !< total number of cell nodes in mesh (including duplicates)
   mesh_Ncells, &                                                                                   !< total number of cells in mesh
   mesh_NipsPerElem, &                                                                              !< number of IPs in per element
   mesh_NcellnodesPerElem, &                                                                        !< number of cell nodes per  element
   mesh_maxNipNeighbors, &                                                                          !< max number of IP neighbors in any CP element
   mesh_maxNsharedElems                                                                             !< max number of CP elements sharing a node
!!!! BEGIN DEPRECATED !!!!!
 integer(pInt), public, protected :: &
   mesh_maxNips, &                                                                                  !< max number of IPs in any CP element
   mesh_maxNcellnodes                                                                               !< max number of cell nodes in any CP element
!!!! BEGIN DEPRECATED !!!!!

 integer(pInt), dimension(:), allocatable, public, protected :: &
   mesh_homogenizationAt, &                                                                         !< homogenization ID of each element
   mesh_microstructureAt                                                                            !< microstructure ID of each element

 integer(pInt), dimension(:,:), allocatable, public, protected :: &
   mesh_CPnodeID, &                                                                                 !< nodes forming an element
   mesh_element, & !DEPRECATED
   mesh_sharedElem, &                                                                               !< entryCount and list of elements containing node
   mesh_nodeTwins                                                                                   !< node twins are surface nodes that lie exactly on opposite sides of the mesh (surfaces nodes with equal coordinate values in two dimensions)

 integer(pInt), dimension(:,:,:,:), allocatable, public, protected :: &
   mesh_ipNeighborhood                                                                              !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]

 real(pReal), public, protected :: &
   mesh_unitlength                                                                                  !< physical length of one unit in mesh

 real(pReal), dimension(:,:), allocatable, public :: &
   mesh_node, &                                                                                     !< node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)
   mesh_cellnode                                                                                    !< cell node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)

 real(pReal), dimension(:,:), allocatable, public, protected :: &
   mesh_ipVolume, &                                                                                 !< volume associated with IP (initially!)
   mesh_node0                                                                                       !< node x,y,z coordinates (initially!)

 real(pReal), dimension(:,:,:), allocatable, public, protected :: &
   mesh_ipArea                                                                                      !< area of interface to neighboring IP (initially!)

 real(pReal), dimension(:,:,:), allocatable, public :: &
   mesh_ipCoordinates                                                                               !< IP x,y,z coordinates (after deformation!)

 real(pReal),dimension(:,:,:,:), allocatable, public, protected :: &
   mesh_ipAreaNormal                                                                                !< area normal of interface to neighboring IP (initially!)

 logical, dimension(3), public, protected :: mesh_periodicSurface                                   !< flag indicating periodic outer surfaces (used for fluxes)

 integer(pInt), private :: &
   mesh_maxNelemInSet, &
   mesh_Nmaterials

 integer(pInt), dimension(2), private :: &
   mesh_maxValStateVar = 0_pInt

integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_cellnodeParent                                                                              !< cellnode's parent element ID, cellnode's intra-element ID

 integer(pInt),dimension(:,:,:), allocatable, private :: &
   mesh_cell                                                                                        !< cell connectivity for each element,ip/cell

 integer(pInt), dimension(:,:,:), allocatable, private :: &
   FE_nodesAtIP, &                                                                                  !< map IP index to node indices in a specific type of element
   FE_ipNeighbor, &                                                                                 !< +x,-x,+y,-y,+z,-z list of intra-element IPs and(negative) neighbor faces per own IP in a specific type of element
   FE_cell, &                                                                                       !< list of intra-element cell node IDs that constitute the cells in a specific type of element geometry
   FE_cellface                                                                                      !< list of intra-cell cell node IDs that constitute the cell faces of a specific type of cell

 real(pReal), dimension(:,:,:), allocatable, private :: &
   FE_cellnodeParentnodeWeights                                                                     !< list of node weights for the generation of cell nodes

 integer(pInt), dimension(:,:,:,:), allocatable, private :: &
   FE_subNodeOnIPFace

! These definitions should actually reside in the FE-solver specific part (different for MARC/ABAQUS)
! Hence, I suggest to prefix with "FE_"

 integer(pInt), parameter, public :: &
   FE_Nelemtypes = 13_pInt, &
   FE_Ngeomtypes = 10_pInt, &
   FE_Ncelltypes = 4_pInt, &
   FE_maxNnodes = 20_pInt, &
   FE_maxNips = 27_pInt, &
   FE_maxNipNeighbors = 6_pInt, &
   FE_maxmaxNnodesAtIP = 8_pInt, &                                                                  !< max number of (equivalent) nodes attached to an IP
   FE_maxNmatchingNodesPerFace = 4_pInt, &
   FE_maxNfaces = 6_pInt, &
   FE_maxNcellnodes = 64_pInt, &
   FE_maxNcellnodesPerCell = 8_pInt, &
   FE_maxNcellfaces = 6_pInt, &
   FE_maxNcellnodesPerCellface = 4_pInt

 integer(pInt), dimension(FE_Nelemtypes), parameter, public :: FE_geomtype = &                      !< geometry type of particular element type
 int([ &
      1, & ! element   6 (2D 3node 1ip)
      2, & ! element 125 (2D 6node 3ip)
      3, & ! element  11 (2D 4node 4ip)
      4, & ! element  27 (2D 8node 9ip)
      3, & ! element  54 (2D 8node 4ip)
      5, & ! element 134 (3D 4node 1ip)
      6, & ! element 157 (3D 5node 4ip)
      6, & ! element 127 (3D 10node 4ip)
      7, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      9, & ! element   7 (3D 8node 8ip)
      9, & ! element  57 (3D 20node 8ip)
     10  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public  :: FE_celltype = &                     !< cell type that is used by each geometry type
 int([ &
      1, & ! element   6 (2D 3node 1ip)
      2, & ! element 125 (2D 6node 3ip)
      2, & ! element  11 (2D 4node 4ip)
      2, & ! element  27 (2D 8node 9ip)
      3, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      4, & ! element 136 (3D 6node 6ip)
      4, & ! element 117 (3D 8node 1ip)
      4, & ! element   7 (3D 8node 8ip)
      4  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_dimension = &                     !< dimension of geometry type
 int([ &
      2, & ! element   6 (2D 3node 1ip)
      2, & ! element 125 (2D 6node 3ip)
      2, & ! element  11 (2D 4node 4ip)
      2, & ! element  27 (2D 8node 9ip)
      3, & ! element 134 (3D 4node 1ip)
      3, & ! element 127 (3D 10node 4ip)
      3, & ! element 136 (3D 6node 6ip)
      3, & ! element 117 (3D 8node 1ip)
      3, & ! element   7 (3D 8node 8ip)
      3  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Nelemtypes), parameter, public :: FE_Nnodes = &                        !< number of nodes that constitute a specific type of element
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      6, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      8, & ! element  27 (2D 8node 9ip)
      8, & ! element  54 (2D 8node 4ip)
      4, & ! element 134 (3D 4node 1ip)
      5, & ! element 157 (3D 5node 4ip)
     10, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      8, & ! element   7 (3D 8node 8ip)
     20, & ! element  57 (3D 20node 8ip)
     20  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_Nfaces = &                        !< number of faces of a specific type of element geometry
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      3, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      4, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      5, & ! element 136 (3D 6node 6ip)
      6, & ! element 117 (3D 8node 1ip)
      6, & ! element   7 (3D 8node 8ip)
      6  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, private :: FE_NmatchingNodes = &               !< number of nodes that are needed for face matching in a specific type of element geometry
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      3, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      4, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      8, & ! element   7 (3D 8node 8ip)
      8  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_maxNfaces,FE_Ngeomtypes), parameter, private :: &
                                                                       FE_NmatchingNodesPerFace = & !< number of matching nodes per face in a specific type of element geometry
 reshape(int([ &
  2,2,2,0,0,0, & ! element   6 (2D 3node 1ip)
  2,2,2,0,0,0, & ! element 125 (2D 6node 3ip)
  2,2,2,2,0,0, & ! element  11 (2D 4node 4ip)
  2,2,2,2,0,0, & ! element  27 (2D 8node 9ip)
  3,3,3,3,0,0, & ! element 134 (3D 4node 1ip)
  3,3,3,3,0,0, & ! element 127 (3D 10node 4ip)
  3,4,4,4,3,0, & ! element 136 (3D 6node 6ip)
  4,4,4,4,4,4, & ! element 117 (3D 8node 1ip)
  4,4,4,4,4,4, & ! element   7 (3D 8node 8ip)
  4,4,4,4,4,4  & ! element  21 (3D 20node 27ip)
  ],pInt),[FE_maxNipNeighbors,FE_Ngeomtypes])

 integer(pInt), dimension(FE_maxNmatchingNodesPerFace,FE_maxNfaces,FE_Ngeomtypes), &
                                                          parameter, private :: FE_face = &         !< List of node indices on each face of a specific type of element geometry
 reshape(int([&
  1,2,0,0 , & ! element   6 (2D 3node 1ip)
  2,3,0,0 , &
  3,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element 125 (2D 6node 3ip)
  2,3,0,0 , &
  3,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element  11 (2D 4node 4ip)
  2,3,0,0 , &
  3,4,0,0 , &
  4,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element  27 (2D 8node 9ip)
  2,3,0,0 , &
  3,4,0,0 , &
  4,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 134 (3D 4node 1ip)
  1,4,2,0 , &
  2,3,4,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 127 (3D 10node 4ip)
  1,4,2,0 , &
  2,4,3,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 136 (3D 6node 6ip)
  1,4,5,2 , &
  2,5,6,3 , &
  1,3,6,4 , &
  4,6,5,0 , &
  0,0,0,0 , &
  1,2,3,4 , & ! element 117 (3D 8node 1ip)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,4 , & ! element   7 (3D 8node 8ip)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,4 , & ! element  21 (3D 20node 27ip)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5   &
  ],pInt),[FE_maxNmatchingNodesPerFace,FE_maxNfaces,FE_Ngeomtypes])

 integer(pInt), dimension(FE_Ngeomtypes), parameter, private :: FE_Ncellnodes = &                   !< number of cell nodes in a specific geometry type
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      7, & ! element 125 (2D 6node 3ip)
      9, & ! element  11 (2D 4node 4ip)
     16, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
     15, & ! element 127 (3D 10node 4ip)
     21, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
     27, & ! element   7 (3D 8node 8ip)
     64  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCell = &             !< number of cell nodes in a specific cell type
 int([ &
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      8  & ! (3D 8node)
  ],pInt)

 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCellface = &        !< number of cell nodes per cell face in a specific cell type
 int([&
      2, & ! (2D 3node)
      2, & ! (2D 4node)
      3, & ! (3D 4node)
      4  & ! (3D 8node)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, public :: FE_Nips = &                          !< number of IPs in a specific type of element
 int([ &
      1, & ! element   6 (2D 3node 1ip)
      3, & ! element 125 (2D 6node 3ip)
      4, & ! element  11 (2D 4node 4ip)
      9, & ! element  27 (2D 8node 9ip)
      1, & ! element 134 (3D 4node 1ip)
      4, & ! element 127 (3D 10node 4ip)
      6, & ! element 136 (3D 6node 6ip)
      1, & ! element 117 (3D 8node 1ip)
      8, & ! element   7 (3D 8node 8ip)
     27  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), dimension(FE_Ncelltypes), parameter, public :: FE_NipNeighbors = &                  !< number of ip neighbors / cell faces in a specific cell type
 int([&
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      6  & ! (3D 8node)
  ],pInt)

 integer(pInt), dimension(FE_Ngeomtypes), parameter, private :: FE_maxNnodesAtIP = &                !< maximum number of parent nodes that belong to an IP for a specific type of element
 int([ &
      3, & ! element   6 (2D 3node 1ip)
      1, & ! element 125 (2D 6node 3ip)
      1, & ! element  11 (2D 4node 4ip)
      2, & ! element  27 (2D 8node 9ip)
      4, & ! element 134 (3D 4node 1ip)
      1, & ! element 127 (3D 10node 4ip)
      1, & ! element 136 (3D 6node 6ip)
      8, & ! element 117 (3D 8node 1ip)
      1, & ! element   7 (3D 8node 8ip)
      4  & ! element  21 (3D 20node 27ip)
  ],pInt)

 integer(pInt), private :: &
   mesh_Nelems, &                                                                                   !< total number of elements in mesh (including non-DAMASK elements)
   mesh_maxNnodes, &                                                                                !< max number of nodes in any CP element
   mesh_NelemSets
 character(len=64), dimension(:), allocatable, private :: &
   mesh_nameElemSet, &                                                                              !< names of elementSet
   mesh_nameMaterial, &                                                                             !< names of material in solid section
   mesh_mapMaterial                                                                                 !< name of elementSet for material
 integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
 integer(pInt), dimension(:,:), allocatable, target, private :: &
   mesh_mapFEtoCPelem, &                                                                            !< [sorted FEid, corresponding CPid]
   mesh_mapFEtoCPnode                                                                               !< [sorted FEid, corresponding CPid]
 logical, private :: noPart                                                                         !< for cases where the ABAQUS input file does not use part/assembly information

 public :: &
   mesh_init, &
   mesh_build_cellnodes, &
   mesh_build_ipVolumes, &
   mesh_build_ipCoordinates, &
   mesh_cellCenterCoordinates, &
   mesh_FEasCP

 private :: &
   mesh_get_damaskOptions, &
   mesh_build_cellconnectivity, &
   mesh_build_ipAreas, &
   mesh_tell_statistics, &
   FE_mapElemtype, &
   mesh_faceMatch, &
   mesh_build_FEdata, &
   mesh_build_nodeTwins, &
   mesh_build_sharedElems, &
   mesh_build_ipNeighborhood, &
   mesh_abaqus_count_nodesAndElements, &
   mesh_abaqus_count_elementSets, &
   mesh_abaqus_count_materials, &
   mesh_abaqus_map_elementSets, &
   mesh_abaqus_map_materials, &
   mesh_abaqus_count_cpElements, &
   mesh_abaqus_map_elements, &
   mesh_abaqus_map_nodes, &
   mesh_abaqus_build_nodes, &
   mesh_abaqus_count_cpSizes, &
   mesh_abaqus_build_elements


 type, public, extends(tMesh) :: tMesh_Abaqus
 
 integer(pInt):: &
   mesh_Nelems, &                                                                                   !< total number of elements in mesh (including non-DAMASK elements)
   mesh_maxNnodes, &                                                                                !< max number of nodes in any CP element
   mesh_NelemSets, &
   mesh_maxNelemInSet, &
   mesh_Nmaterials
 character(len=64), dimension(:), allocatable :: &
   mesh_nameElemSet, &                                                                              !< names of elementSet
   mesh_nameMaterial, &                                                                             !< names of material in solid section
   mesh_mapMaterial                                                                                 !< name of elementSet for material
 integer(pInt), dimension(:,:), allocatable :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
 logical:: noPart                                                                         !< for cases where the ABAQUS input file does not use part/assembly information

 contains 
   procedure :: init=>tMesh_abaqus_init
 end type tMesh_Abaqus
 
 type(tMesh_Abaqus), public, protected :: theMesh
 
contains

subroutine tMesh_abaqus_init(self)
 implicit none
 class(tMesh_abaqus) :: self
 
end subroutine tMesh_abaqus_init

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use DAMASK_interface
 use IO, only: &
   IO_open_InputFile, &
   IO_timeStamp, &
   IO_error, &
   IO_write_jobFile
 use debug, only: &
   debug_e, &
   debug_i, &
   debug_level, &
   debug_mesh, &
   debug_levelBasic
 use numerics, only: &
   usePingPong, &
   numerics_unitlength, &
   worldrank
 use FEsolving, only: &
   modelName, &
   calcMode, &   FEsolving_execElem, &
   FEsolving_execIP

 implicit none
 integer(pInt), parameter :: FILEUNIT = 222_pInt
 integer(pInt), intent(in), optional :: el, ip
 integer(pInt) :: j
 logical :: myDebug

 write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 call mesh_build_FEdata                                                                             ! get properties of the different types of elements
 mesh_unitlength = numerics_unitlength                                                              ! set physical extent of a length unit in mesh

 myDebug = (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt)

 call IO_open_inputFile(FILEUNIT,modelName)                                                         ! parse info from input file...
 if (myDebug) write(6,'(a)') ' Opened input file'; flush(6)
 noPart = hasNoPart(FILEUNIT)
 call mesh_abaqus_count_nodesAndElements(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted nodes/elements'; flush(6)
 call mesh_abaqus_count_elementSets(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted element sets'; flush(6)
 call mesh_abaqus_count_materials(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted materials'; flush(6)
 call mesh_abaqus_map_elementSets(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Mapped element sets'; flush(6)
 call mesh_abaqus_map_materials(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Mapped materials'; flush(6)
 call mesh_abaqus_count_cpElements(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted CP elements'; flush(6)
 call mesh_abaqus_map_elements(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Mapped elements'; flush(6)
 call mesh_abaqus_map_nodes(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Mapped nodes'; flush(6)
 call mesh_abaqus_build_nodes(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Built nodes'; flush(6)
 call mesh_abaqus_count_cpSizes(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted CP sizes'; flush(6)
 call mesh_abaqus_build_elements(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Built elements'; flush(6)
 call mesh_get_damaskOptions(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Got DAMASK options'; flush(6)
 call mesh_build_cellconnectivity
 if (myDebug) write(6,'(a)') ' Built cell connectivity'; flush(6)
 mesh_cellnode = mesh_build_cellnodes(mesh_node,mesh_Ncellnodes)
 if (myDebug) write(6,'(a)') ' Built cell nodes'; flush(6)
 call mesh_build_ipCoordinates
 if (myDebug) write(6,'(a)') ' Built IP coordinates'; flush(6)
 call mesh_build_ipVolumes
 if (myDebug) write(6,'(a)') ' Built IP volumes'; flush(6)
 call mesh_build_ipAreas
 if (myDebug) write(6,'(a)') ' Built IP areas'; flush(6)
 close (FILEUNIT)
 call mesh_build_nodeTwins
 if (myDebug) write(6,'(a)') ' Built node twins'; flush(6)
 call mesh_build_sharedElems
 if (myDebug) write(6,'(a)') ' Built shared elements'; flush(6)
 call mesh_build_ipNeighborhood
 if (myDebug) write(6,'(a)') ' Built IP neighborhood'; flush(6)

 if (worldrank == 0_pInt) then
   call mesh_tell_statistics
 endif

 if (usePingPong .and. (mesh_Nelems /= mesh_NcpElems)) &
   call IO_error(600_pInt)                                                                          ! ping-pong must be disabled when having non-DAMASK elements
 if (debug_e < 1 .or. debug_e > mesh_NcpElems) &
   call IO_error(602_pInt,ext_msg='element')                                                        ! selected element does not exist
 if (debug_i < 1 .or. debug_i > FE_Nips(FE_geomtype(mesh_element(2_pInt,debug_e)))) &
   call IO_error(602_pInt,ext_msg='IP')                                                             ! selected element does not have requested IP
 FEsolving_execElem = [ 1_pInt,mesh_NcpElems ]                                                      ! parallel loop bounds set to comprise all DAMASK elements
 allocate(FEsolving_execIP(2_pInt,mesh_NcpElems), source=1_pInt)                                    ! parallel loop bounds set to comprise from first IP...
 forall (j = 1_pInt:mesh_NcpElems) FEsolving_execIP(2,j) = FE_Nips(FE_geomtype(mesh_element(2,j)))  ! ...up to own IP count for each element
 allocate(calcMode(mesh_maxNips,mesh_NcpElems))
 calcMode = .false.                                                                                 ! pretend to have collected what first call is asking (F = I)
 calcMode(ip,mesh_FEasCP('elem',el)) = .true.                                                       ! first ip,el needs to be already pingponged to "calc"

!!!! COMPATIBILITY HACK !!!!
! for a homogeneous mesh, all elements have the same number of IPs and and cell nodes.
! hence, xxPerElem instead of maxXX
 mesh_NipsPerElem       = mesh_maxNips
 mesh_NcellnodesPerElem = mesh_maxNcellnodes
! better name
 mesh_homogenizationAt  = mesh_element(3,:)
 mesh_microstructureAt  = mesh_element(4,:)
 mesh_CPnodeID          = mesh_element(5:4+mesh_NipsPerElem,:)
!!!!!!!!!!!!!!!!!!!!!!!!

contains
!--------------------------------------------------------------------------------------------------
!> @brief check if the input file for Abaqus contains part info
!--------------------------------------------------------------------------------------------------
logical function hasNoPart(fileUnit)
 use IO, only: &
   IO_stringPos, &
   IO_stringValue, &
   IO_lc
 
 implicit none
 integer(pInt),    intent(in)                :: fileUnit

 integer(pInt), allocatable, dimension(:)    :: chunkPos
 character(len=65536)                        :: line

 hasNoPart = .true.

610 FORMAT(A65536)
 rewind(fileUnit)
 do
   read(fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if (IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) then
     hasNoPart = .false.
     exit
   endif
 enddo

620 end function hasNoPart

end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Gives the FE to CP ID mapping by binary search through lookup array
!! valid questions (what) are 'elem', 'node'
!--------------------------------------------------------------------------------------------------
integer(pInt) function mesh_FEasCP(what,myID)
 use IO, only: &
   IO_lc

 implicit none
 character(len=*), intent(in) :: what
 integer(pInt),    intent(in) :: myID

 integer(pInt), dimension(:,:), pointer :: lookupMap
 integer(pInt) :: lower,upper,center

 mesh_FEasCP = 0_pInt
 select case(IO_lc(what(1:4)))
   case('elem')
     lookupMap => mesh_mapFEtoCPelem
   case('node')
     lookupMap => mesh_mapFEtoCPnode
   case default
     return
 endselect

 lower = 1_pInt
 upper = int(size(lookupMap,2_pInt),pInt)

 if (lookupMap(1_pInt,lower) == myID) then                                                          ! check at bounds QUESTION is it valid to extend bounds by 1 and just do binary search w/o init check at bounds?
   mesh_FEasCP = lookupMap(2_pInt,lower)
   return
 elseif (lookupMap(1_pInt,upper) == myID) then
   mesh_FEasCP = lookupMap(2_pInt,upper)
   return
 endif
 binarySearch: do while (upper-lower > 1_pInt)
   center = (lower+upper)/2_pInt
   if (lookupMap(1_pInt,center) < myID) then
     lower = center
   elseif (lookupMap(1_pInt,center) > myID) then
     upper = center
   else
     mesh_FEasCP = lookupMap(2_pInt,center)
     exit
   endif
 enddo binarySearch

end function mesh_FEasCP


!--------------------------------------------------------------------------------------------------
!> @brief Split CP elements into cells.
!> @details Build a mapping between cells and the corresponding cell nodes ('mesh_cell').
!> Cell nodes that are also matching nodes are unique in the list of cell nodes,
!> all others (currently) might be stored more than once.
!> Also allocates the 'mesh_node' array.
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_cellconnectivity

 implicit none
 integer(pInt), dimension(:), allocatable :: &
   matchingNode2cellnode
 integer(pInt), dimension(:,:), allocatable :: &
   cellnodeParent
 integer(pInt), dimension(mesh_maxNcellnodes) :: &
   localCellnode2globalCellnode
 integer(pInt) :: &
   e,t,g,c,n,i, &
   matchingNodeID, &
   localCellnodeID

 allocate(mesh_cell(FE_maxNcellnodesPerCell,mesh_maxNips,mesh_NcpElems), source=0_pInt)
 allocate(matchingNode2cellnode(mesh_Nnodes),                            source=0_pInt)
 allocate(cellnodeParent(2_pInt,mesh_maxNcellnodes*mesh_NcpElems),       source=0_pInt)

!--------------------------------------------------------------------------------------------------
! Count cell nodes (including duplicates) and generate cell connectivity list
 mesh_Ncellnodes = 0_pInt
 mesh_Ncells = 0_pInt
 do e = 1_pInt,mesh_NcpElems                                                                        ! loop over cpElems
   t = mesh_element(2_pInt,e)                                                                       ! get element type
   g = FE_geomtype(t)                                                                               ! get geometry type
   c = FE_celltype(g)                                                                               ! get cell type
   localCellnode2globalCellnode = 0_pInt
   mesh_Ncells = mesh_Ncells + FE_Nips(g)
   do i = 1_pInt,FE_Nips(g)                                                                         ! loop over ips=cells in this element
     do n = 1_pInt,FE_NcellnodesPerCell(c)                                                          ! loop over cell nodes in this cell
       localCellnodeID = FE_cell(n,i,g)
       if (localCellnodeID <= FE_NmatchingNodes(g)) then                                            ! this cell node is a matching node
         matchingNodeID = mesh_element(4_pInt+localCellnodeID,e)
         if (matchingNode2cellnode(matchingNodeID) == 0_pInt) then                                  ! if this matching node does not yet exist in the glbal cell node list ...
           mesh_Ncellnodes = mesh_Ncellnodes + 1_pInt                                               ! ... count it as cell node ...
           matchingNode2cellnode(matchingNodeID) = mesh_Ncellnodes                                  ! ... and remember its global ID
           cellnodeParent(1_pInt,mesh_Ncellnodes) = e                                               ! ... and where it belongs to
           cellnodeParent(2_pInt,mesh_Ncellnodes) = localCellnodeID
         endif
         mesh_cell(n,i,e) = matchingNode2cellnode(matchingNodeID)
       else                                                                                         ! this cell node is no matching node
         if (localCellnode2globalCellnode(localCellnodeID) == 0_pInt) then                          ! if this local cell node does not yet exist in the  global cell node list ...
           mesh_Ncellnodes = mesh_Ncellnodes + 1_pInt                                               ! ... count it as cell node ...
           localCellnode2globalCellnode(localCellnodeID) = mesh_Ncellnodes                          ! ... and remember its global ID ...
           cellnodeParent(1_pInt,mesh_Ncellnodes) = e                                               ! ... and it belongs to
           cellnodeParent(2_pInt,mesh_Ncellnodes) = localCellnodeID
         endif
         mesh_cell(n,i,e) = localCellnode2globalCellnode(localCellnodeID)
       endif
     enddo
   enddo
 enddo

 allocate(mesh_cellnodeParent(2_pInt,mesh_Ncellnodes))
 allocate(mesh_cellnode(3_pInt,mesh_Ncellnodes))
 forall(n = 1_pInt:mesh_Ncellnodes)
   mesh_cellnodeParent(1,n) = cellnodeParent(1,n)
   mesh_cellnodeParent(2,n) = cellnodeParent(2,n)
 endforall

end subroutine mesh_build_cellconnectivity


!--------------------------------------------------------------------------------------------------
!> @brief Calculate position of cellnodes from the given position of nodes
!> Build list of cellnodes' coordinates.
!> Cellnode coordinates are calculated from a weighted sum of node coordinates.
!--------------------------------------------------------------------------------------------------
function mesh_build_cellnodes(nodes,Ncellnodes)

 implicit none
 integer(pInt),                         intent(in) :: Ncellnodes                                    !< requested number of cellnodes
 real(pReal), dimension(3,mesh_Nnodes), intent(in) :: nodes
 real(pReal), dimension(3,Ncellnodes) :: mesh_build_cellnodes

 integer(pInt) :: &
   e,t,n,m, &
   localCellnodeID
 real(pReal), dimension(3) :: &
   myCoords

 mesh_build_cellnodes = 0.0_pReal
!$OMP PARALLEL DO PRIVATE(e,localCellnodeID,t,myCoords)
 do n = 1_pInt,Ncellnodes                                                                           ! loop over cell nodes
   e = mesh_cellnodeParent(1,n)
   localCellnodeID = mesh_cellnodeParent(2,n)
   t = mesh_element(2,e)                                                                            ! get element type
   myCoords = 0.0_pReal
   do m = 1_pInt,FE_Nnodes(t)
     myCoords = myCoords + nodes(1:3,mesh_element(4_pInt+m,e)) &
                         * FE_cellnodeParentnodeWeights(m,localCellnodeID,t)
   enddo
   mesh_build_cellnodes(1:3,n) = myCoords / sum(FE_cellnodeParentnodeWeights(:,localCellnodeID,t))
 enddo
!$OMP END PARALLEL DO

end function mesh_build_cellnodes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume. Allocates global array 'mesh_ipVolume'
!> @details The IP volume is calculated differently depending on the cell type.
!> 2D cells assume an element depth of one in order to calculate the volume.
!> For the hexahedral cell we subdivide the cell into subvolumes of pyramidal
!> shape with a cell face as basis and the central ip at the tip. This subvolume is
!> calculated as an average of four tetrahedals with three corners on the cell face
!> and one corner at the central ip.
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipVolumes
 use math, only: &
   math_volTetrahedron, &
   math_areaTriangle

 implicit none
 integer(pInt) ::                                e,t,g,c,i,m,f,n
 real(pReal), dimension(FE_maxNcellnodesPerCellface,FE_maxNcellfaces) :: subvolume

 if (.not. allocated(mesh_ipVolume)) then
   allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems))
   mesh_ipVolume = 0.0_pReal
 endif

 !$OMP PARALLEL DO PRIVATE(t,g,c,m,subvolume)
   do e = 1_pInt,mesh_NcpElems                                                                      ! loop over cpElems
     t = mesh_element(2_pInt,e)                                                                     ! get element type
     g = FE_geomtype(t)                                                                             ! get geometry type
     c = FE_celltype(g)                                                                             ! get cell type
     select case (c)

       case (1_pInt)                                                                                ! 2D 3node
         forall (i = 1_pInt:FE_Nips(g)) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_areaTriangle(mesh_cellnode(1:3,mesh_cell(1,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(3,i,e)))

       case (2_pInt)                                                                                ! 2D 4node
         forall (i = 1_pInt:FE_Nips(g)) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_areaTriangle(mesh_cellnode(1:3,mesh_cell(1,i,e)), &            ! here we assume a planar shape, so division in two triangles suffices
                                                  mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(3,i,e))) &
                              + math_areaTriangle(mesh_cellnode(1:3,mesh_cell(3,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(4,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(1,i,e)))

       case (3_pInt)                                                                                ! 3D 4node
         forall (i = 1_pInt:FE_Nips(g)) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_volTetrahedron(mesh_cellnode(1:3,mesh_cell(1,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(3,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(4,i,e)))

       case (4_pInt)                                                                                ! 3D 8node
         m = FE_NcellnodesPerCellface(c)
         do i = 1_pInt,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           subvolume = 0.0_pReal
           forall(f = 1_pInt:FE_NipNeighbors(c), n = 1_pInt:FE_NcellnodesPerCellface(c)) &
             subvolume(n,f) = math_volTetrahedron(&
                                mesh_cellnode(1:3,mesh_cell(FE_cellface(      n     ,f,c),i,e)), &
                                mesh_cellnode(1:3,mesh_cell(FE_cellface(1+mod(n  ,m),f,c),i,e)), &
                                mesh_cellnode(1:3,mesh_cell(FE_cellface(1+mod(n+1,m),f,c),i,e)), &
                                mesh_ipCoordinates(1:3,i,e))
           mesh_ipVolume(i,e) = 0.5_pReal * sum(subvolume)                                         ! each subvolume is based on four tetrahedrons, altough the face consists of only two triangles -> averaging factor two
         enddo

     end select
   enddo
 !$OMP END PARALLEL DO

end subroutine mesh_build_ipVolumes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP Coordinates. Allocates global array 'mesh_ipCoordinates'
! Called by all solvers in mesh_init in order to initialize the ip coordinates.
! Later on the current ip coordinates are directly prvided by the spectral solver and by Abaqus,
! so no need to use this subroutine anymore; Marc however only provides nodal displacements,
! so in this case the ip coordinates are always calculated on the basis of this subroutine.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR THE MOMENT THIS SUBROUTINE ACTUALLY CALCULATES THE CELL CENTER AND NOT THE IP COORDINATES,
! AS THE IP IS NOT (ALWAYS) LOCATED IN THE CENTER OF THE IP VOLUME.
! HAS TO BE CHANGED IN A LATER VERSION.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipCoordinates

 implicit none
 integer(pInt) :: e,t,g,c,i,n
 real(pReal), dimension(3) :: myCoords

 if (.not. allocated(mesh_ipCoordinates)) &
   allocate(mesh_ipCoordinates(3,mesh_maxNips,mesh_NcpElems),source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,myCoords)
 do e = 1_pInt,mesh_NcpElems                                                                        ! loop over cpElems
   t = mesh_element(2_pInt,e)                                                                       ! get element type
   g = FE_geomtype(t)                                                                               ! get geometry type
   c = FE_celltype(g)                                                                               ! get cell type
   do i = 1_pInt,FE_Nips(g)                                                                         ! loop over ips=cells in this element
     myCoords = 0.0_pReal
     do n = 1_pInt,FE_NcellnodesPerCell(c)                                                          ! loop over cell nodes in this cell
       myCoords = myCoords + mesh_cellnode(1:3,mesh_cell(n,i,e))
     enddo
     mesh_ipCoordinates(1:3,i,e) = myCoords / real(FE_NcellnodesPerCell(c),pReal)
   enddo
 enddo
 !$OMP END PARALLEL DO

end subroutine mesh_build_ipCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief Calculates cell center coordinates.
!--------------------------------------------------------------------------------------------------
pure function mesh_cellCenterCoordinates(ip,el)

 implicit none
 integer(pInt), intent(in) :: el, &                                                                  !< element number
                              ip                                                                     !< integration point number
 real(pReal), dimension(3) :: mesh_cellCenterCoordinates                                             !< x,y,z coordinates of the cell center of the requested IP cell
 integer(pInt) :: t,g,c,n

 t = mesh_element(2_pInt,el)                                                                         ! get element type
 g = FE_geomtype(t)                                                                                  ! get geometry type
 c = FE_celltype(g)                                                                                  ! get cell type
 mesh_cellCenterCoordinates = 0.0_pReal
 do n = 1_pInt,FE_NcellnodesPerCell(c)                                                               ! loop over cell nodes in this cell
   mesh_cellCenterCoordinates = mesh_cellCenterCoordinates + mesh_cellnode(1:3,mesh_cell(n,ip,el))
 enddo
 mesh_cellCenterCoordinates = mesh_cellCenterCoordinates / real(FE_NcellnodesPerCell(c),pReal)

 end function mesh_cellCenterCoordinates



!--------------------------------------------------------------------------------------------------
!> @brief builds mesh of (distorted) cubes for given coordinates (= center of the cubes)
!--------------------------------------------------------------------------------------------------
function mesh_nodesAroundCentres(gDim,Favg,centres) result(nodes)
 use debug, only: &
   debug_mesh, &
   debug_level, &
   debug_levelBasic
 use math, only: &
   math_mul33x3

 implicit none
 real(pReal), intent(in), dimension(:,:,:,:) :: &
   centres
 real(pReal),             dimension(3,size(centres,2)+1,size(centres,3)+1,size(centres,4)+1) :: &
   nodes
 real(pReal), intent(in), dimension(3) :: &
   gDim
 real(pReal), intent(in), dimension(3,3) :: &
   Favg
 real(pReal),             dimension(3,size(centres,2)+2,size(centres,3)+2,size(centres,4)+2) :: &
   wrappedCentres

 integer(pInt) :: &
   i,j,k,n
 integer(pInt),           dimension(3), parameter :: &
   diag = 1_pInt
 integer(pInt),           dimension(3) :: &
   shift = 0_pInt, &
   lookup = 0_pInt, &
   me = 0_pInt, &
   iRes = 0_pInt
 integer(pInt),           dimension(3,8) :: &
   neighbor = reshape([ &
                       0_pInt, 0_pInt, 0_pInt, &
                       1_pInt, 0_pInt, 0_pInt, &
                       1_pInt, 1_pInt, 0_pInt, &
                       0_pInt, 1_pInt, 0_pInt, &
                       0_pInt, 0_pInt, 1_pInt, &
                       1_pInt, 0_pInt, 1_pInt, &
                       1_pInt, 1_pInt, 1_pInt, &
                       0_pInt, 1_pInt, 1_pInt  ], [3,8])

!--------------------------------------------------------------------------------------------------
! initializing variables
 iRes =  [size(centres,2),size(centres,3),size(centres,4)]
 nodes = 0.0_pReal
 wrappedCentres = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! report
 if (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt) then
   write(6,'(a)')          ' Meshing cubes around centroids'
   write(6,'(a,3(e12.5))') ' Dimension: ', gDim
   write(6,'(a,3(i5))')    ' Resolution:', iRes
 endif

!--------------------------------------------------------------------------------------------------
! building wrappedCentres = centroids + ghosts
 wrappedCentres(1:3,2_pInt:iRes(1)+1_pInt,2_pInt:iRes(2)+1_pInt,2_pInt:iRes(3)+1_pInt) = centres
 do k = 0_pInt,iRes(3)+1_pInt
   do j = 0_pInt,iRes(2)+1_pInt
     do i = 0_pInt,iRes(1)+1_pInt
       if (k==0_pInt .or. k==iRes(3)+1_pInt .or. &                                                  ! z skin
           j==0_pInt .or. j==iRes(2)+1_pInt .or. &                                                  ! y skin
           i==0_pInt .or. i==iRes(1)+1_pInt      ) then                                             ! x skin
         me = [i,j,k]                                                                               ! me on skin
         shift = sign(abs(iRes+diag-2_pInt*me)/(iRes+diag),iRes+diag-2_pInt*me)
         lookup = me-diag+shift*iRes
         wrappedCentres(1:3,i+1_pInt,        j+1_pInt,        k+1_pInt) = &
                centres(1:3,lookup(1)+1_pInt,lookup(2)+1_pInt,lookup(3)+1_pInt) &
                - math_mul33x3(Favg, real(shift,pReal)*gDim)
       endif
 enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! averaging
 do k = 0_pInt,iRes(3); do j = 0_pInt,iRes(2); do i = 0_pInt,iRes(1)
   do n = 1_pInt,8_pInt
    nodes(1:3,i+1_pInt,j+1_pInt,k+1_pInt) = &
    nodes(1:3,i+1_pInt,j+1_pInt,k+1_pInt) + wrappedCentres(1:3,i+1_pInt+neighbor(1,n), &
                                                               j+1_pInt+neighbor(2,n), &
                                                               k+1_pInt+neighbor(3,n) )
   enddo
 enddo; enddo; enddo
 nodes = nodes/8.0_pReal

end function mesh_nodesAroundCentres


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of nodes and elements in mesh and stores them in
!! 'mesh_Nelems' and 'mesh_Nnodes'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_nodesAndElements(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 logical :: inPart

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if (inPart .or. noPart) then
     select case ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)))
       case('*node')
          if( &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'print'    .and. &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'file'     .and. &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'response' &
             ) &
            mesh_Nnodes = mesh_Nnodes + IO_countDataLines(fileUnit)
       case('*element')
          if( &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'matrix'   .and. &
              IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'response' &
             ) then
            mesh_Nelems = mesh_Nelems + IO_countDataLines(fileUnit)
          endif
     endselect
   endif
 enddo

620 if (mesh_Nnodes < 2_pInt)  call IO_error(error_ID=900_pInt)
 if (mesh_Nelems == 0_pInt) call IO_error(error_ID=901_pInt)

end subroutine mesh_abaqus_count_nodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief count overall number of element sets in mesh and write 'mesh_NelemSets' and
!! 'mesh_maxNelemInSet'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_elementSets(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 logical :: inPart

 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = mesh_Nelems                                                                   ! have to be conservative, since Abaqus allows for recursive definitons

610 FORMAT(A300)

 inPart = .false.
 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*elset' ) &
     mesh_NelemSets = mesh_NelemSets + 1_pInt
 enddo

620 continue
 if (mesh_NelemSets == 0) call IO_error(error_ID=902_pInt)

end subroutine mesh_abaqus_count_elementSets


!--------------------------------------------------------------------------------------------------
! count overall number of solid sections sets in mesh (Abaqus only)
!
! mesh_Nmaterials
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_materials(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 logical inPart

 mesh_Nmaterials = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,chunkPos,1_pInt)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,chunkPos,2_pInt)) == 'section' ) &
     mesh_Nmaterials = mesh_Nmaterials + 1_pInt
 enddo

620 if (mesh_Nmaterials == 0_pInt) call IO_error(error_ID=903_pInt)

end subroutine mesh_abaqus_count_materials


!--------------------------------------------------------------------------------------------------
! Build element set mapping
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_elementSets(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_continuousIntValues, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) :: elemSet = 0_pInt,i
 logical :: inPart = .false.

 allocate (mesh_nameElemSet(mesh_NelemSets)); mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1_pInt+mesh_maxNelemInSet,mesh_NelemSets),source=0_pInt)

610 FORMAT(A300)


 rewind(fileUnit)
 do
   read (fileUnit,610,END=640) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*elset' ) then
     elemSet = elemSet + 1_pInt
     mesh_nameElemSet(elemSet)  = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2_pInt)),'elset'))
     mesh_mapElemSet(:,elemSet) = IO_continuousIntValues(fileUnit,mesh_Nelems,mesh_nameElemSet,&
                                          mesh_mapElemSet,elemSet-1_pInt)
   endif
 enddo

640 do i = 1_pInt,elemSet
   if (mesh_mapElemSet(1,i) == 0_pInt) call IO_error(error_ID=904_pInt,ext_msg=mesh_nameElemSet(i))
 enddo

end subroutine mesh_abaqus_map_elementSets


!--------------------------------------------------------------------------------------------------
! map solid section (Abaqus only)
!
! allocate globals: mesh_nameMaterial, mesh_mapMaterial
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_materials(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line

 integer(pInt) :: i,c = 0_pInt
 logical :: inPart = .false.
 character(len=64) :: elemSetName,materialName

 allocate (mesh_nameMaterial(mesh_Nmaterials)); mesh_nameMaterial = ''
 allocate (mesh_mapMaterial(mesh_Nmaterials));  mesh_mapMaterial = ''

610 FORMAT(A300)

 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,chunkPos,1_pInt)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,chunkPos,2_pInt)) == 'section' ) then

     elemSetName = ''
     materialName = ''

     do i = 3_pInt,chunkPos(1_pInt)
       if (IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'elset') /= '') &
         elemSetName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'elset'))
       if (IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'material') /= '') &
         materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'material'))
     enddo

     if (elemSetName /= '' .and. materialName /= '') then
       c = c + 1_pInt
       mesh_nameMaterial(c) = materialName                                                          ! name of material used for this section
       mesh_mapMaterial(c)  = elemSetName                                                           ! mapped to respective element set
     endif
   endif
 enddo

620 if (c==0_pInt) call IO_error(error_ID=905_pInt)
 do i=1_pInt,c
   if (mesh_nameMaterial(i)=='' .or. mesh_mapMaterial(i)=='') call IO_error(error_ID=905_pInt)
 enddo

 end subroutine mesh_abaqus_map_materials


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of CP elements in mesh and stores them in 'mesh_NcpElems'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_cpElements(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_extractValue

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line
 integer(pInt) :: i,k
 logical :: materialFound = .false.
 character(len=64) ::materialName,elemSetName

 mesh_NcpElems = 0_pInt

610 FORMAT(A300)

 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   select case ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if (IO_lc(IO_StringValue(line,chunkPos,2_pInt)) == 'material' .and. materialFound) then
         do i = 1_pInt,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) &                                            ! matched?
                 mesh_NcpElems = mesh_NcpElems + mesh_mapElemSet(1,k)                               ! add those elem count
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo

620 if (mesh_NcpElems == 0_pInt) call IO_error(error_ID=906_pInt)

end subroutine mesh_abaqus_count_cpElements


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPelem'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_elements(fileUnit)

 use math, only: math_qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) ::i,j,k,cpElem = 0_pInt
 logical :: materialFound = .false.
 character (len=64) materialName,elemSetName                                                        ! why limited to 64? ABAQUS?

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems), source = 0_pInt)

610 FORMAT(A300)

 rewind(fileUnit)
 do
   read (fileUnit,610,END=660) line
   chunkPos = IO_stringPos(line)
   select case ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if (IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'material' .and. materialFound) then
         do i = 1_pInt,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) then                                         ! matched?
                 do j = 1_pInt,mesh_mapElemSet(1,k)
                   cpElem = cpElem + 1_pInt
                   mesh_mapFEtoCPelem(1,cpElem) = mesh_mapElemSet(1_pInt+j,k)                       ! store FE id
                   mesh_mapFEtoCPelem(2,cpElem) = cpElem                                            ! store our id
                 enddo
               endif
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo

660 call math_qsort(mesh_mapFEtoCPelem,1_pInt,int(size(mesh_mapFEtoCPelem,2_pInt),pInt))            ! should be mesh_NcpElems

 if (int(size(mesh_mapFEtoCPelem),pInt) < 2_pInt) call IO_error(error_ID=907_pInt)

end subroutine mesh_abaqus_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPnode'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_nodes(fileUnit)

 use math, only: math_qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_intValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line

 integer(pInt) :: i,c,cpNode = 0_pInt
 logical :: inPart = .false.

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes), source=0_pInt)

610 FORMAT(A300)

 rewind(fileUnit)
 do
   read (fileUnit,610,END=650) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'response' ) &
   ) then
     c = IO_countDataLines(fileUnit)
     do i = 1_pInt,c
       backspace(fileUnit)
     enddo
     do i = 1_pInt,c
       read (fileUnit,610,END=650) line
       chunkPos = IO_stringPos(line)
       cpNode = cpNode + 1_pInt
       mesh_mapFEtoCPnode(1_pInt,cpNode) = IO_intValue(line,chunkPos,1_pInt)
       mesh_mapFEtoCPnode(2_pInt,cpNode) = cpNode
     enddo
   endif
 enddo

650 call math_qsort(mesh_mapFEtoCPnode,1_pInt,int(size(mesh_mapFEtoCPnode,2_pInt),pInt))

 if (int(size(mesh_mapFEtoCPnode),pInt) == 0_pInt) call IO_error(error_ID=908_pInt)

end subroutine mesh_abaqus_map_nodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!! Allocates global arrays 'mesh_node0' and 'mesh_node'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_build_nodes(fileUnit)
 use IO, only: &
   IO_lc, &
   IO_stringValue, &
   IO_floatValue, &
   IO_stringPos, &
   IO_error, &
   IO_countDataLines, &
   IO_intValue

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) :: i,j,m,c
 logical :: inPart

 allocate ( mesh_node0 (3,mesh_Nnodes), source=0.0_pReal)
 allocate ( mesh_node  (3,mesh_Nnodes), source=0.0_pReal)

610 FORMAT(A300)

 inPart = .false.
 rewind(fileUnit)
 do
   read (fileUnit,610,END=670) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'response' ) &
   ) then
     c = IO_countDataLines(fileUnit)                                                                  ! how many nodes are defined here?
     do i = 1_pInt,c
       backspace(fileUnit)                                                                            ! rewind to first entry
     enddo
     do i = 1_pInt,c
       read (fileUnit,610,END=670) line
       chunkPos = IO_stringPos(line)
       m = mesh_FEasCP('node',IO_intValue(line,chunkPos,1_pInt))
       do j=1_pInt, 3_pInt
         mesh_node0(j,m) = mesh_unitlength * IO_floatValue(line,chunkPos,j+1_pInt)
       enddo
     enddo
   endif
 enddo

670 if (int(size(mesh_node0,2_pInt),pInt) /= mesh_Nnodes) call IO_error(error_ID=909_pInt)
 mesh_node = mesh_node0

end subroutine mesh_abaqus_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets maximum count of nodes, IPs, IP neighbors, and subNodes among cpElements.
!! Sets global values 'mesh_maxNnodes', 'mesh_maxNips', 'mesh_maxNipNeighbors',
!! and 'mesh_maxNcellnodes'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_count_cpSizes(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue ,&
                 IO_error, &
                 IO_countDataLines, &
                 IO_intValue

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) :: i,c,t,g
 logical :: inPart

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNcellnodes   = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2_pInt)),'type'))           ! remember elem type
     g = FE_geomtype(t)
     c = FE_celltype(g)
     mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(t))
     mesh_maxNips =         max(mesh_maxNips,FE_Nips(g))
     mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(c))
     mesh_maxNcellnodes =   max(mesh_maxNcellnodes,FE_Ncellnodes(g))
   endif
 enddo

620 end subroutine mesh_abaqus_count_cpSizes


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, mat, tex, and node list per elemen.
!! Allocates global array 'mesh_element'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_build_elements(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_extractValue, &
                 IO_floatValue, &
                 IO_countDataLines, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos

 integer(pInt) :: i,j,k,c,e,t,homog,micro, nNodesAlreadyRead
 logical inPart,materialFound
 character (len=64) :: materialName,elemSetName
 character(len=300) :: line

 allocate(mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems), source=0_pInt)
 mesh_elemType = -1_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2_pInt)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2_pInt)),'type'))        ! remember elem type
     c = IO_countDataLines(fileUnit)
     do i = 1_pInt,c
       backspace(fileUnit)
     enddo
     do i = 1_pInt,c
       read (fileUnit,610,END=620) line
       chunkPos = IO_stringPos(line)                                                                ! limit to 64 nodes max
       e = mesh_FEasCP('elem',IO_intValue(line,chunkPos,1_pInt))
       if (e /= 0_pInt) then                                                                        ! disregard non CP elems
         mesh_element(1,e) = -1_pInt                                                                ! DEPRECATED
         if (mesh_elemType /= t .and. mesh_elemType /= -1_pInt) &
           call IO_error(191,el=t,ip=mesh_elemType)
         mesh_elemType = t
         mesh_element(2,e) = t                                                                     ! elem type
         nNodesAlreadyRead = 0_pInt
         do j = 1_pInt,chunkPos(1)-1_pInt
           mesh_element(4_pInt+j,e) = mesh_FEasCP('node',IO_intValue(line,chunkPos,1_pInt+j))         ! put CP ids of nodes to position 5:
         enddo
         nNodesAlreadyRead = chunkPos(1) - 1_pInt
         do while(nNodesAlreadyRead < FE_Nnodes(t))                                                ! read on if not all nodes in one line
           read (fileUnit,610,END=620) line
           chunkPos = IO_stringPos(line)
           do j = 1_pInt,chunkPos(1)
             mesh_element(4_pInt+nNodesAlreadyRead+j,e) &
               = mesh_FEasCP('node',IO_IntValue(line,chunkPos,j))                                     ! CP ids of nodes
           enddo
           nNodesAlreadyRead = nNodesAlreadyRead + chunkPos(1)
         enddo
       endif
     enddo
   endif
 enddo


620 rewind(fileUnit)                                                                                 ! just in case "*material" definitions apear before "*element"

 materialFound = .false.
 do
   read (fileUnit,610,END=630) line
   chunkPos = IO_stringPos(line)
   select case ( IO_lc(IO_StringValue(line,chunkPos,1_pInt)))
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_StringValue(line,chunkPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if ( IO_lc(IO_StringValue(line,chunkPos,2_pInt)) == 'material' .and. &
            materialFound ) then
         read (fileUnit,610,END=630) line                                                           ! read homogenization and microstructure
         chunkPos = IO_stringPos(line)
         homog = nint(IO_floatValue(line,chunkPos,1_pInt),pInt)
         micro = nint(IO_floatValue(line,chunkPos,2_pInt),pInt)
         do i = 1_pInt,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) then                                         ! matched?
                 do j = 1_pInt,mesh_mapElemSet(1,k)
                   e = mesh_FEasCP('elem',mesh_mapElemSet(1+j,k))
                   mesh_element(3,e) = homog                                                        ! store homogenization
                   mesh_element(4,e) = micro                                                        ! store microstructure
                   mesh_maxValStateVar(1) = max(mesh_maxValStateVar(1),homog)
                   mesh_maxValStateVar(2) = max(mesh_maxValStateVar(2),micro)
                 enddo
               endif
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo

630 end subroutine mesh_abaqus_build_elements



!--------------------------------------------------------------------------------------------------
!> @brief get any additional damask options from input file, sets mesh_periodicSurface
!--------------------------------------------------------------------------------------------------
subroutine mesh_get_damaskOptions(fileUnit)

use IO, only: &
  IO_lc, &
  IO_stringValue, &
  IO_stringPos

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) chunk, Nchunks
 character(len=300) :: line, damaskOption, v
 character(len=300) :: keyword

 mesh_periodicSurface = .false.
 keyword = '**damask'


 rewind(fileUnit)
 do
   read (fileUnit,610,END=620) line
   chunkPos = IO_stringPos(line)
   Nchunks = chunkPos(1)
   if (IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == keyword .and. Nchunks > 1_pInt) then             ! found keyword for damask option and there is at least one more chunk to read
     damaskOption = IO_lc(IO_stringValue(line,chunkPos,2_pInt))
     select case(damaskOption)
       case('periodic')                                                                             ! damask Option that allows to specify periodic fluxes
         do chunk = 3_pInt,Nchunks                                                                  ! loop through chunks (skipping the keyword)
            v = IO_lc(IO_stringValue(line,chunkPos,chunk))                                             ! chunk matches keyvalues x,y, or z?
            mesh_periodicSurface(1) = mesh_periodicSurface(1) .or. v == 'x'
            mesh_periodicSurface(2) = mesh_periodicSurface(2) .or. v == 'y'
            mesh_periodicSurface(3) = mesh_periodicSurface(3) .or. v == 'z'
         enddo
     endselect
   endif
 enddo

610 FORMAT(A300)

620 end subroutine mesh_get_damaskOptions



!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipAreas
 use math, only: &
   math_crossproduct

 implicit none
 integer(pInt) :: e,t,g,c,i,f,n,m
 real(pReal), dimension (3,FE_maxNcellnodesPerCellface) :: nodePos, normals
 real(pReal), dimension(3) :: normal

 allocate(mesh_ipArea(mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems), source=0.0_pReal)
 allocate(mesh_ipAreaNormal(3_pInt,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems), source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,nodePos,normal,normals)
   do e = 1_pInt,mesh_NcpElems                                                                      ! loop over cpElems
     t = mesh_element(2_pInt,e)                                                                     ! get element type
     g = FE_geomtype(t)                                                                             ! get geometry type
     c = FE_celltype(g)                                                                             ! get cell type
     select case (c)

       case (1_pInt,2_pInt)                                                                         ! 2D 3 or 4 node
         do i = 1_pInt,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal(1) =   nodePos(2,2) - nodePos(2,1)                                              ! x_normal =  y_connectingVector
             normal(2) = -(nodePos(1,2) - nodePos(1,1))                                             ! y_normal = -x_connectingVector
             normal(3) = 0.0_pReal
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (3_pInt)                                                                                ! 3D 4node
         do i = 1_pInt,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal = math_crossproduct(nodePos(1:3,2) - nodePos(1:3,1), &
                                         nodePos(1:3,3) - nodePos(1:3,1))
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (4_pInt)                                                                                ! 3D 8node
         ! for this cell type we get the normal of the quadrilateral face as an average of
         ! four normals of triangular subfaces; since the face consists only of two triangles,
         ! the sum has to be divided by two; this whole prcedure tries to compensate for
         ! probable non-planar cell surfaces
         m = FE_NcellnodesPerCellface(c)
         do i = 1_pInt,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               normals(1:3,n) = 0.5_pReal &
                              * math_crossproduct(nodePos(1:3,1+mod(n  ,m)) - nodePos(1:3,n), &
                                                   nodePos(1:3,1+mod(n+1,m)) - nodePos(1:3,n))
             normal = 0.5_pReal * sum(normals,2)
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)
           enddo
         enddo

     end select
   enddo
 !$OMP END PARALLEL DO

end subroutine mesh_build_ipAreas


!--------------------------------------------------------------------------------------------------
!> @brief assignment of twin nodes for each cp node, allocate globals '_nodeTwins'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_nodeTwins

 implicit none
 integer(pInt) dir, &      ! direction of periodicity
               node, &
               minimumNode, &
               maximumNode, &
               n1, &
               n2
 integer(pInt), dimension(mesh_Nnodes+1) :: minimumNodes, maximumNodes                              ! list of surface nodes (minimum and maximum coordinate value) with first entry giving the number of nodes
 real(pReal)   minCoord, maxCoord, &                                                                ! extreme positions in one dimension
               tolerance                                                                            ! tolerance below which positions are assumed identical
 real(pReal), dimension(3) ::  distance                                                             ! distance between two nodes in all three coordinates
 logical, dimension(mesh_Nnodes) :: unpaired

 allocate(mesh_nodeTwins(3,mesh_Nnodes))
 mesh_nodeTwins = 0_pInt

 tolerance = 0.001_pReal * minval(mesh_ipVolume) ** 0.333_pReal

 do dir = 1_pInt,3_pInt                                                                             ! check periodicity in directions of x,y,z
   if (mesh_periodicSurface(dir)) then                                                              ! only if periodicity is requested


     !*** find out which nodes sit on the surface
     !*** and have a minimum or maximum position in this dimension

     minimumNodes = 0_pInt
     maximumNodes = 0_pInt
     minCoord = minval(mesh_node0(dir,:))
     maxCoord = maxval(mesh_node0(dir,:))
     do node = 1_pInt,mesh_Nnodes                                                                   ! loop through all nodes and find surface nodes
       if (abs(mesh_node0(dir,node) - minCoord) <= tolerance) then
         minimumNodes(1) = minimumNodes(1) + 1_pInt
         minimumNodes(minimumNodes(1)+1_pInt) = node
       elseif (abs(mesh_node0(dir,node) - maxCoord) <= tolerance) then
         maximumNodes(1) = maximumNodes(1) + 1_pInt
         maximumNodes(maximumNodes(1)+1_pInt) = node
       endif
     enddo


     !*** find the corresponding node on the other side with the same position in this dimension

     unpaired = .true.
     do n1 = 1_pInt,minimumNodes(1)
       minimumNode = minimumNodes(n1+1_pInt)
       if (unpaired(minimumNode)) then
         do n2 = 1_pInt,maximumNodes(1)
           maximumNode = maximumNodes(n2+1_pInt)
           distance = abs(mesh_node0(:,minimumNode) - mesh_node0(:,maximumNode))
           if (sum(distance) - distance(dir) <= tolerance) then                                     ! minimum possible distance (within tolerance)
             mesh_nodeTwins(dir,minimumNode) = maximumNode
             mesh_nodeTwins(dir,maximumNode) = minimumNode
             unpaired(maximumNode) = .false.                                                        ! remember this node, we don't have to look for his partner again
             exit
           endif
         enddo
       endif
     enddo

   endif
 enddo

end subroutine mesh_build_nodeTwins


!--------------------------------------------------------------------------------------------------
!> @brief  get maximum count of shared elements among cpElements and build list of elements shared
!! by each node in mesh. Allocate globals '_maxNsharedElems' and '_sharedElem'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_sharedElems

 implicit none
 integer(pint)   e, &                                                                                ! element index
                 g, &                                                                                ! element type
                 node, &                                                                             ! CP node index
                 n, &                                                                                ! node index per element
                 myDim, &                                                                            ! dimension index
                 nodeTwin                                                                            ! node twin in the specified dimension
 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt), dimension(:), allocatable :: node_seen

 allocate(node_seen(maxval(FE_NmatchingNodes)))

 node_count = 0_pInt

 do e = 1_pInt,mesh_NcpElems
   g = FE_geomtype(mesh_element(2,e))                                                                ! get elemGeomType
   node_seen = 0_pInt                                                                                ! reset node duplicates
   do n = 1_pInt,FE_NmatchingNodes(g)                                                                ! check each node of element
     node = mesh_element(4+n,e)
     if (all(node_seen /= node)) then
       node_count(node) = node_count(node) + 1_pInt                                                  ! if FE node not yet encountered -> count it
       do myDim = 1_pInt,3_pInt                                                                      ! check in each dimension...
         nodeTwin = mesh_nodeTwins(myDim,node)
         if (nodeTwin > 0_pInt) &                                                                    ! if I am a twin of some node...
           node_count(nodeTwin) = node_count(nodeTwin) + 1_pInt                                      ! -> count me again for the twin node
       enddo
     endif
     node_seen(n) = node                                                                             ! remember this node to be counted already
   enddo
 enddo

 mesh_maxNsharedElems = int(maxval(node_count),pInt)                                                 ! most shared node

 allocate(mesh_sharedElem(1+mesh_maxNsharedElems,mesh_Nnodes),source=0_pInt)

 do e = 1_pInt,mesh_NcpElems
   g = FE_geomtype(mesh_element(2,e))                                                                ! get elemGeomType
   node_seen = 0_pInt
   do n = 1_pInt,FE_NmatchingNodes(g)
     node = mesh_element(4_pInt+n,e)
     if (all(node_seen /= node)) then
       mesh_sharedElem(1,node) = mesh_sharedElem(1,node) + 1_pInt                                    ! count for each node the connected elements
       mesh_sharedElem(mesh_sharedElem(1,node)+1_pInt,node) = e                                      ! store the respective element id
       do myDim = 1_pInt,3_pInt                                                                      ! check in each dimension...
         nodeTwin = mesh_nodeTwins(myDim,node)
         if (nodeTwin > 0_pInt) then                                                                 ! if i am a twin of some node...
           mesh_sharedElem(1,nodeTwin) = mesh_sharedElem(1,nodeTwin) + 1_pInt                        ! ...count me again for the twin
           mesh_sharedElem(mesh_sharedElem(1,nodeTwin)+1,nodeTwin) = e                               ! store the respective element id
         endif
       enddo
     endif
     node_seen(n) = node
   enddo
 enddo

end subroutine mesh_build_sharedElems


!--------------------------------------------------------------------------------------------------
!> @brief build up of IP neighborhood, allocate globals '_ipNeighborhood'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipNeighborhood
 use math, only: &
   math_mul3x3

 implicit none
 integer(pInt)      ::           myElem, &                                                           ! my CP element index
                                 myIP, &
                                 myType, &                                                           ! my element type
                                 myFace, &
                                 neighbor, &                                                         ! neighor index
                                 neighboringIPkey, &                                                 ! positive integer indicating the neighboring IP (for intra-element) and negative integer indicating the face towards neighbor (for neighboring element)
                                 candidateIP, &
                                 neighboringType, &                                                  ! element type of neighbor
                                 NlinkedNodes, &                                                     ! number of linked nodes
                                 twin_of_linkedNode, &                                               ! node twin of a specific linkedNode
                                 NmatchingNodes, &                                                   ! number of matching nodes
                                 dir, &                                                              ! direction of periodicity
                                 matchingElem, &                                                     ! CP elem number of matching element
                                 matchingFace, &                                                     ! face ID of matching element
                                 a, anchor, &
                                 neighboringIP, &
                                 neighboringElem, &
                                 pointingToMe
 integer(pInt), dimension(FE_maxmaxNnodesAtIP) :: &
                                 linkedNodes = 0_pInt, &
                                 matchingNodes
 logical checkTwins

 allocate(mesh_ipNeighborhood(3,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems))
 mesh_ipNeighborhood = 0_pInt


 do myElem = 1_pInt,mesh_NcpElems                                                                    ! loop over cpElems
   myType = FE_geomtype(mesh_element(2,myElem))                                                      ! get elemGeomType
   do myIP = 1_pInt,FE_Nips(myType)                                                                  ! loop over IPs of elem

     do neighbor = 1_pInt,FE_NipNeighbors(FE_celltype(myType))                                       ! loop over neighbors of IP
       neighboringIPkey = FE_ipNeighbor(neighbor,myIP,myType)

       !*** if the key is positive, the neighbor is inside the element
       !*** that means, we have already found our neighboring IP

       if (neighboringIPkey > 0_pInt) then
         mesh_ipNeighborhood(1,neighbor,myIP,myElem) = myElem
         mesh_ipNeighborhood(2,neighbor,myIP,myElem) = neighboringIPkey


       !*** if the key is negative, the neighbor resides in a neighboring element
       !*** that means, we have to look through the face indicated by the key and see which element is behind that face

       elseif (neighboringIPkey < 0_pInt) then                                                       ! neighboring element's IP
         myFace = -neighboringIPkey
         call mesh_faceMatch(myElem, myFace, matchingElem, matchingFace)                             ! get face and CP elem id of face match
         if (matchingElem > 0_pInt) then                                                             ! found match?
           neighboringType = FE_geomtype(mesh_element(2,matchingElem))

           !*** trivial solution if neighbor has only one IP

           if (FE_Nips(neighboringType) == 1_pInt) then
             mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
             mesh_ipNeighborhood(2,neighbor,myIP,myElem) = 1_pInt
             cycle
           endif

           !*** find those nodes which build the link to the neighbor

           NlinkedNodes = 0_pInt
           linkedNodes = 0_pInt
           do a = 1_pInt,FE_maxNnodesAtIP(myType)                                                    ! figure my anchor nodes on connecting face
             anchor = FE_nodesAtIP(a,myIP,myType)
             if (anchor /= 0_pInt) then                                                              ! valid anchor node
               if (any(FE_face(:,myFace,myType) == anchor)) then                                     ! ip anchor sits on face?
                 NlinkedNodes = NlinkedNodes + 1_pInt
                 linkedNodes(NlinkedNodes) = mesh_element(4_pInt+anchor,myElem)                      ! CP id of anchor node
               else                                                                                  ! something went wrong with the linkage, since not all anchors sit on my face
                 NlinkedNodes = 0_pInt
                 linkedNodes = 0_pInt
                 exit
               endif
             endif
           enddo

           !*** loop through the ips of my neighbor
           !*** and try to find an ip with matching nodes
           !*** also try to match with node twins

 checkCandidateIP: do candidateIP = 1_pInt,FE_Nips(neighboringType)
             NmatchingNodes = 0_pInt
             matchingNodes = 0_pInt
             do a = 1_pInt,FE_maxNnodesAtIP(neighboringType)                                         ! check each anchor node of that ip
               anchor = FE_nodesAtIP(a,candidateIP,neighboringType)
               if (anchor /= 0_pInt) then                                                            ! valid anchor node
                 if (any(FE_face(:,matchingFace,neighboringType) == anchor)) then                    ! sits on matching face?
                   NmatchingNodes = NmatchingNodes + 1_pInt
                   matchingNodes(NmatchingNodes) = mesh_element(4+anchor,matchingElem)               ! CP id of neighbor's anchor node
                 else                                                                                ! no matching, because not all nodes sit on the matching face
                   NmatchingNodes = 0_pInt
                   matchingNodes = 0_pInt
                   exit
                 endif
               endif
             enddo

             if (NmatchingNodes /= NlinkedNodes) &                                                   ! this ip has wrong count of anchors on face
               cycle checkCandidateIP

             !*** check "normal" nodes whether they match or not

             checkTwins = .false.
             do a = 1_pInt,NlinkedNodes
               if (all(matchingNodes /= linkedNodes(a))) then                                        ! this linkedNode does not match any matchingNode
                 checkTwins = .true.
                 exit                                                                                ! no need to search further
               endif
             enddo

             !*** if no match found, then also check node twins

             if(checkTwins) then
               dir = int(maxloc(abs(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem)),1),pInt)            ! check for twins only in direction of the surface normal
               do a = 1_pInt,NlinkedNodes
                 twin_of_linkedNode = mesh_nodeTwins(dir,linkedNodes(a))
                 if (twin_of_linkedNode == 0_pInt .or. &                                             ! twin of linkedNode does not exist...
                     all(matchingNodes /= twin_of_linkedNode)) then                                  ! ... or it does not match any matchingNode
                   cycle checkCandidateIP                                                            ! ... then check next candidateIP
                 endif
               enddo
             endif

             !*** we found a match !!!

             mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
             mesh_ipNeighborhood(2,neighbor,myIP,myElem) = candidateIP
             exit checkCandidateIP
           enddo checkCandidateIP
         endif                                                                                       ! end of valid external matching
       endif                                                                                         ! end of internal/external matching
     enddo
   enddo
 enddo
 do myElem = 1_pInt,mesh_NcpElems                                                                    ! loop over cpElems
   myType = FE_geomtype(mesh_element(2,myElem))                                                      ! get elemGeomType
   do myIP = 1_pInt,FE_Nips(myType)                                                                  ! loop over IPs of elem
     do neighbor = 1_pInt,FE_NipNeighbors(FE_celltype(myType))                                       ! loop over neighbors of IP
       neighboringElem = mesh_ipNeighborhood(1,neighbor,myIP,myElem)
       neighboringIP   = mesh_ipNeighborhood(2,neighbor,myIP,myElem)
       if (neighboringElem > 0_pInt .and. neighboringIP > 0_pInt) then                               ! if neighbor exists ...
         neighboringType = FE_geomtype(mesh_element(2,neighboringElem))
         do pointingToMe = 1_pInt,FE_NipNeighbors(FE_celltype(neighboringType))                      ! find neighboring index that points from my neighbor to myself
           if (    myElem == mesh_ipNeighborhood(1,pointingToMe,neighboringIP,neighboringElem) &
               .and. myIP == mesh_ipNeighborhood(2,pointingToMe,neighboringIP,neighboringElem)) then ! possible candidate
             if (math_mul3x3(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem),&
                             mesh_ipAreaNormal(1:3,pointingToMe,neighboringIP,neighboringElem)) < 0.0_pReal) then ! area normals have opposite orientation (we have to check that because of special case for single element with two ips and periodicity. In this case the neighbor is identical in two different directions.)
               mesh_ipNeighborhood(3,neighbor,myIP,myElem) = pointingToMe                            ! found match
               exit                                                                                  ! so no need to search further
             endif
           endif
         enddo
       endif
     enddo
   enddo
 enddo

end subroutine mesh_build_ipNeighborhood



!--------------------------------------------------------------------------------------------------
!> @brief write statistics regarding input file parsing to the output file
!--------------------------------------------------------------------------------------------------
subroutine mesh_tell_statistics
 use math, only: &
   math_range
 use IO, only: &
   IO_error
 use debug, only: &
   debug_level, &
   debug_MESH, &
   debug_LEVELBASIC, &
   debug_LEVELEXTENSIVE, &
   debug_LEVELSELECTIVE, &
   debug_e, &
   debug_i

 implicit none
 integer(pInt), dimension (:,:), allocatable :: mesh_HomogMicro
 character(len=64) :: myFmt
 integer(pInt) :: i,e,n,f,t,g,c, myDebug

 myDebug = debug_level(debug_mesh)

 if (mesh_maxValStateVar(1) < 1_pInt) call IO_error(error_ID=170_pInt)                              ! no homogenization specified
 if (mesh_maxValStateVar(2) < 1_pInt) call IO_error(error_ID=180_pInt)                              ! no microstructure specified

 allocate (mesh_HomogMicro(mesh_maxValStateVar(1),mesh_maxValStateVar(2)),source = 0_pInt)
 do e = 1_pInt,mesh_NcpElems
   if (mesh_element(3,e) < 1_pInt) call IO_error(error_ID=170_pInt,el=e)                             ! no homogenization specified
   if (mesh_element(4,e) < 1_pInt) call IO_error(error_ID=180_pInt,el=e)                             ! no microstructure specified
   mesh_HomogMicro(mesh_element(3,e),mesh_element(4,e)) = &
   mesh_HomogMicro(mesh_element(3,e),mesh_element(4,e)) + 1_pInt                                     ! count combinations of homogenization and microstructure
 enddo
!$OMP CRITICAL (write2out)
  if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
    write(6,'(/,a,/)') ' Input Parser: STATISTICS'
    write(6,*) mesh_NcpElems,         ' : total number of CP elements in mesh'
    write(6,*) mesh_Nnodes,           ' : total number of nodes in mesh'
    write(6,'(/,a,/)') ' Input Parser: HOMOGENIZATION/MICROSTRUCTURE'
    write(6,*) mesh_maxValStateVar(1), ' : maximum homogenization index'
    write(6,*) mesh_maxValStateVar(2), ' : maximum microstructure index'
    write(6,*)
    write (myFmt,'(a,i32.32,a)') '(9x,a2,1x,',mesh_maxValStateVar(2),'(i8))'
    write(6,myFmt) '+-',math_range(mesh_maxValStateVar(2))
    write (myFmt,'(a,i32.32,a)') '(i8,1x,a2,1x,',mesh_maxValStateVar(2),'(i8))'
    do i=1_pInt,mesh_maxValStateVar(1)                                                              ! loop over all (possibly assigned) homogenizations
      write(6,myFmt) i,'| ',mesh_HomogMicro(i,:)                                                    ! loop over all (possibly assigned) microstructures
    enddo
    write(6,'(/,a,/)') ' Input Parser: ADDITIONAL MPIE OPTIONS'
    write(6,*) 'periodic surface : ', mesh_periodicSurface
    write(6,*)
    flush(6)
  endif

  if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
    write(6,'(/,a,/)') 'Input Parser: ELEMENT TYPE'
    write(6,'(a8,3(1x,a8))') 'elem','elemtype','geomtype','celltype'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = mesh_element(2,e)                      ! get elemType
      g = FE_geomtype(t)                         ! get elemGeomType
      c = FE_celltype(g)                         ! get cellType
      write(6,'(i8,3(1x,i8))') e,t,g,c
    enddo
    write(6,'(/,a)') 'Input Parser: ELEMENT VOLUME'
    write(6,'(/,a13,1x,e15.8)') 'total volume', sum(mesh_ipVolume)
    write(6,'(/,a8,1x,a5,1x,a15,1x,a5,1x,a15,1x,a16)') 'elem','IP','volume','face','area','-- normal --'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = mesh_element(2,e)                                                              ! get element type
      g = FE_geomtype(t)                                                                 ! get geometry type
      c = FE_celltype(g)                                                                 ! get cell type
      do i = 1_pInt,FE_Nips(g)
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write(6,'(i8,1x,i5,1x,e15.8)') e,i,mesh_IPvolume(i,e)
        do f = 1_pInt,FE_NipNeighbors(c)
          write(6,'(i33,1x,e15.8,1x,3(f6.3,1x))') f,mesh_ipArea(f,i,e),mesh_ipAreaNormal(:,f,i,e)
        enddo
      enddo
    enddo
    write(6,'(/,a,/)') 'Input Parser: CELLNODE COORDINATES'
    write(6,'(a8,1x,a2,1x,a8,3(1x,a12))') 'elem','IP','cellnode','x','y','z'
    do e = 1_pInt,mesh_NcpElems                                                          ! loop over cpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = mesh_element(2,e)                                                              ! get element type
      g = FE_geomtype(t)                                                                 ! get geometry type
      c = FE_celltype(g)                                                                 ! get cell type
      do i = 1_pInt,FE_Nips(g)                                                           ! loop over IPs of elem
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write(6,'(i8,1x,i2)') e,i
        do n = 1_pInt,FE_NcellnodesPerCell(c)                                            ! loop over cell nodes in the cell
          write(6,'(12x,i8,3(1x,f12.8))')  mesh_cell(n,i,e), &
                                           mesh_cellnode(1:3,mesh_cell(n,i,e))
        enddo
      enddo
    enddo
    write(6,'(/,a)') 'Input Parser: IP COORDINATES'
    write(6,'(a8,1x,a5,3(1x,a12))') 'elem','IP','x','y','z'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write(6,'(i8,1x,i5,3(1x,f12.8))') e, i, mesh_ipCoordinates(:,i,e)
      enddo
    enddo
    write(6,'(/,a,/)') 'Input Parser: IP NEIGHBORHOOD'
    write(6,'(a8,1x,a10,1x,a10,1x,a3,1x,a13,1x,a13)') 'elem','IP','neighbor','','elemNeighbor','ipNeighbor'
    do e = 1_pInt,mesh_NcpElems                                                          ! loop over cpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = mesh_element(2,e)                                                              ! get element type
      g = FE_geomtype(t)                                                                 ! get geometry type
      c = FE_celltype(g)                                                                 ! get cell type
      do i = 1_pInt,FE_Nips(g)                                                           ! loop over IPs of elem
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        do n = 1_pInt,FE_NipNeighbors(c)                                                 ! loop over neighbors of IP
          write(6,'(i8,1x,i10,1x,i10,1x,a3,1x,i13,1x,i13)') e,i,n,'-->',mesh_ipNeighborhood(1,n,i,e),mesh_ipNeighborhood(2,n,i,e)
        enddo
      enddo
    enddo
  endif
!$OMP END CRITICAL (write2out)

end subroutine mesh_tell_statistics


!--------------------------------------------------------------------------------------------------
!> @brief mapping of FE element types to internal representation
!--------------------------------------------------------------------------------------------------
integer(pInt) function FE_mapElemtype(what)
 use IO, only: IO_lc, IO_error

 implicit none
 character(len=*), intent(in) :: what

 select case (IO_lc(what))
    case ( 'cpe4', &
           'cpe4t')
      FE_mapElemtype = 3_pInt            ! Arbitrary Quadrilateral Plane-strain
    case ( 'cpe8', &
           'cpe8t')
      FE_mapElemtype = 4_pInt            ! Plane Strain, Eight-node Distorted Quadrilateral
    case ( 'c3d4', &
           'c3d4t')
      FE_mapElemtype = 6_pInt            ! Three-dimensional Four-node Tetrahedron
    case ( 'c3d6', &
           'c3d6t')
      FE_mapElemtype = 9_pInt            ! Three-dimensional Arbitrarily Distorted Pentahedral
    case ( 'c3d8r', &
           'c3d8rt')
      FE_mapElemtype = 10_pInt           ! Three-dimensional Arbitrarily Distorted linear hexahedral with reduced integration
    case ( 'c3d8', &
           'c3d8t')
      FE_mapElemtype = 11_pInt           ! Three-dimensional Arbitrarily Distorted Brick
    case ( 'c3d20r', &
           'c3d20rt')
      FE_mapElemtype = 12_pInt           ! Three-dimensional Arbitrarily Distorted quad hexahedral with reduced integration
    case ( 'c3d20', &
           'c3d20t')
      FE_mapElemtype = 13_pInt           ! Three-dimensional Arbitrarily Distorted quadratic hexahedral
    case default
      call IO_error(error_ID=190_pInt,ext_msg=IO_lc(what))
 end select

end function FE_mapElemtype


!--------------------------------------------------------------------------------------------------
!> @brief find face-matching element of same type
!--------------------------------------------------------------------------------------------------
subroutine mesh_faceMatch(elem, face ,matchingElem, matchingFace)

implicit none
integer(pInt), intent(out) ::     matchingElem, &                                                   ! matching CP element ID
                                  matchingFace                                                      ! matching face ID
integer(pInt), intent(in) ::      face, &                                                           ! face ID
                                  elem                                                              ! CP elem ID
integer(pInt), dimension(FE_NmatchingNodesPerFace(face,FE_geomtype(mesh_element(2,elem)))) :: &
                                  myFaceNodes                                                       ! global node ids on my face
integer(pInt)        ::           myType, &
                                  candidateType, &
                                  candidateElem, &
                                  candidateFace, &
                                  candidateFaceNode, &
                                  minNsharedElems, &
                                  NsharedElems, &
                                  lonelyNode = 0_pInt, &
                                  i, &
                                  n, &
                                  dir                                                               ! periodicity direction
integer(pInt), dimension(:), allocatable :: element_seen
logical checkTwins

matchingElem = 0_pInt
matchingFace = 0_pInt
minNsharedElems = mesh_maxNsharedElems + 1_pInt                                                     ! init to worst case
myType = FE_geomtype(mesh_element(2_pInt,elem))                                                     ! figure elemGeomType

do n = 1_pInt,FE_NmatchingNodesPerFace(face,myType)                                                 ! loop over nodes on face
  myFaceNodes(n) = mesh_element(4_pInt+FE_face(n,face,myType),elem)                                 ! CP id of face node
  NsharedElems = mesh_sharedElem(1_pInt,myFaceNodes(n))                                             ! figure # shared elements for this node
  if (NsharedElems < minNsharedElems) then
    minNsharedElems = NsharedElems                                                                  ! remember min # shared elems
    lonelyNode = n                                                                                  ! remember most lonely node
  endif
enddo

allocate(element_seen(minNsharedElems))
element_seen = 0_pInt

checkCandidate: do i = 1_pInt,minNsharedElems                                                       ! iterate over lonelyNode's shared elements
  candidateElem = mesh_sharedElem(1_pInt+i,myFaceNodes(lonelyNode))                                 ! present candidate elem
  if (all(element_seen /= candidateElem)) then                                                      ! element seen for the first time?
    element_seen(i) = candidateElem
    candidateType = FE_geomtype(mesh_element(2_pInt,candidateElem))                                 ! figure elemGeomType of candidate
checkCandidateFace: do candidateFace = 1_pInt,FE_maxNipNeighbors                                    ! check each face of candidate
      if (FE_NmatchingNodesPerFace(candidateFace,candidateType) &
          /= FE_NmatchingNodesPerFace(face,myType) &                                                ! incompatible face
          .or. (candidateElem == elem .and. candidateFace == face)) then                            ! this is my face
        cycle checkCandidateFace
      endif
      checkTwins = .false.
      do n = 1_pInt,FE_NmatchingNodesPerFace(candidateFace,candidateType)                           ! loop through nodes on face
        candidateFaceNode = mesh_element(4_pInt+FE_face(n,candidateFace,candidateType),candidateElem)
        if (all(myFaceNodes /= candidateFaceNode)) then                                             ! candidate node does not match any of my face nodes
          checkTwins = .true.                                                                       ! perhaps the twin nodes do match
          exit
        endif
      enddo
      if(checkTwins) then
checkCandidateFaceTwins: do dir = 1_pInt,3_pInt
          do n = 1_pInt,FE_NmatchingNodesPerFace(candidateFace,candidateType)                       ! loop through nodes on face
            candidateFaceNode = mesh_element(4+FE_face(n,candidateFace,candidateType),candidateElem)
            if (all(myFaceNodes /= mesh_nodeTwins(dir,candidateFaceNode))) then                     ! node twin does not match either
              if (dir == 3_pInt) then
                cycle checkCandidateFace
              else
                cycle checkCandidateFaceTwins                                                       ! try twins in next dimension
              endif
            endif
          enddo
          exit checkCandidateFaceTwins
        enddo checkCandidateFaceTwins
      endif
      matchingFace = candidateFace
      matchingElem = candidateElem
      exit checkCandidate                                                                           ! found my matching candidate
    enddo checkCandidateFace
  endif
enddo checkCandidate

end subroutine mesh_faceMatch


!--------------------------------------------------------------------------------------------------
!> @brief get properties of different types of finite elements
!> @details assign globals: FE_nodesAtIP, FE_ipNeighbor, FE_cellnodeParentnodeWeights, FE_subNodeOnIPFace
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_FEdata

 implicit none
 integer(pInt) :: me
 allocate(FE_nodesAtIP(FE_maxmaxNnodesAtIP,FE_maxNips,FE_Ngeomtypes), source=0_pInt)
 allocate(FE_ipNeighbor(FE_maxNipNeighbors,FE_maxNips,FE_Ngeomtypes), source=0_pInt)
 allocate(FE_cell(FE_maxNcellnodesPerCell,FE_maxNips,FE_Ngeomtypes),  source=0_pInt)
 allocate(FE_cellnodeParentnodeWeights(FE_maxNnodes,FE_maxNcellnodes,FE_Nelemtypes), source=0.0_pReal)
 allocate(FE_cellface(FE_maxNcellnodesPerCellface,FE_maxNcellfaces,FE_Ncelltypes),  source=0_pInt)


 !*** fill FE_nodesAtIP with data ***

 me = 0_pInt

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    1,2,3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
    1,  &
    2,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  27 (2D 8node 9ip)
    reshape(int([&
    1,0,  &
    1,2,  &
    2,0,  &
    1,4,  &
    0,0,  &
    2,3,  &
    4,0,  &
    3,4,  &
    3,0   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
    1,2,3,4   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4,  &
    5,  &
    6   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
    1,2,3,4,5,6,7,8   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3,  &
    5,  &
    6,  &
    8,  &
    7   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1_pInt
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  21 (3D 20node 27ip)
    reshape(int([&
    1,0, 0,0,  &
    1,2, 0,0,  &
    2,0, 0,0,  &
    1,4, 0,0,  &
    1,3, 2,4,  &
    2,3, 0,0,  &
    4,0, 0,0,  &
    3,4, 0,0,  &
    3,0, 0,0,  &
    1,5, 0,0,  &
    1,6, 2,5,  &
    2,6, 0,0,  &
    1,8, 4,5,  &
    0,0, 0,0,  &
    2,7, 3,6,  &
    4,8, 0,0,  &
    3,8, 4,7,  &
    3,7, 0,0,  &
    5,0, 0,0,  &
    5,6, 0,0,  &
    6,0, 0,0,  &
    5,8, 0,0,  &
    5,7, 6,8,  &
    6,7, 0,0,  &
    8,0, 0,0,  &
    7,8, 0,0,  &
    7,0, 0,0   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])


 ! *** FE_ipNeighbor ***
 ! is a list of the neighborhood of each IP.
 ! It is sorted in (local) +x,-x, +y,-y, +z,-z direction.
 ! Positive integers denote an intra-FE IP identifier.
 ! Negative integers denote the interface behind which the neighboring (extra-FE) IP will be located.
 me = 0_pInt

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    -2,-3,-1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
     2,-3, 3,-1,  &
    -2, 1, 3,-1,  &
     2,-3,-2, 1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
     2,-4, 3,-1,  &
    -2, 1, 4,-1,  &
     4,-4,-3, 1,  &
    -2, 3,-3, 2   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element  27 (2D 8node 9ip)
    reshape(int([&
     2,-4, 4,-1,  &
     3, 1, 5,-1,  &
    -2, 2, 6,-1,  &
     5,-4, 7, 1,  &
     6, 4, 8, 2,  &
    -2, 5, 9, 3,  &
     8,-4,-3, 4,  &
     9, 7,-3, 5,  &
    -2, 8,-3, 6   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
    -1,-2,-3,-4   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -2, 1, 3,-2, 4,-1,  &
     2,-4,-3, 1, 4,-1,  &
     2,-4, 3,-2,-3, 1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -3, 1, 3,-2, 5,-1,  &
     2,-4,-3, 1, 6,-1,  &
     5,-4, 6,-2,-5, 1,  &
    -3, 4, 6,-2,-5, 2,  &
     5,-4,-3, 4,-5, 3   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
    -3,-5,-4,-2,-6,-1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
     2,-5, 3,-2, 5,-1,  &
    -3, 1, 4,-2, 6,-1,  &
     4,-5,-4, 1, 7,-1,  &
    -3, 3,-4, 2, 8,-1,  &
     6,-5, 7,-2,-6, 1,  &
    -3, 5, 8,-2,-6, 2,  &
     8,-5,-4, 5,-6, 3,  &
    -3, 7,-4, 6,-6, 4   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element  21 (3D 20node 27ip)
    reshape(int([&
     2,-5, 4,-2,10,-1,  &
     3, 1, 5,-2,11,-1,  &
    -3, 2, 6,-2,12,-1,  &
     5,-5, 7, 1,13,-1,  &
     6, 4, 8, 2,14,-1,  &
    -3, 5, 9, 3,15,-1,  &
     8,-5,-4, 4,16,-1,  &
     9, 7,-4, 5,17,-1,  &
    -3, 8,-4, 6,18,-1,  &
    11,-5,13,-2,19, 1,  &
    12,10,14,-2,20, 2,  &
    -3,11,15,-2,21, 3,  &
    14,-5,16,10,22, 4,  &
    15,13,17,11,23, 5,  &
    -3,14,18,12,24, 6,  &
    17,-5,-4,13,25, 7,  &
    18,16,-4,14,26, 8,  &
    -3,17,-4,15,27, 9,  &
    20,-5,22,-2,-6,10,  &
    21,19,23,-2,-6,11,  &
    -3,20,24,-2,-6,12,  &
    23,-5,25,19,-6,13,  &
    24,22,26,20,-6,14,  &
    -3,23,27,21,-6,15,  &
    26,-5,-4,22,-6,16,  &
    27,25,-4,23,-6,17,  &
    -3,26,-4,24,-6,18   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])


 ! *** FE_cell ***
 me = 0_pInt

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    1,2,3   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   125 (2D 6node 3ip)
    reshape(int([&
    1, 4, 7, 6,   &
    2, 5, 7, 4,   &
    3, 6, 7, 5    &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   11 (2D 4node 4ip)
    reshape(int([&
    1, 5, 9, 8,   &
    5, 2, 6, 9,   &
    8, 9, 7, 4,   &
    9, 6, 3, 7    &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   27 (2D 8node 9ip)
    reshape(int([&
    1, 5,13,12,   &
    5, 6,14,13,   &
    6, 2, 7,14,   &
   12,13,16,11,   &
   13,14,15,16,   &
   14, 7, 8,15,   &
   11,16,10, 4,   &
   16,15, 9,10,   &
   15, 8, 3, 9    &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   134 (3D 4node 1ip)
    reshape(int([&
    1, 2, 3, 4   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   127 (3D 10node 4ip)
    reshape(int([&
    1, 5,11, 7, 8,12,15,14,  &
    5, 2, 6,11,12, 9,13,15,  &
    7,11, 6, 3,14,15,13,10,  &
    8,12,15, 4, 4, 9,13,10   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   136 (3D 6node 6ip)
    reshape(int([&
    1, 7,16, 9,10,17,21,19,  &
    7, 2, 8,16,17,11,18,21,  &
    9,16, 8, 3,19,21,18,12,  &
   10,17,21,19, 4,13,20,15,  &
   17,11,18,21,13, 5,14,20,  &
   19,21,18,12,15,20,14, 6   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   117 (3D 8node 1ip)
    reshape(int([&
    1, 2, 3, 4, 5, 6, 7, 8   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   7 (3D 8node 8ip)
    reshape(int([&
    1, 9,21,12,13,22,27,25,  &
    9, 2,10,21,22,14,23,27,  &
   12,21,11, 4,25,27,24,16,  &
   21,10, 3,11,27,23,15,24,  &
   13,22,27,25, 5,17,26,20,  &
   22,14,23,27,17, 6,18,26,  &
   25,27,24,16,20,26,19, 8,  &
   27,23,15,24,26,18, 7,19   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1_pInt
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   21 (3D 20node 27ip)
    reshape(int([&
    1, 9,33,16,17,37,57,44,  &
    9,10,34,33,37,38,58,57,  &
   10, 2,11,34,38,18,39,58,  &
   16,33,36,15,44,57,60,43,  &
   33,34,35,36,57,58,59,60,  &
   34,11,12,35,58,39,40,59,  &
   15,36,14, 4,43,60,42,20,  &
   36,35,13,14,60,59,41,42,  &
   35,12, 3,13,59,40,19,41,  &
   17,37,57,44,21,45,61,52,  &
   37,38,58,57,45,46,62,61,  &
   38,18,39,58,46,22,47,62,  &
   44,57,60,43,52,61,64,51,  &
   57,58,59,60,61,62,63,64,  &
   58,39,40,59,62,47,48,63,  &
   43,60,42,20,51,64,50,24,  &
   60,59,41,42,64,63,49,50,  &
   59,40,19,41,63,48,23,49,  &
   21,45,61,52, 5,25,53,32,  &
   45,46,62,61,25,26,54,53,  &
   46,22,47,62,26, 6,27,54,  &
   52,61,64,51,32,53,56,31,  &
   61,62,63,64,53,54,55,56,  &
   62,47,48,63,54,27,28,55,  &
   51,64,50,24,31,56,30, 8,  &
   64,63,49,50,56,55,29,30,  &
   63,48,23,49,55,28, 7,29   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])


 ! *** FE_cellnodeParentnodeWeights ***
 ! center of gravity of the weighted nodes gives the position of the cell node.
 ! fill with 0.
 ! example: face-centered cell node with face nodes 1,2,5,6 to be used in,
 !          e.g., an 8 node element, would be encoded:
 !          1, 1, 0, 0, 1, 1, 0, 0
 me = 0_pInt

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element   6 (2D 3node 1ip)
    reshape(real([&
    1, 0, 0,  &
    0, 1, 0,  &
    0, 0, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 125 (2D 6node 3ip)
    reshape(real([&
    1, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0,  &
    0, 0, 0, 0, 1, 0,  &
    0, 0, 0, 0, 0, 1,  &
    1, 1, 1, 2, 2, 2   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element  11 (2D 4node 4ip)
    reshape(real([&
    1, 0, 0, 0,  &
    0, 1, 0, 0,  &
    0, 0, 1, 0,  &
    0, 0, 0, 1,  &
    1, 1, 0, 0,  &
    0, 1, 1, 0,  &
    0, 0, 1, 1,  &
    1, 0, 0, 1,  &
    1, 1, 1, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element  27 (2D 8node 9ip)
    reshape(real([&
    1, 0, 0, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0, 0, 0,  &
    1, 0, 0, 0, 2, 0, 0, 0,  &
    0, 1, 0, 0, 2, 0, 0, 0,  &
    0, 1, 0, 0, 0, 2, 0, 0,  &
    0, 0, 1, 0, 0, 2, 0, 0,  &
    0, 0, 1, 0, 0, 0, 2, 0,  &
    0, 0, 0, 1, 0, 0, 2, 0,  &
    0, 0, 0, 1, 0, 0, 0, 2,  &
    1, 0, 0, 0, 0, 0, 0, 2,  &
    4, 1, 1, 1, 8, 2, 2, 8,  &
    1, 4, 1, 1, 8, 8, 2, 2,  &
    1, 1, 4, 1, 2, 8, 8, 2,  &
    1, 1, 1, 4, 2, 2, 8, 8   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element  54 (2D 8node 4ip)
    reshape(real([&
    1, 0, 0, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0, 0, 0,  &
    0, 0, 0, 0, 1, 0, 0, 0,  &
    0, 0, 0, 0, 0, 1, 0, 0,  &
    0, 0, 0, 0, 0, 0, 1, 0,  &
    0, 0, 0, 0, 0, 0, 0, 1,  &
    1, 1, 1, 1, 2, 2, 2, 2   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 134 (3D 4node 1ip)
    reshape(real([&
    1, 0, 0, 0,  &
    0, 1, 0, 0,  &
    0, 0, 1, 0,  &
    0, 0, 0, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 157 (3D 5node 4ip)
    reshape(real([&
    1, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0,  &
    0, 0, 1, 0, 0,  &
    0, 0, 0, 1, 0,  &
    1, 1, 0, 0, 0,  &
    0, 1, 1, 0, 0,  &
    1, 0, 1, 0, 0,  &
    1, 0, 0, 1, 0,  &
    0, 1, 0, 1, 0,  &
    0, 0, 1, 1, 0,  &
    1, 1, 1, 0, 0,  &
    1, 1, 0, 1, 0,  &
    0, 1, 1, 1, 0,  &
    1, 0, 1, 1, 0,  &
    0, 0, 0, 0, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 127 (3D 10node 4ip)
    reshape(real([&
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  &
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  &
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  &
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  &
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  &
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  &
    1, 1, 1, 0, 2, 2, 2, 0, 0, 0,  &
    1, 1, 0, 1, 2, 0, 0, 2, 2, 0,  &
    0, 1, 1, 1, 0, 2, 0, 0, 2, 2,  &
    1, 0, 1, 1, 0, 0, 2, 2, 0, 2,  &
    3, 3, 3, 3, 4, 4, 4, 4, 4, 4   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 136 (3D 6node 6ip)
    reshape(real([&
    1, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0,  &
    0, 0, 0, 0, 1, 0,  &
    0, 0, 0, 0, 0, 1,  &
    1, 1, 0, 0, 0, 0,  &
    0, 1, 1, 0, 0, 0,  &
    1, 0, 1, 0, 0, 0,  &
    1, 0, 0, 1, 0, 0,  &
    0, 1, 0, 0, 1, 0,  &
    0, 0, 1, 0, 0, 1,  &
    0, 0, 0, 1, 1, 0,  &
    0, 0, 0, 0, 1, 1,  &
    0, 0, 0, 1, 0, 1,  &
    1, 1, 1, 0, 0, 0,  &
    1, 1, 0, 1, 1, 0,  &
    0, 1, 1, 0, 1, 1,  &
    1, 0, 1, 1, 0, 1,  &
    0, 0, 0, 1, 1, 1,  &
    1, 1, 1, 1, 1, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 117 (3D 8node 1ip)
    reshape(real([&
    1, 0, 0, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0, 0, 0,  &
    0, 0, 0, 0, 1, 0, 0, 0,  &
    0, 0, 0, 0, 0, 1, 0, 0,  &
    0, 0, 0, 0, 0, 0, 1, 0,  &
    0, 0, 0, 0, 0, 0, 0, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element   7 (3D 8node 8ip)
    reshape(real([&
    1, 0, 0, 0,  0, 0, 0, 0,  &   !
    0, 1, 0, 0,  0, 0, 0, 0,  &   !
    0, 0, 1, 0,  0, 0, 0, 0,  &   !
    0, 0, 0, 1,  0, 0, 0, 0,  &   !
    0, 0, 0, 0,  1, 0, 0, 0,  &   !  5
    0, 0, 0, 0,  0, 1, 0, 0,  &   !
    0, 0, 0, 0,  0, 0, 1, 0,  &   !
    0, 0, 0, 0,  0, 0, 0, 1,  &   !
    1, 1, 0, 0,  0, 0, 0, 0,  &   !
    0, 1, 1, 0,  0, 0, 0, 0,  &   ! 10
    0, 0, 1, 1,  0, 0, 0, 0,  &   !
    1, 0, 0, 1,  0, 0, 0, 0,  &   !
    1, 0, 0, 0,  1, 0, 0, 0,  &   !
    0, 1, 0, 0,  0, 1, 0, 0,  &   !
    0, 0, 1, 0,  0, 0, 1, 0,  &   ! 15
    0, 0, 0, 1,  0, 0, 0, 1,  &   !
    0, 0, 0, 0,  1, 1, 0, 0,  &   !
    0, 0, 0, 0,  0, 1, 1, 0,  &   !
    0, 0, 0, 0,  0, 0, 1, 1,  &   !
    0, 0, 0, 0,  1, 0, 0, 1,  &   ! 20
    1, 1, 1, 1,  0, 0, 0, 0,  &   !
    1, 1, 0, 0,  1, 1, 0, 0,  &   !
    0, 1, 1, 0,  0, 1, 1, 0,  &   !
    0, 0, 1, 1,  0, 0, 1, 1,  &   !
    1, 0, 0, 1,  1, 0, 0, 1,  &   ! 25
    0, 0, 0, 0,  1, 1, 1, 1,  &   !
    1, 1, 1, 1,  1, 1, 1, 1   &   !
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element  57 (3D 20node 8ip)
    reshape(real([&
    1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !  5
    0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   ! 10
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 1, 0, &   ! 15
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 1, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0, &   ! 20
    1, 1, 1, 1,  0, 0, 0, 0,  2, 2, 2, 2,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    1, 1, 0, 0,  1, 1, 0, 0,  2, 0, 0, 0,  2, 0, 0, 0,  2, 2, 0, 0, &   !
    0, 1, 1, 0,  0, 1, 1, 0,  0, 2, 0, 0,  0, 2, 0, 0,  0, 2, 2, 0, &   !
    0, 0, 1, 1,  0, 0, 1, 1,  0, 0, 2, 0,  0, 0, 2, 0,  0, 0, 2, 2, &   !
    1, 0, 0, 1,  1, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 2,  2, 0, 0, 2, &   ! 25
    0, 0, 0, 0,  1, 1, 1, 1,  0, 0, 0, 0,  2, 2, 2, 2,  0, 0, 0, 0, &   !
    3, 3, 3, 3,  3, 3, 3, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4  &   !
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1_pInt
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element  21 (3D 20node 27ip)
    reshape(real([&
    1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !  5
    0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    1, 0, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 1, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   ! 10
    0, 1, 0, 0,  0, 0, 0, 0,  0, 2, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 1, 0,  0, 0, 0, 0,  0, 2, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0, &   ! 15
    1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0, &   !
    0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 2, 0, 0, &   !
    0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 2, 0, &   !
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2, &   ! 20
    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 2, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 2, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2, &   !
    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0, &   ! 25
    0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  0, 2, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 2, 0, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0, &   ! 30
    0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0, &   !
    4, 1, 1, 1,  0, 0, 0, 0,  8, 2, 2, 8,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    1, 4, 1, 1,  0, 0, 0, 0,  8, 8, 2, 2,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    1, 1, 4, 1,  0, 0, 0, 0,  2, 8, 8, 2,  0, 0, 0, 0,  0, 0, 0, 0, &   ! 35
    1, 1, 1, 4,  0, 0, 0, 0,  2, 2, 8, 8,  0, 0, 0, 0,  0, 0, 0, 0, &   !
    4, 1, 0, 0,  1, 1, 0, 0,  8, 0, 0, 0,  2, 0, 0, 0,  8, 2, 0, 0, &   !
    1, 4, 0, 0,  1, 1, 0, 0,  8, 0, 0, 0,  2, 0, 0, 0,  2, 8, 0, 0, &   !
    0, 4, 1, 0,  0, 1, 1, 0,  0, 8, 0, 0,  0, 2, 0, 0,  0, 8, 2, 0, &   !
    0, 1, 4, 0,  0, 1, 1, 0,  0, 8, 0, 0,  0, 2, 0, 0,  0, 2, 8, 0, &   ! 40
    0, 0, 4, 1,  0, 0, 1, 1,  0, 0, 8, 0,  0, 0, 2, 0,  0, 0, 8, 2, &   !
    0, 0, 1, 4,  0, 0, 1, 1,  0, 0, 8, 0,  0, 0, 2, 0,  0, 0, 2, 8, &   !
    1, 0, 0, 4,  1, 0, 0, 1,  0, 0, 0, 8,  0, 0, 0, 2,  2, 0, 0, 8, &   !
    4, 0, 0, 1,  1, 0, 0, 1,  0, 0, 0, 8,  0, 0, 0, 2,  8, 0, 0, 2, &   !
    1, 1, 0, 0,  4, 1, 0, 0,  2, 0, 0, 0,  8, 0, 0, 0,  8, 2, 0, 0, &   ! 45
    1, 1, 0, 0,  1, 4, 0, 0,  2, 0, 0, 0,  8, 0, 0, 0,  2, 8, 0, 0, &   !
    0, 1, 1, 0,  0, 4, 1, 0,  0, 2, 0, 0,  0, 8, 0, 0,  0, 8, 2, 0, &   !
    0, 1, 1, 0,  0, 1, 4, 0,  0, 2, 0, 0,  0, 8, 0, 0,  0, 2, 8, 0, &   !
    0, 0, 1, 1,  0, 0, 4, 1,  0, 0, 2, 0,  0, 0, 8, 0,  0, 0, 8, 2, &   !
    0, 0, 1, 1,  0, 0, 1, 4,  0, 0, 2, 0,  0, 0, 8, 0,  0, 0, 2, 8, &   ! 50
    1, 0, 0, 1,  1, 0, 0, 4,  0, 0, 0, 2,  0, 0, 0, 8,  2, 0, 0, 8, &   !
    1, 0, 0, 1,  4, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 8,  8, 0, 0, 2, &   !
    0, 0, 0, 0,  4, 1, 1, 1,  0, 0, 0, 0,  8, 2, 2, 8,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  1, 4, 1, 1,  0, 0, 0, 0,  8, 8, 2, 2,  0, 0, 0, 0, &   !
    0, 0, 0, 0,  1, 1, 4, 1,  0, 0, 0, 0,  2, 8, 8, 2,  0, 0, 0, 0, &   ! 55
    0, 0, 0, 0,  1, 1, 1, 4,  0, 0, 0, 0,  2, 2, 8, 8,  0, 0, 0, 0, &   !
   24, 8, 4, 8,  8, 4, 3, 4, 32,12,12,32, 12, 4, 4,12, 32,12, 4,12, &   !
    8,24, 8, 4,  4, 8, 4, 3, 32,32,12,12, 12,12, 4, 4, 12,32,12, 4, &   !
    4, 8,24, 8,  3, 4, 8, 4, 12,32,32,12,  4,12,12, 4,  4,12,32,12, &   !
    8, 4, 8,24,  4, 3, 4, 8, 12,12,32,32,  4, 4,12,12, 12, 4,12,32, &   ! 60
    8, 4, 3, 4, 24, 8, 4, 8, 12, 4, 4,12, 32,12,12,32, 32,12, 4,12, &   !
    4, 8, 4, 3,  8,24, 8, 4, 12,12, 4, 4, 32,32,12,12, 12,32,12, 4, &   !
    3, 4, 8, 4,  4, 8,24, 8,  4,12,12, 4, 12,32,32,12,  4,12,32,12, &   !
    4, 3, 4, 8,  8, 4, 8,24,  4, 4,12,12, 12,12,32,32, 12, 4,12,32  &   !
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])



 ! *** FE_cellface ***
 me = 0_pInt

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 2D 3node, VTK_TRIANGLE (5)
    reshape(int([&
    2,3,  &
    3,1,  &
    1,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 2D 4node, VTK_QUAD (9)
    reshape(int([&
    2,3,  &
    4,1,  &
    3,4,  &
    1,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 3D 4node, VTK_TETRA (10)
    reshape(int([&
    1,3,2,  &
    1,2,4,  &
    2,3,4,  &
    1,4,3   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1_pInt
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 3D 8node, VTK_HEXAHEDRON (12)
    reshape(int([&
    2,3,7,6,  &
    4,1,5,8,  &
    3,4,8,7,  &
    1,2,6,5,  &
    5,6,7,8,  &
    1,4,3,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])


end subroutine mesh_build_FEdata

end module mesh
