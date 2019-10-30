!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solvers MSC.Marc, Abaqus and the spectral solver
!--------------------------------------------------------------------------------------------------
module mesh
 use prec
 use mesh_base
 use geometry_plastic_nonlocal
 use discretization
 use math

 implicit none
 private
 
 integer, public, protected :: &
   mesh_NcpElems, &                                                                                 !< total number of CP elements in local mesh
   mesh_elemType, &                                                                                 !< Element type of the mesh (only support homogeneous meshes)
   mesh_Nnodes, &                                                                                   !< total number of nodes in mesh
   mesh_Ncellnodes, &                                                                               !< total number of cell nodes in mesh (including duplicates)
   mesh_Ncells, &                                                                                   !< total number of cells in mesh
   mesh_maxNipNeighbors, &                                                                          !< max number of IP neighbors in any CP element
   mesh_maxNsharedElems                                                                             !< max number of CP elements sharing a node
!!!! BEGIN DEPRECATED !!!!!
 integer, public, protected :: &
   mesh_maxNips, &                                                                                  !< max number of IPs in any CP element
   mesh_maxNcellnodes                                                                               !< max number of cell nodes in any CP element
!!!! BEGIN DEPRECATED !!!!!

 integer, dimension(:,:), allocatable, public, protected :: &
   mesh_element, & !DEPRECATED
   mesh_sharedElem, &                                                                               !< entryCount and list of elements containing node
   mesh_nodeTwins                                                                                   !< node twins are surface nodes that lie exactly on opposite sides of the mesh (surfaces nodes with equal coordinate values in two dimensions)

 integer, dimension(:,:,:,:), allocatable, public, protected :: &
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

 integer, private :: &
   mesh_maxNelemInSet, &
   mesh_Nmaterials

 integer, dimension(2), private :: &
   mesh_maxValStateVar = 0

integer, dimension(:,:), allocatable, private :: &
   mesh_cellnodeParent                                                                              !< cellnode's parent element ID, cellnode's intra-element ID

 integer,dimension(:,:,:), allocatable, private :: &
   mesh_cell                                                                                        !< cell connectivity for each element,ip/cell

 integer, dimension(:,:,:), allocatable, private :: &
   FE_nodesAtIP, &                                                                                  !< map IP index to node indices in a specific type of element
   FE_ipNeighbor, &                                                                                 !< +x,-x,+y,-y,+z,-z list of intra-element IPs and(negative) neighbor faces per own IP in a specific type of element
   FE_cell, &                                                                                       !< list of intra-element cell node IDs that constitute the cells in a specific type of element geometry
   FE_cellface                                                                                      !< list of intra-cell cell node IDs that constitute the cell faces of a specific type of cell

 real(pReal), dimension(:,:,:), allocatable, private :: &
   FE_cellnodeParentnodeWeights                                                                     !< list of node weights for the generation of cell nodes

 integer, dimension(:,:,:,:), allocatable, private :: &
   FE_subNodeOnIPFace

! These definitions should actually reside in the FE-solver specific part (different for MARC/ABAQUS)
! Hence, I suggest to prefix with "FE_"

 integer, parameter, public :: &
   FE_Nelemtypes = 13, &
   FE_Ngeomtypes = 10, &
   FE_Ncelltypes = 4, &
   FE_maxNnodes = 20, &
   FE_maxNips = 27, &
   FE_maxNipNeighbors = 6, &
   FE_maxmaxNnodesAtIP = 8, &                                                                       !< max number of (equivalent) nodes attached to an IP
   FE_maxNmatchingNodesPerFace = 4, &
   FE_maxNfaces = 6, &
   FE_maxNcellnodes = 64, &
   FE_maxNcellnodesPerCell = 8, &
   FE_maxNcellfaces = 6, &
   FE_maxNcellnodesPerCellface = 4

 integer, dimension(FE_Nelemtypes), parameter, public :: FE_geomtype = &                           !< geometry type of particular element type
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

 integer, dimension(FE_Ngeomtypes), parameter, public  :: FE_celltype = &                          !< cell type that is used by each geometry type
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

 integer, dimension(FE_Ngeomtypes), parameter, public :: FE_dimension = &                          !< dimension of geometry type
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

 integer, dimension(FE_Nelemtypes), parameter, public :: FE_Nnodes = &                             !< number of nodes that constitute a specific type of element
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

 integer, dimension(FE_Ngeomtypes), parameter, public :: FE_Nfaces = &                             !< number of faces of a specific type of element geometry
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

 integer, dimension(FE_Ngeomtypes), parameter, private :: FE_NmatchingNodes = &                    !< number of nodes that are needed for face matching in a specific type of element geometry
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

 integer, dimension(FE_maxNfaces,FE_Ngeomtypes), parameter, private :: FE_NmatchingNodesPerFace = & !< number of matching nodes per face in a specific type of element geometry
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

 integer, dimension(FE_maxNmatchingNodesPerFace,FE_maxNfaces,FE_Ngeomtypes), &
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

 integer, dimension(FE_Ngeomtypes), parameter, private :: FE_Ncellnodes = &                   !< number of cell nodes in a specific geometry type
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

 integer, dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCell = &             !< number of cell nodes in a specific cell type
 int([ &
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      8  & ! (3D 8node)
  ],pInt)

 integer, dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCellface = &        !< number of cell nodes per cell face in a specific cell type
 int([&
      2, & ! (2D 3node)
      2, & ! (2D 4node)
      3, & ! (3D 4node)
      4  & ! (3D 8node)
  ],pInt)

 integer, dimension(FE_Ngeomtypes), parameter, public :: FE_Nips = &                          !< number of IPs in a specific type of element
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

 integer, dimension(FE_Ncelltypes), parameter, public :: FE_NipNeighbors = &                  !< number of ip neighbors / cell faces in a specific cell type
 int([&
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      6  & ! (3D 8node)
  ],pInt)

 integer, dimension(FE_Ngeomtypes), parameter, private :: FE_maxNnodesAtIP = &                !< maximum number of parent nodes that belong to an IP for a specific type of element
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

 integer, private :: &
   mesh_Nelems, &                                                                                   !< total number of elements in mesh (including non-DAMASK elements)
   mesh_maxNnodes, &                                                                                !< max number of nodes in any CP element
   mesh_NelemSets
 character(len=64), dimension(:), allocatable, private :: &
   mesh_nameElemSet, &                                                                              !< names of elementSet
   mesh_nameMaterial, &                                                                             !< names of material in solid section
   mesh_mapMaterial                                                                                 !< name of elementSet for material
 integer, dimension(:,:), allocatable, private :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
 integer, dimension(:,:), allocatable, target, private :: &
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
   FE_mapElemtype, &
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


 type, public, extends(tMesh) :: tMesh_abaqus
 
 integer:: &
   mesh_Nelems, &                                                                                   !< total number of elements in mesh (including non-DAMASK elements)
   mesh_maxNnodes, &                                                                                !< max number of nodes in any CP element
   mesh_NelemSets, &
   mesh_maxNelemInSet, &
   mesh_Nmaterials
 character(len=64), dimension(:), allocatable :: &
   mesh_nameElemSet, &                                                                              !< names of elementSet
   mesh_nameMaterial, &                                                                             !< names of material in solid section
   mesh_mapMaterial                                                                                 !< name of elementSet for material
 integer, dimension(:,:), allocatable :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
 logical:: noPart                                                                         !< for cases where the ABAQUS input file does not use part/assembly information

 contains 
   procedure, pass(self) :: tMesh_abaqus_init
   generic, public :: init => tMesh_abaqus_init
 end type tMesh_abaqus
 
 type(tMesh_abaqus), public, protected :: theMesh
 
contains

subroutine tMesh_abaqus_init(self,elemType,nodes)
 
 
 class(tMesh_abaqus) :: self
 real(pReal), dimension(:,:), intent(in) :: nodes
 integer, intent(in) :: elemType
 
 call self%tMesh%init('mesh',elemType,nodes)
 
end subroutine tMesh_abaqus_init

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)
 use DAMASK_interface
 use IO, only: &
   IO_open_InputFile, &
   IO_error
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

 
 integer, parameter :: FILEUNIT = 222
 integer, intent(in), optional :: el, ip
 integer :: j
 logical :: myDebug

 write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'

 mesh_unitlength = numerics_unitlength                                                              ! set physical extent of a length unit in mesh

 myDebug = (iand(debug_level(debug_mesh),debug_levelBasic) /= 0)

 call IO_open_inputFile(FILEUNIT)                                                                   ! parse info from input file...
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
 call mesh_get_damaskOptions(mesh_periodicSurface,FILEUNIT)
 if (myDebug) write(6,'(a)') ' Got DAMASK options'; flush(6)
 close (FILEUNIT)
  
 call theMesh%init(mesh_element(2,1),mesh_node0)
 call theMesh%setNelems(mesh_NcpElems)
 call mesh_build_FEdata                                                                             ! get properties of the different types of elements
 
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
 call mesh_build_nodeTwins
 if (myDebug) write(6,'(a)') ' Built node twins'; flush(6)
 call mesh_build_sharedElems
 if (myDebug) write(6,'(a)') ' Built shared elements'; flush(6)
 call mesh_build_ipNeighborhood
 if (myDebug) write(6,'(a)') ' Built IP neighborhood'; flush(6)

 if (usePingPong .and. (mesh_Nelems /= mesh_NcpElems)) &
   call IO_error(600)                                                                          ! ping-pong must be disabled when having non-DAMASK elements
 if (debug_e < 1 .or. debug_e > mesh_NcpElems) &
   call IO_error(602,ext_msg='element')                                                        ! selected element does not exist
 if (debug_i < 1 .or. debug_i > FE_Nips(FE_geomtype(mesh_element(2,debug_e)))) &
   call IO_error(602,ext_msg='IP')                                                             ! selected element does not have requested IP
 FEsolving_execElem = [ 1,mesh_NcpElems ]                                                      ! parallel loop bounds set to comprise all DAMASK elements
 allocate(FEsolving_execIP(2,mesh_NcpElems), source=1)                                    ! parallel loop bounds set to comprise from first IP...
 forall (j = 1:mesh_NcpElems) FEsolving_execIP(2,j) = FE_Nips(FE_geomtype(mesh_element(2,j)))  ! ...up to own IP count for each element
 allocate(calcMode(mesh_maxNips,mesh_NcpElems))
 calcMode = .false.                                                                                 ! pretend to have collected what first call is asking (F = I)
 calcMode(ip,mesh_FEasCP('elem',el)) = .true.                                                       ! first ip,el needs to be already pingponged to "calc"


! better name
 theMesh%homogenizationAt  = mesh_element(3,:)
 theMesh%microstructureAt  = mesh_element(4,:)

   call discretization_init(mesh_element(3,:),mesh_element(4,:),&
                            reshape(mesh_ipCoordinates,[3,theMesh%elem%nIPs*theMesh%nElems]),&
                            mesh_node0)
 call geometry_plastic_nonlocal_setIPvolume(mesh_ipVolume)
 call geometry_plastic_nonlocal_setIPneighborhood(mesh_ipNeighborhood)
 call geometry_plastic_nonlocal_setIParea(mesh_IParea)
 call geometry_plastic_nonlocal_setIPareaNormal(mesh_IPareaNormal)
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief check if the input file for Abaqus contains part info
!--------------------------------------------------------------------------------------------------
logical function hasNoPart(fileUnit)
 use IO, only: &
   IO_stringPos, &
   IO_stringValue, &
   IO_lc
 
 
 integer,    intent(in)                :: fileUnit

 integer, allocatable, dimension(:)    :: chunkPos
 character(len=65536)                        :: line

 hasNoPart = .true.

 rewind(fileUnit)
 do
   read(fileUnit,'(a65536)',END=620) line
   chunkPos = IO_stringPos(line)
   if (IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) then
     hasNoPart = .false.
     exit
   endif
 enddo

620 end function hasNoPart

end subroutine mesh_init








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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart

 mesh_Nnodes = 0
 mesh_Nelems = 0
 
 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if (inPart .or. noPart) then
     select case ( IO_lc(IO_stringValue(line,chunkPos,1)))
       case('*node')
          if( &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'print'    .and. &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'file'     .and. &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'response' &
             ) &
            mesh_Nnodes = mesh_Nnodes + IO_countDataLines(fileUnit)
       case('*element')
          if( &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'matrix'   .and. &
              IO_lc(IO_stringValue(line,chunkPos,2)) /= 'response' &
             ) then
            mesh_Nelems = mesh_Nelems + IO_countDataLines(fileUnit)
          endif
     endselect
   endif
 enddo

 if (mesh_Nnodes < 2)  call IO_error(error_ID=900)
 if (mesh_Nelems == 0) call IO_error(error_ID=901)

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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart

 mesh_NelemSets     = 0
 mesh_maxNelemInSet = mesh_Nelems                                                                   ! have to be conservative, since Abaqus allows for recursive definitons

 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,chunkPos,1)) == '*elset' ) &
     mesh_NelemSets = mesh_NelemSets + 1
 enddo

 if (mesh_NelemSets == 0) call IO_error(error_ID=902)

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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart

 mesh_Nmaterials = 0

 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,chunkPos,1)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,chunkPos,2)) == 'section' ) &
     mesh_Nmaterials = mesh_Nmaterials + 1
 enddo

 if (mesh_Nmaterials == 0) call IO_error(error_ID=903)

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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart
 integer :: elemSet,i

 allocate (mesh_nameElemSet(mesh_NelemSets)); mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1+mesh_maxNelemInSet,mesh_NelemSets),source=0)


 elemSet = 0
 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,chunkPos,1)) == '*elset' ) then
     elemSet = elemSet + 1
     mesh_nameElemSet(elemSet)  = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2)),'elset'))
     mesh_mapElemSet(:,elemSet) = IO_continuousIntValues(fileUnit,mesh_Nelems,mesh_nameElemSet,&
                                          mesh_mapElemSet,elemSet-1)
   endif
 enddo

 do i = 1,elemSet
   if (mesh_mapElemSet(1,i) == 0) call IO_error(error_ID=904,ext_msg=mesh_nameElemSet(i))
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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart
 integer :: i,c
 character(len=64) :: elemSetName,materialName

 allocate (mesh_nameMaterial(mesh_Nmaterials)); mesh_nameMaterial = ''
 allocate (mesh_mapMaterial(mesh_Nmaterials));  mesh_mapMaterial = ''

 c = 0
 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,chunkPos,1)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,chunkPos,2)) == 'section' ) then

     elemSetName = ''
     materialName = ''

     do i = 3,chunkPos(1)
       if (IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'elset') /= '') &
         elemSetName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'elset'))
       if (IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'material') /= '') &
         materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,i)),'material'))
     enddo

     if (elemSetName /= '' .and. materialName /= '') then
       c = c + 1
       mesh_nameMaterial(c) = materialName                                                          ! name of material used for this section
       mesh_mapMaterial(c)  = elemSetName                                                           ! mapped to respective element set
     endif
   endif
 enddo

 if (c==0) call IO_error(error_ID=905)
 do i=1,c
   if (mesh_nameMaterial(i)=='' .or. mesh_mapMaterial(i)=='') call IO_error(error_ID=905)
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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: materialFound
 integer :: i,k
 character(len=64) ::materialName,elemSetName

 mesh_NcpElems = 0
 materialFound = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   select case ( IO_lc(IO_stringValue(line,chunkPos,1)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if (IO_lc(IO_StringValue(line,chunkPos,2)) == 'material' .and. materialFound) then
         do i = 1,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) &                                            ! matched?
                 mesh_NcpElems = mesh_NcpElems + mesh_mapElemSet(1,k)                               ! add those elem count
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo

 if (mesh_NcpElems == 0) call IO_error(error_ID=906)

end subroutine mesh_abaqus_count_cpElements


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPelem'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_elements(fileUnit)

 use math, only: math_sort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: materialFound
 integer ::i,j,k,cpElem
 character (len=64) materialName,elemSetName                                                        ! why limited to 64? ABAQUS?

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems), source = 0)

 cpElem = 0
 materialFound = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   select case ( IO_lc(IO_stringValue(line,chunkPos,1)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if (IO_lc(IO_stringValue(line,chunkPos,2)) == 'material' .and. materialFound) then
         do i = 1,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) then                                         ! matched?
                 do j = 1,mesh_mapElemSet(1,k)
                   cpElem = cpElem + 1
                   mesh_mapFEtoCPelem(1,cpElem) = mesh_mapElemSet(1+j,k)                       ! store FE id
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

 call math_sort(mesh_mapFEtoCPelem,1,int(size(mesh_mapFEtoCPelem,2),pInt))               ! should be mesh_NcpElems

 if (int(size(mesh_mapFEtoCPelem),pInt) < 2) call IO_error(error_ID=907)

end subroutine mesh_abaqus_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!! Allocates global array 'mesh_mapFEtoCPnode'
!--------------------------------------------------------------------------------------------------
subroutine mesh_abaqus_map_nodes(fileUnit)

 use math, only: math_sort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_intValue, &
                 IO_error

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart
 integer :: i,c,cpNode

 allocate (mesh_mapFEtoCPnode(2,mesh_Nnodes), source=0)
 
 cpNode = 0
 inPart = .false. 
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'response' ) &
   ) then
     c = IO_countDataLines(fileUnit)
     do i = 1,c
       backspace(fileUnit)
     enddo
     do i = 1,c
       read (fileUnit,'(a300)') line
       chunkPos = IO_stringPos(line)
       cpNode = cpNode + 1
       mesh_mapFEtoCPnode(1,cpNode) = IO_intValue(line,chunkPos,1)
       mesh_mapFEtoCPnode(2,cpNode) = cpNode
     enddo
   endif
 enddo

 call math_sort(mesh_mapFEtoCPnode,1,int(size(mesh_mapFEtoCPnode,2),pInt))

 if (int(size(mesh_mapFEtoCPnode),pInt) == 0) call IO_error(error_ID=908)

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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart
 integer :: i,j,m,c

 allocate ( mesh_node0 (3,mesh_Nnodes), source=0.0_pReal)
 allocate ( mesh_node  (3,mesh_Nnodes), source=0.0_pReal)

 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'response' ) &
   ) then
     c = IO_countDataLines(fileUnit)                                                                  ! how many nodes are defined here?
     do i = 1,c
       backspace(fileUnit)                                                                            ! rewind to first entry
     enddo
     do i = 1,c
       read (fileUnit,'(a300)') line
       chunkPos = IO_stringPos(line)
       m = mesh_FEasCP('node',IO_intValue(line,chunkPos,1))
       do j=1, 3
         mesh_node0(j,m) = mesh_unitlength * IO_floatValue(line,chunkPos,j+1)
       enddo
     enddo
   endif
 enddo

 if (int(size(mesh_node0,2),pInt) /= mesh_Nnodes) call IO_error(error_ID=909)
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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart
 integer :: i,c,t,g

 mesh_maxNnodes       = 0
 mesh_maxNips         = 0
 mesh_maxNipNeighbors = 0
 mesh_maxNcellnodes   = 0


 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2)),'type'))           ! remember elem type
     g = FE_geomtype(t)
     c = FE_celltype(g)
     mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(t))
     mesh_maxNips =         max(mesh_maxNips,FE_Nips(g))
     mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(c))
     mesh_maxNcellnodes =   max(mesh_maxNcellnodes,FE_Ncellnodes(g))
   endif
 enddo

end subroutine mesh_abaqus_count_cpSizes


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

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 logical :: inPart, materialFound
 integer :: i,j,k,c,e,t,homog,micro, nNodesAlreadyRead
 character (len=64) :: materialName,elemSetName

 allocate(mesh_element (4+mesh_maxNnodes,mesh_NcpElems), source=0)
 mesh_elemType = -1

 inPart = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,chunkPos,2)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,chunkPos,1)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,chunkPos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,chunkPos,2)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,chunkPos,2)),'type'))        ! remember elem type
     c = IO_countDataLines(fileUnit)
     do i = 1,c
       backspace(fileUnit)
     enddo
     do i = 1,c
       read (fileUnit,'(a300)') line
       chunkPos = IO_stringPos(line)                                                       ! limit to 64 nodes max
       e = mesh_FEasCP('elem',IO_intValue(line,chunkPos,1))
       if (e /= 0) then                                                                        ! disregard non CP elems
         mesh_element(1,e) = -1                                                                ! DEPRECATED
         if (mesh_elemType /= t .and. mesh_elemType /= -1) &
           call IO_error(191,el=t,ip=mesh_elemType)
         mesh_elemType = t
         mesh_element(2,e) = t                                                                     ! elem type
         nNodesAlreadyRead = 0
         do j = 1,chunkPos(1)-1
           mesh_element(4+j,e) = mesh_FEasCP('node',IO_intValue(line,chunkPos,1+j))         ! put CP ids of nodes to position 5:
         enddo
         nNodesAlreadyRead = chunkPos(1) - 1
         do while(nNodesAlreadyRead < FE_Nnodes(t))                                                ! read on if not all nodes in one line
           read (fileUnit,'(a300)') line
           chunkPos = IO_stringPos(line)
           do j = 1,chunkPos(1)
             mesh_element(4+nNodesAlreadyRead+j,e) &
               = mesh_FEasCP('node',IO_IntValue(line,chunkPos,j))                                     ! CP ids of nodes
           enddo
           nNodesAlreadyRead = nNodesAlreadyRead + chunkPos(1)
         enddo
       endif
     enddo
   endif
 enddo


 rewind(fileUnit)                                                                                 ! just in case "*material" definitions apear before "*element"

 materialFound = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   select case ( IO_lc(IO_StringValue(line,chunkPos,1)))
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_StringValue(line,chunkPos,2)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if ( IO_lc(IO_StringValue(line,chunkPos,2)) == 'material' .and. &
            materialFound ) then
         read (fileUnit,'(a300)') line                                                              ! read homogenization and microstructure
         chunkPos = IO_stringPos(line)
         homog = nint(IO_floatValue(line,chunkPos,1),pInt)
         micro = nint(IO_floatValue(line,chunkPos,2),pInt)
         do i = 1,mesh_Nmaterials                                                              ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                                           ! found one
             elemSetName = mesh_mapMaterial(i)                                                      ! take corresponding elemSet
             do k = 1,mesh_NelemSets                                                           ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) then                                         ! matched?
                 do j = 1,mesh_mapElemSet(1,k)
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

end subroutine mesh_abaqus_build_elements


!--------------------------------------------------------------------------------------------------
!> @brief get any additional damask options from input file, sets mesh_periodicSurface
!--------------------------------------------------------------------------------------------------
subroutine mesh_get_damaskOptions(periodic_surface,fileUnit)

use IO, only: &
  IO_lc, &
  IO_stringValue, &
  IO_stringPos

 
 integer, intent(in) :: fileUnit

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 integer :: chunk, Nchunks
 character(len=300) ::  v
 logical, dimension(3) :: periodic_surface
 

 periodic_surface = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   Nchunks = chunkPos(1)
   if (IO_lc(IO_stringValue(line,chunkPos,1)) == '**damask' .and. Nchunks > 1) then          ! found keyword for damask option and there is at least one more chunk to read
     select case(IO_lc(IO_stringValue(line,chunkPos,2)))
       case('periodic')                                                                             ! damask Option that allows to specify periodic fluxes
         do chunk = 3,Nchunks                                                                  ! loop through chunks (skipping the keyword)
            v = IO_lc(IO_stringValue(line,chunkPos,chunk))                                          ! chunk matches keyvalues x,y, or z?
            mesh_periodicSurface(1) = mesh_periodicSurface(1) .or. v == 'x'
            mesh_periodicSurface(2) = mesh_periodicSurface(2) .or. v == 'y'
            mesh_periodicSurface(3) = mesh_periodicSurface(3) .or. v == 'z'
         enddo
     endselect
   endif
 enddo

end subroutine mesh_get_damaskOptions


!--------------------------------------------------------------------------------------------------
!> @brief Split CP elements into cells.
!> @details Build a mapping between cells and the corresponding cell nodes ('mesh_cell').
!> Cell nodes that are also matching nodes are unique in the list of cell nodes,
!> all others (currently) might be stored more than once.
!> Also allocates the 'mesh_node' array.
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_cellconnectivity

 
 integer, dimension(:), allocatable :: &
   matchingNode2cellnode
 integer, dimension(:,:), allocatable :: &
   cellnodeParent
 integer, dimension(mesh_maxNcellnodes) :: &
   localCellnode2globalCellnode
 integer :: &
   e,t,g,c,n,i, &
   matchingNodeID, &
   localCellnodeID

 allocate(mesh_cell(FE_maxNcellnodesPerCell,mesh_maxNips,mesh_NcpElems), source=0)
 allocate(matchingNode2cellnode(mesh_Nnodes),                            source=0)
 allocate(cellnodeParent(2,mesh_maxNcellnodes*mesh_NcpElems),       source=0)

!--------------------------------------------------------------------------------------------------
! Count cell nodes (including duplicates) and generate cell connectivity list
 mesh_Ncellnodes = 0
 mesh_Ncells = 0
 do e = 1,mesh_NcpElems                                                                        ! loop over cpElems
   t = mesh_element(2,e)                                                                       ! get element type
   g = FE_geomtype(t)                                                                               ! get geometry type
   c = FE_celltype(g)                                                                               ! get cell type
   localCellnode2globalCellnode = 0
   mesh_Ncells = mesh_Ncells + FE_Nips(g)
   do i = 1,FE_Nips(g)                                                                         ! loop over ips=cells in this element
     do n = 1,FE_NcellnodesPerCell(c)                                                          ! loop over cell nodes in this cell
       localCellnodeID = FE_cell(n,i,g)
       if (localCellnodeID <= FE_NmatchingNodes(g)) then                                            ! this cell node is a matching node
         matchingNodeID = mesh_element(4+localCellnodeID,e)
         if (matchingNode2cellnode(matchingNodeID) == 0) then                                  ! if this matching node does not yet exist in the glbal cell node list ...
           mesh_Ncellnodes = mesh_Ncellnodes + 1                                               ! ... count it as cell node ...
           matchingNode2cellnode(matchingNodeID) = mesh_Ncellnodes                                  ! ... and remember its global ID
           cellnodeParent(1,mesh_Ncellnodes) = e                                               ! ... and where it belongs to
           cellnodeParent(2,mesh_Ncellnodes) = localCellnodeID
         endif
         mesh_cell(n,i,e) = matchingNode2cellnode(matchingNodeID)
       else                                                                                         ! this cell node is no matching node
         if (localCellnode2globalCellnode(localCellnodeID) == 0) then                          ! if this local cell node does not yet exist in the  global cell node list ...
           mesh_Ncellnodes = mesh_Ncellnodes + 1                                               ! ... count it as cell node ...
           localCellnode2globalCellnode(localCellnodeID) = mesh_Ncellnodes                          ! ... and remember its global ID ...
           cellnodeParent(1,mesh_Ncellnodes) = e                                               ! ... and it belongs to
           cellnodeParent(2,mesh_Ncellnodes) = localCellnodeID
         endif
         mesh_cell(n,i,e) = localCellnode2globalCellnode(localCellnodeID)
       endif
     enddo
   enddo
 enddo

 allocate(mesh_cellnodeParent(2,mesh_Ncellnodes))
 allocate(mesh_cellnode(3,mesh_Ncellnodes))
 forall(n = 1:mesh_Ncellnodes)
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

 
 integer,                         intent(in) :: Ncellnodes                                    !< requested number of cellnodes
 real(pReal), dimension(3,mesh_Nnodes), intent(in) :: nodes
 real(pReal), dimension(3,Ncellnodes) :: mesh_build_cellnodes

 integer :: &
   e,t,n,m, &
   localCellnodeID
 real(pReal), dimension(3) :: &
   myCoords

 mesh_build_cellnodes = 0.0_pReal
!$OMP PARALLEL DO PRIVATE(e,localCellnodeID,t,myCoords)
 do n = 1,Ncellnodes                                                                           ! loop over cell nodes
   e = mesh_cellnodeParent(1,n)
   localCellnodeID = mesh_cellnodeParent(2,n)
   t = mesh_element(2,e)                                                                            ! get element type
   myCoords = 0.0_pReal
   do m = 1,FE_Nnodes(t)
     myCoords = myCoords + nodes(1:3,mesh_element(4+m,e)) &
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

 
 integer ::                                e,t,g,c,i,m,f,n
 real(pReal), dimension(FE_maxNcellnodesPerCellface,FE_maxNcellfaces) :: subvolume

 allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems),source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,m,subvolume)
   do e = 1,mesh_NcpElems                                                                      ! loop over cpElems
     t = mesh_element(2,e)                                                                     ! get element type
     g = FE_geomtype(t)                                                                             ! get geometry type
     c = FE_celltype(g)                                                                             ! get cell type
     select case (c)

       case (1)                                                                                ! 2D 3node
         forall (i = 1:FE_Nips(g)) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_areaTriangle(mesh_cellnode(1:3,mesh_cell(1,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(3,i,e)))

       case (2)                                                                                ! 2D 4node
         forall (i = 1:FE_Nips(g)) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_areaTriangle(mesh_cellnode(1:3,mesh_cell(1,i,e)), &            ! here we assume a planar shape, so division in two triangles suffices
                                                  mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(3,i,e))) &
                              + math_areaTriangle(mesh_cellnode(1:3,mesh_cell(3,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(4,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(1,i,e)))

       case (3)                                                                                ! 3D 4node
         forall (i = 1:FE_Nips(g)) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_volTetrahedron(mesh_cellnode(1:3,mesh_cell(1,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(3,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(4,i,e)))

       case (4)                                                                                ! 3D 8node
         m = FE_NcellnodesPerCellface(c)
         do i = 1,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           subvolume = 0.0_pReal
           forall(f = 1:FE_NipNeighbors(c), n = 1:FE_NcellnodesPerCellface(c)) &
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

 
 integer :: e,t,g,c,i,n
 real(pReal), dimension(3) :: myCoords

 if (.not. allocated(mesh_ipCoordinates)) &
   allocate(mesh_ipCoordinates(3,mesh_maxNips,mesh_NcpElems),source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,myCoords)
 do e = 1,mesh_NcpElems                                                                        ! loop over cpElems
   t = mesh_element(2,e)                                                                       ! get element type
   g = FE_geomtype(t)                                                                               ! get geometry type
   c = FE_celltype(g)                                                                               ! get cell type
   do i = 1,FE_Nips(g)                                                                         ! loop over ips=cells in this element
     myCoords = 0.0_pReal
     do n = 1,FE_NcellnodesPerCell(c)                                                          ! loop over cell nodes in this cell
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

 
 integer, intent(in) :: el, &                                                                  !< element number
                              ip                                                                     !< integration point number
 real(pReal), dimension(3) :: mesh_cellCenterCoordinates                                             !< x,y,z coordinates of the cell center of the requested IP cell
 integer :: t,g,c,n

 t = mesh_element(2,el)                                                                         ! get element type
 g = FE_geomtype(t)                                                                                  ! get geometry type
 c = FE_celltype(g)                                                                                  ! get cell type
 mesh_cellCenterCoordinates = 0.0_pReal
 do n = 1,FE_NcellnodesPerCell(c)                                                               ! loop over cell nodes in this cell
   mesh_cellCenterCoordinates = mesh_cellCenterCoordinates + mesh_cellnode(1:3,mesh_cell(n,ip,el))
 enddo
 mesh_cellCenterCoordinates = mesh_cellCenterCoordinates / real(FE_NcellnodesPerCell(c),pReal)

 end function mesh_cellCenterCoordinates






!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipAreas
 use math, only: &
   math_cross

 
 integer :: e,t,g,c,i,f,n,m
 real(pReal), dimension (3,FE_maxNcellnodesPerCellface) :: nodePos, normals
 real(pReal), dimension(3) :: normal

 allocate(mesh_ipArea(mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems), source=0.0_pReal)
 allocate(mesh_ipAreaNormal(3,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems), source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,nodePos,normal,normals)
   do e = 1,mesh_NcpElems                                                                      ! loop over cpElems
     t = mesh_element(2,e)                                                                     ! get element type
     g = FE_geomtype(t)                                                                             ! get geometry type
     c = FE_celltype(g)                                                                             ! get cell type
     select case (c)

       case (1,2)                                                                         ! 2D 3 or 4 node
         do i = 1,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           do f = 1,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal(1) =   nodePos(2,2) - nodePos(2,1)                                              ! x_normal =  y_connectingVector
             normal(2) = -(nodePos(1,2) - nodePos(1,1))                                             ! y_normal = -x_connectingVector
             normal(3) = 0.0_pReal
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (3)                                                                                ! 3D 4node
         do i = 1,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           do f = 1,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal = math_cross(nodePos(1:3,2) - nodePos(1:3,1), &
                                         nodePos(1:3,3) - nodePos(1:3,1))
             mesh_ipArea(f,i,e) = norm2(normal)
             mesh_ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (4)                                                                                ! 3D 8node
         ! for this cell type we get the normal of the quadrilateral face as an average of
         ! four normals of triangular subfaces; since the face consists only of two triangles,
         ! the sum has to be divided by two; this whole prcedure tries to compensate for
         ! probable non-planar cell surfaces
         m = FE_NcellnodesPerCellface(c)
         do i = 1,FE_Nips(g)                                                                   ! loop over ips=cells in this element
           do f = 1,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             forall(n = 1:FE_NcellnodesPerCellface(c)) &
               normals(1:3,n) = 0.5_pReal &
                              * math_cross(nodePos(1:3,1+mod(n  ,m)) - nodePos(1:3,n), &
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

 
 integer dir, &      ! direction of periodicity
               node, &
               minimumNode, &
               maximumNode, &
               n1, &
               n2
 integer, dimension(mesh_Nnodes+1) :: minimumNodes, maximumNodes                              ! list of surface nodes (minimum and maximum coordinate value) with first entry giving the number of nodes
 real(pReal)   minCoord, maxCoord, &                                                                ! extreme positions in one dimension
               tolerance                                                                            ! tolerance below which positions are assumed identical
 real(pReal), dimension(3) ::  distance                                                             ! distance between two nodes in all three coordinates
 logical, dimension(mesh_Nnodes) :: unpaired

 allocate(mesh_nodeTwins(3,mesh_Nnodes))
 mesh_nodeTwins = 0

 tolerance = 0.001_pReal * minval(mesh_ipVolume) ** 0.333_pReal

 do dir = 1,3                                                                             ! check periodicity in directions of x,y,z
   if (mesh_periodicSurface(dir)) then                                                              ! only if periodicity is requested


     !*** find out which nodes sit on the surface
     !*** and have a minimum or maximum position in this dimension

     minimumNodes = 0
     maximumNodes = 0
     minCoord = minval(mesh_node0(dir,:))
     maxCoord = maxval(mesh_node0(dir,:))
     do node = 1,mesh_Nnodes                                                                   ! loop through all nodes and find surface nodes
       if (abs(mesh_node0(dir,node) - minCoord) <= tolerance) then
         minimumNodes(1) = minimumNodes(1) + 1
         minimumNodes(minimumNodes(1)+1) = node
       elseif (abs(mesh_node0(dir,node) - maxCoord) <= tolerance) then
         maximumNodes(1) = maximumNodes(1) + 1
         maximumNodes(maximumNodes(1)+1) = node
       endif
     enddo


     !*** find the corresponding node on the other side with the same position in this dimension

     unpaired = .true.
     do n1 = 1,minimumNodes(1)
       minimumNode = minimumNodes(n1+1)
       if (unpaired(minimumNode)) then
         do n2 = 1,maximumNodes(1)
           maximumNode = maximumNodes(n2+1)
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

 
 integer(pint)   e, &                                                                                ! element index
                 g, &                                                                                ! element type
                 node, &                                                                             ! CP node index
                 n, &                                                                                ! node index per element
                 myDim, &                                                                            ! dimension index
                 nodeTwin                                                                            ! node twin in the specified dimension
 integer, dimension (mesh_Nnodes) :: node_count
 integer, dimension(:), allocatable :: node_seen

 allocate(node_seen(maxval(FE_NmatchingNodes)))

 node_count = 0

 do e = 1,mesh_NcpElems
   g = FE_geomtype(mesh_element(2,e))                                                                ! get elemGeomType
   node_seen = 0                                                                                ! reset node duplicates
   do n = 1,FE_NmatchingNodes(g)                                                                ! check each node of element
     node = mesh_element(4+n,e)
     if (all(node_seen /= node)) then
       node_count(node) = node_count(node) + 1                                                  ! if FE node not yet encountered -> count it
       do myDim = 1,3                                                                      ! check in each dimension...
         nodeTwin = mesh_nodeTwins(myDim,node)
         if (nodeTwin > 0) &                                                                    ! if I am a twin of some node...
           node_count(nodeTwin) = node_count(nodeTwin) + 1                                      ! -> count me again for the twin node
       enddo
     endif
     node_seen(n) = node                                                                             ! remember this node to be counted already
   enddo
 enddo

 mesh_maxNsharedElems = int(maxval(node_count),pInt)                                                 ! most shared node

 allocate(mesh_sharedElem(1+mesh_maxNsharedElems,mesh_Nnodes),source=0)

 do e = 1,mesh_NcpElems
   g = FE_geomtype(mesh_element(2,e))                                                                ! get elemGeomType
   node_seen = 0
   do n = 1,FE_NmatchingNodes(g)
     node = mesh_element(4+n,e)
     if (all(node_seen /= node)) then
       mesh_sharedElem(1,node) = mesh_sharedElem(1,node) + 1                                    ! count for each node the connected elements
       mesh_sharedElem(mesh_sharedElem(1,node)+1,node) = e                                      ! store the respective element id
       do myDim = 1,3                                                                      ! check in each dimension...
         nodeTwin = mesh_nodeTwins(myDim,node)
         if (nodeTwin > 0) then                                                                 ! if i am a twin of some node...
           mesh_sharedElem(1,nodeTwin) = mesh_sharedElem(1,nodeTwin) + 1                        ! ...count me again for the twin
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
 
 integer      ::           myElem, &                                                           ! my CP element index
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
 integer, dimension(FE_maxmaxNnodesAtIP) :: &
                                 linkedNodes = 0, &
                                 matchingNodes
 logical checkTwins

 allocate(mesh_ipNeighborhood(3,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems))
 mesh_ipNeighborhood = 0


 do myElem = 1,mesh_NcpElems                                                                    ! loop over cpElems
   myType = FE_geomtype(mesh_element(2,myElem))                                                      ! get elemGeomType
   do myIP = 1,FE_Nips(myType)                                                                  ! loop over IPs of elem

     do neighbor = 1,FE_NipNeighbors(FE_celltype(myType))                                       ! loop over neighbors of IP
       neighboringIPkey = FE_ipNeighbor(neighbor,myIP,myType)

       !*** if the key is positive, the neighbor is inside the element
       !*** that means, we have already found our neighboring IP

       if (neighboringIPkey > 0) then
         mesh_ipNeighborhood(1,neighbor,myIP,myElem) = myElem
         mesh_ipNeighborhood(2,neighbor,myIP,myElem) = neighboringIPkey


       !*** if the key is negative, the neighbor resides in a neighboring element
       !*** that means, we have to look through the face indicated by the key and see which element is behind that face

       elseif (neighboringIPkey < 0) then                                                       ! neighboring element's IP
         myFace = -neighboringIPkey
         call mesh_faceMatch(myElem, myFace, matchingElem, matchingFace)                             ! get face and CP elem id of face match
         if (matchingElem > 0) then                                                             ! found match?
           neighboringType = FE_geomtype(mesh_element(2,matchingElem))

           !*** trivial solution if neighbor has only one IP

           if (FE_Nips(neighboringType) == 1) then
             mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
             mesh_ipNeighborhood(2,neighbor,myIP,myElem) = 1
             cycle
           endif

           !*** find those nodes which build the link to the neighbor

           NlinkedNodes = 0
           linkedNodes = 0
           do a = 1,FE_maxNnodesAtIP(myType)                                                    ! figure my anchor nodes on connecting face
             anchor = FE_nodesAtIP(a,myIP,myType)
             if (anchor /= 0) then                                                              ! valid anchor node
               if (any(FE_face(:,myFace,myType) == anchor)) then                                     ! ip anchor sits on face?
                 NlinkedNodes = NlinkedNodes + 1
                 linkedNodes(NlinkedNodes) = mesh_element(4+anchor,myElem)                      ! CP id of anchor node
               else                                                                                  ! something went wrong with the linkage, since not all anchors sit on my face
                 NlinkedNodes = 0
                 linkedNodes = 0
                 exit
               endif
             endif
           enddo

           !*** loop through the ips of my neighbor
           !*** and try to find an ip with matching nodes
           !*** also try to match with node twins

 checkCandidateIP: do candidateIP = 1,FE_Nips(neighboringType)
             NmatchingNodes = 0
             matchingNodes = 0
             do a = 1,FE_maxNnodesAtIP(neighboringType)                                         ! check each anchor node of that ip
               anchor = FE_nodesAtIP(a,candidateIP,neighboringType)
               if (anchor /= 0) then                                                            ! valid anchor node
                 if (any(FE_face(:,matchingFace,neighboringType) == anchor)) then                    ! sits on matching face?
                   NmatchingNodes = NmatchingNodes + 1
                   matchingNodes(NmatchingNodes) = mesh_element(4+anchor,matchingElem)               ! CP id of neighbor's anchor node
                 else                                                                                ! no matching, because not all nodes sit on the matching face
                   NmatchingNodes = 0
                   matchingNodes = 0
                   exit
                 endif
               endif
             enddo

             if (NmatchingNodes /= NlinkedNodes) &                                                   ! this ip has wrong count of anchors on face
               cycle checkCandidateIP

             !*** check "normal" nodes whether they match or not

             checkTwins = .false.
             do a = 1,NlinkedNodes
               if (all(matchingNodes /= linkedNodes(a))) then                                        ! this linkedNode does not match any matchingNode
                 checkTwins = .true.
                 exit                                                                                ! no need to search further
               endif
             enddo

             !*** if no match found, then also check node twins

             if(checkTwins) then
               dir = int(maxloc(abs(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem)),1),pInt)            ! check for twins only in direction of the surface normal
               do a = 1,NlinkedNodes
                 twin_of_linkedNode = mesh_nodeTwins(dir,linkedNodes(a))
                 if (twin_of_linkedNode == 0 .or. &                                             ! twin of linkedNode does not exist...
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
 do myElem = 1,mesh_NcpElems                                                                    ! loop over cpElems
   myType = FE_geomtype(mesh_element(2,myElem))                                                      ! get elemGeomType
   do myIP = 1,FE_Nips(myType)                                                                  ! loop over IPs of elem
     do neighbor = 1,FE_NipNeighbors(FE_celltype(myType))                                       ! loop over neighbors of IP
       neighboringElem = mesh_ipNeighborhood(1,neighbor,myIP,myElem)
       neighboringIP   = mesh_ipNeighborhood(2,neighbor,myIP,myElem)
       if (neighboringElem > 0 .and. neighboringIP > 0) then                               ! if neighbor exists ...
         neighboringType = FE_geomtype(mesh_element(2,neighboringElem))
         do pointingToMe = 1,FE_NipNeighbors(FE_celltype(neighboringType))                      ! find neighboring index that points from my neighbor to myself
           if (    myElem == mesh_ipNeighborhood(1,pointingToMe,neighboringIP,neighboringElem) &
               .and. myIP == mesh_ipNeighborhood(2,pointingToMe,neighboringIP,neighboringElem)) then ! possible candidate
             if (math_inner(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem),&
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
 
  call geometry_plastic_nonlocal_set_IPneighborhood(mesh_ipNeighborhood)
 
 contains
 !--------------------------------------------------------------------------------------------------
!> @brief find face-matching element of same type
!--------------------------------------------------------------------------------------------------
subroutine mesh_faceMatch(elem, face ,matchingElem, matchingFace)


integer, intent(out) ::     matchingElem, &                                                   ! matching CP element ID
                                  matchingFace                                                      ! matching face ID
integer, intent(in) ::      face, &                                                           ! face ID
                                  elem                                                              ! CP elem ID
integer, dimension(FE_NmatchingNodesPerFace(face,FE_geomtype(mesh_element(2,elem)))) :: &
                                  myFaceNodes                                                       ! global node ids on my face
integer        ::           myType, &
                                  candidateType, &
                                  candidateElem, &
                                  candidateFace, &
                                  candidateFaceNode, &
                                  minNsharedElems, &
                                  NsharedElems, &
                                  lonelyNode = 0, &
                                  i, &
                                  n, &
                                  dir                                                               ! periodicity direction
integer, dimension(:), allocatable :: element_seen
logical checkTwins

matchingElem = 0
matchingFace = 0
minNsharedElems = mesh_maxNsharedElems + 1                                                     ! init to worst case
myType = FE_geomtype(mesh_element(2,elem))                                                     ! figure elemGeomType

do n = 1,FE_NmatchingNodesPerFace(face,myType)                                                 ! loop over nodes on face
  myFaceNodes(n) = mesh_element(4+FE_face(n,face,myType),elem)                                 ! CP id of face node
  NsharedElems = mesh_sharedElem(1,myFaceNodes(n))                                             ! figure # shared elements for this node
  if (NsharedElems < minNsharedElems) then
    minNsharedElems = NsharedElems                                                                  ! remember min # shared elems
    lonelyNode = n                                                                                  ! remember most lonely node
  endif
enddo

allocate(element_seen(minNsharedElems))
element_seen = 0

checkCandidate: do i = 1,minNsharedElems                                                       ! iterate over lonelyNode's shared elements
  candidateElem = mesh_sharedElem(1+i,myFaceNodes(lonelyNode))                                 ! present candidate elem
  if (all(element_seen /= candidateElem)) then                                                      ! element seen for the first time?
    element_seen(i) = candidateElem
    candidateType = FE_geomtype(mesh_element(2,candidateElem))                                 ! figure elemGeomType of candidate
checkCandidateFace: do candidateFace = 1,FE_maxNipNeighbors                                    ! check each face of candidate
      if (FE_NmatchingNodesPerFace(candidateFace,candidateType) &
          /= FE_NmatchingNodesPerFace(face,myType) &                                                ! incompatible face
          .or. (candidateElem == elem .and. candidateFace == face)) then                            ! this is my face
        cycle checkCandidateFace
      endif
      checkTwins = .false.
      do n = 1,FE_NmatchingNodesPerFace(candidateFace,candidateType)                           ! loop through nodes on face
        candidateFaceNode = mesh_element(4+FE_face(n,candidateFace,candidateType),candidateElem)
        if (all(myFaceNodes /= candidateFaceNode)) then                                             ! candidate node does not match any of my face nodes
          checkTwins = .true.                                                                       ! perhaps the twin nodes do match
          exit
        endif
      enddo
      if(checkTwins) then
checkCandidateFaceTwins: do dir = 1,3
          do n = 1,FE_NmatchingNodesPerFace(candidateFace,candidateType)                       ! loop through nodes on face
            candidateFaceNode = mesh_element(4+FE_face(n,candidateFace,candidateType),candidateElem)
            if (all(myFaceNodes /= mesh_nodeTwins(dir,candidateFaceNode))) then                     ! node twin does not match either
              if (dir == 3) then
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

end subroutine mesh_build_ipNeighborhood


!--------------------------------------------------------------------------------------------------
!> @brief mapping of FE element types to internal representation
!--------------------------------------------------------------------------------------------------
integer function FE_mapElemtype(what)
 use IO, only: IO_lc, IO_error

 
 character(len=*), intent(in) :: what

 select case (IO_lc(what))
    case ( 'cpe4', &
           'cpe4t')
      FE_mapElemtype = 3            ! Arbitrary Quadrilateral Plane-strain
    case ( 'cpe8', &
           'cpe8t')
      FE_mapElemtype = 4            ! Plane Strain, Eight-node Distorted Quadrilateral
    case ( 'c3d4', &
           'c3d4t')
      FE_mapElemtype = 6            ! Three-dimensional Four-node Tetrahedron
    case ( 'c3d6', &
           'c3d6t')
      FE_mapElemtype = 9            ! Three-dimensional Arbitrarily Distorted Pentahedral
    case ( 'c3d8r', &
           'c3d8rt')
      FE_mapElemtype = 10           ! Three-dimensional Arbitrarily Distorted linear hexahedral with reduced integration
    case ( 'c3d8', &
           'c3d8t')
      FE_mapElemtype = 11           ! Three-dimensional Arbitrarily Distorted Brick
    case ( 'c3d20r', &
           'c3d20rt')
      FE_mapElemtype = 12           ! Three-dimensional Arbitrarily Distorted quad hexahedral with reduced integration
    case ( 'c3d20', &
           'c3d20t')
      FE_mapElemtype = 13           ! Three-dimensional Arbitrarily Distorted quadratic hexahedral
    case default
      call IO_error(error_ID=190,ext_msg=IO_lc(what))
 end select

end function FE_mapElemtype





!--------------------------------------------------------------------------------------------------
!> @brief get properties of different types of finite elements
!> @details assign globals: FE_nodesAtIP, FE_ipNeighbor, FE_cellnodeParentnodeWeights, FE_subNodeOnIPFace
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_FEdata

 
 integer :: me
 allocate(FE_nodesAtIP(FE_maxmaxNnodesAtIP,FE_maxNips,FE_Ngeomtypes), source=0)
 allocate(FE_ipNeighbor(FE_maxNipNeighbors,FE_maxNips,FE_Ngeomtypes), source=0)
 allocate(FE_cell(FE_maxNcellnodesPerCell,FE_maxNips,FE_Ngeomtypes),  source=0)
 allocate(FE_cellnodeParentnodeWeights(FE_maxNnodes,FE_maxNcellnodes,FE_Nelemtypes), source=0.0_pReal)
 allocate(FE_cellface(FE_maxNcellnodesPerCellface,FE_maxNcellfaces,FE_Ncelltypes),  source=0)


 !*** fill FE_nodesAtIP with data ***

 me = 0

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    1,2,3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
    1,  &
    2,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
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

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
    1,2,3,4   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4,  &
    5,  &
    6   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
 FE_nodesAtIP(1:FE_maxNnodesAtIP(me),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
    1,2,3,4,5,6,7,8   &
    ],pInt),[FE_maxNnodesAtIP(me),FE_Nips(me)])

 me = me + 1
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

 me = me + 1
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
 me = 0

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    -2,-3,-1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 125 (2D 6node 3ip)
    reshape(int([&
     2,-3, 3,-1,  &
    -2, 1, 3,-1,  &
     2,-3,-2, 1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element  11 (2D 4node 4ip)
    reshape(int([&
     2,-4, 3,-1,  &
    -2, 1, 4,-1,  &
     4,-4,-3, 1,  &
    -2, 3,-3, 2   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
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

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 134 (3D 4node 1ip)
    reshape(int([&
    -1,-2,-3,-4   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 127 (3D 10node 4ip)
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -2, 1, 3,-2, 4,-1,  &
     2,-4,-3, 1, 4,-1,  &
     2,-4, 3,-2,-3, 1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 136 (3D 6node 6ip)
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -3, 1, 3,-2, 5,-1,  &
     2,-4,-3, 1, 6,-1,  &
     5,-4, 6,-2,-5, 1,  &
    -3, 4, 6,-2,-5, 2,  &
     5,-4,-3, 4,-5, 3   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_ipNeighbor(1:FE_NipNeighbors(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element 117 (3D 8node 1ip)
    reshape(int([&
    -3,-5,-4,-2,-6,-1   &
    ],pInt),[FE_NipNeighbors(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
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

 me = me + 1
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
 me = 0

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   6 (2D 3node 1ip)
    reshape(int([&
    1,2,3   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   125 (2D 6node 3ip)
    reshape(int([&
    1, 4, 7, 6,   &
    2, 5, 7, 4,   &
    3, 6, 7, 5    &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   11 (2D 4node 4ip)
    reshape(int([&
    1, 5, 9, 8,   &
    5, 2, 6, 9,   &
    8, 9, 7, 4,   &
    9, 6, 3, 7    &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
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

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   134 (3D 4node 1ip)
    reshape(int([&
    1, 2, 3, 4   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   127 (3D 10node 4ip)
    reshape(int([&
    1, 5,11, 7, 8,12,15,14,  &
    5, 2, 6,11,12, 9,13,15,  &
    7,11, 6, 3,14,15,13,10,  &
    8,12,15, 4, 4, 9,13,10   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   136 (3D 6node 6ip)
    reshape(int([&
    1, 7,16, 9,10,17,21,19,  &
    7, 2, 8,16,17,11,18,21,  &
    9,16, 8, 3,19,21,18,12,  &
   10,17,21,19, 4,13,20,15,  &
   17,11,18,21,13, 5,14,20,  &
   19,21,18,12,15,20,14, 6   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
 FE_cell(1:FE_NcellnodesPerCell(FE_celltype(me)),1:FE_Nips(me),me) = &  ! element   117 (3D 8node 1ip)
    reshape(int([&
    1, 2, 3, 4, 5, 6, 7, 8   &
    ],pInt),[FE_NcellnodesPerCell(FE_celltype(me)),FE_Nips(me)])

 me = me + 1
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

 me = me + 1
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
 me = 0

 me = me + 1
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element   6 (2D 3node 1ip)
    reshape(real([&
    1, 0, 0,  &
    0, 1, 0,  &
    0, 0, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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

 me = me + 1
 FE_cellnodeParentnodeWeights(1:FE_Nnodes(me),1:FE_Ncellnodes(FE_geomtype(me)),me) = &  ! element 134 (3D 4node 1ip)
    reshape(real([&
    1, 0, 0, 0,  &
    0, 1, 0, 0,  &
    0, 0, 1, 0,  &
    0, 0, 0, 1   &
    ],pReal),[FE_Nnodes(me),FE_Ncellnodes(FE_geomtype(me))])

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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

 me = me + 1
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
 me = 0

 me = me + 1
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 2D 3node, VTK_TRIANGLE (5)
    reshape(int([&
    2,3,  &
    3,1,  &
    1,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 2D 4node, VTK_QUAD (9)
    reshape(int([&
    2,3,  &
    4,1,  &
    3,4,  &
    1,2   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1
 FE_cellface(1:FE_NcellnodesPerCellface(me),1:FE_NipNeighbors(me),me) = &                           ! 3D 4node, VTK_TETRA (10)
    reshape(int([&
    1,3,2,  &
    1,2,4,  &
    2,3,4,  &
    1,4,3   &
    ],pInt),[FE_NcellnodesPerCellface(me),FE_NipNeighbors(me)])

 me = me + 1
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


!--------------------------------------------------------------------------------------------------
!> @brief Gives the FE to CP ID mapping by binary search through lookup array
!! valid questions (what) are 'elem', 'node'
!--------------------------------------------------------------------------------------------------
integer function mesh_FEasCP(what,myID)
 use IO, only: &
   IO_lc

 character(len=*), intent(in) :: what
 integer,    intent(in) :: myID

 integer, dimension(:,:), pointer :: lookupMap
 integer :: lower,upper,center

 mesh_FEasCP = 0
 select case(IO_lc(what(1:4)))
   case('elem')
     lookupMap => mesh_mapFEtoCPelem
   case('node')
     lookupMap => mesh_mapFEtoCPnode
   case default
     return
 endselect

 lower = 1
 upper = int(size(lookupMap,2),pInt)

 if (lookupMap(1,lower) == myID) then                                                          ! check at bounds QUESTION is it valid to extend bounds by 1 and just do binary search w/o init check at bounds?
   mesh_FEasCP = lookupMap(2,lower)
   return
 elseif (lookupMap(1,upper) == myID) then
   mesh_FEasCP = lookupMap(2,upper)
   return
 endif
 binarySearch: do while (upper-lower > 1)
   center = (lower+upper)/2
   if (lookupMap(1,center) < myID) then
     lower = center
   elseif (lookupMap(1,center) > myID) then
     upper = center
   else
     mesh_FEasCP = lookupMap(2,center)
     exit
   endif
 enddo binarySearch

end function mesh_FEasCP

end module mesh
