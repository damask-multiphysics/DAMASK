!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solver MSC.Marc
!--------------------------------------------------------------------------------------------------
module mesh
 use prec, only: pReal, pInt
 use mesh_base

 implicit none
 private
 integer(pInt), public, protected :: &
   mesh_elemType, &                                                                                 !< Element type of the mesh (only support homogeneous meshes)
   mesh_Nnodes, &                                                                                   !< total number of nodes in mesh
   mesh_Ncellnodes, &                                                                               !< total number of cell nodes in mesh (including duplicates)
   mesh_Ncells, &                                                                                   !< total number of cells in mesh
   mesh_maxNsharedElems                                                                             !< max number of CP elements sharing a node

 integer(pInt), dimension(:,:), allocatable, public, protected :: &
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


integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_cellnodeParent                                                                              !< cellnode's parent element ID, cellnode's intra-element ID

 integer(pInt),dimension(:,:,:), allocatable, private :: &
   mesh_cell                                                                                        !< cell connectivity for each element,ip/cell

 integer(pInt), dimension(:,:,:), allocatable, private :: &
   FE_cellface                                                                                      !< list of intra-cell cell node IDs that constitute the cell faces of a specific type of cell


! These definitions should actually reside in the FE-solver specific part (different for MARC/ABAQUS)
! Hence, I suggest to prefix with "FE_"

 integer(pInt), parameter, public :: &
   FE_Nelemtypes = 13_pInt, &
   FE_Ngeomtypes = 10_pInt, &
   FE_Ncelltypes = 4_pInt, &
   FE_maxNipNeighbors = 6_pInt, &
   FE_maxmaxNnodesAtIP = 8_pInt, &                                                                  !< max number of (equivalent) nodes attached to an IP
   FE_maxNmatchingNodesPerFace = 4_pInt, &
   FE_maxNfaces = 6_pInt, &
   FE_maxNcellnodes = 64_pInt, &
   FE_maxNcellnodesPerCell = 8_pInt, &
   FE_maxNcellfaces = 6_pInt, &
   FE_maxNcellnodesPerCellface = 4_pInt


 integer(pInt), dimension(FE_Ngeomtypes), parameter, private :: FE_Nfaces = &                        !< number of faces of a specific type of element geometry
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


 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NcellnodesPerCellface = &        !< number of cell nodes per cell face in a specific cell type
 int([&
      2, & ! (2D 3node)
      2, & ! (2D 4node)
      3, & ! (3D 4node)
      4  & ! (3D 8node)
  ],pInt)

 integer(pInt), dimension(FE_Ncelltypes), parameter, private :: FE_NipNeighbors = &                  !< number of ip neighbors / cell faces in a specific cell type
 int([&
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      6  & ! (3D 8node)
  ],pInt)


 integer(pInt), private :: &
   mesh_Nelems, &                                                                                   !< total number of elements in mesh (including non-DAMASK elements)
   mesh_NelemSets
 character(len=64), dimension(:), allocatable, private :: &
   mesh_nameElemSet
   
 integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
 integer(pInt), dimension(:,:), allocatable, target, private :: &
   mesh_mapFEtoCPelem, &                                                                            !< [sorted FEid, corresponding CPid]
   mesh_mapFEtoCPnode                                                                               !< [sorted FEid, corresponding CPid]


 integer(pInt), private :: &
   MarcVersion, &                                                                                   !< Version of input file format (Marc only)
   hypoelasticTableStyle, &                                                                         !< Table style (Marc only)
   initialcondTableStyle                                                                            !< Table style (Marc only)
 integer(pInt), dimension(:), allocatable, private :: &
   Marc_matNumber                                                                                   !< array of material numbers for hypoelastic material (Marc only)

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
   mesh_marc_get_fileFormat, &
   mesh_marc_get_tableStyles, &
   mesh_marc_get_matNumber, &
   mesh_marc_count_nodesAndElements, &
   mesh_marc_count_elementSets, &
   mesh_marc_map_elementSets, &
   mesh_marc_map_Elements, &
   mesh_marc_map_nodes, &
   mesh_marc_build_nodes, &
   mesh_marc_build_elements

type, public, extends(tMesh) :: tMesh_marc
 
 contains 
   procedure, pass(self) :: tMesh_marc_init
   generic, public :: init => tMesh_marc_init
end type tMesh_marc

 type(tMesh_marc), public, protected :: theMesh
 
 
contains

subroutine tMesh_marc_init(self,elemType,nodes)
 
 implicit none
 class(tMesh_marc) :: self
 real(pReal), dimension(:,:), intent(in) :: nodes
 integer(pInt), intent(in) :: elemType
 
 call self%tMesh%init('mesh',elemType,nodes)
 
end subroutine tMesh_marc_init

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
   calcMode, &
   FEsolving_execElem, &
   FEsolving_execIP

 implicit none
 integer(pInt), intent(in) :: el, ip
  
 integer(pInt), parameter  :: FILEUNIT = 222_pInt
 integer(pInt) :: j, fileFormatVersion, elemType
 integer(pInt) :: &
   mesh_maxNelemInSet, &
   mesh_NcpElems
 logical :: myDebug

 write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'

 mesh_unitlength = numerics_unitlength                                                              ! set physical extent of a length unit in mesh

 myDebug = (iand(debug_level(debug_mesh),debug_levelBasic) /= 0_pInt)

 call IO_open_inputFile(FILEUNIT,modelName)                                                         ! parse info from input file...
 if (myDebug) write(6,'(a)') ' Opened input file'; flush(6)
 
 MarcVersion = mesh_marc_get_fileFormat(FILEUNIT)
 fileFormatVersion  = MarcVersion
 if (myDebug) write(6,'(a)') ' Got input file format'; flush(6)
 
 call mesh_marc_get_tableStyles(initialcondTableStyle,hypoelasticTableStyle,FILEUNIT)
 if (myDebug) write(6,'(a)') ' Got table styles'; flush(6)
 
 if (fileFormatVersion > 12) then
   Marc_matNumber = mesh_marc_get_matNumber(FILEUNIT,hypoelasticTableStyle)
   if (myDebug) write(6,'(a)') ' Got hypoleastic material number'; flush(6)
 endif
 
 call mesh_marc_count_nodesAndElements(mesh_nNodes, mesh_nElems, FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted nodes/elements'; flush(6)
 
 call mesh_marc_count_elementSets(mesh_NelemSets,mesh_maxNelemInSet,FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted element sets'; flush(6)

 allocate(mesh_nameElemSet(mesh_NelemSets)); mesh_nameElemSet = 'n/a'
 allocate(mesh_mapElemSet(1_pInt+mesh_maxNelemInSet,mesh_NelemSets),source=0_pInt)
 call mesh_marc_map_elementSets(mesh_nameElemSet,mesh_mapElemSet,FILEUNIT)
 if (myDebug) write(6,'(a)') ' Mapped element sets'; flush(6)
 
 mesh_NcpElems =  mesh_nElems
 if (myDebug) write(6,'(a)') ' Counted CP elements'; flush(6)
 
 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems), source = 0_pInt)
 call mesh_marc_map_elements(hypoelasticTableStyle,mesh_nameElemSet,mesh_mapElemSet,mesh_NcpElems,FILEUNIT)
 if (myDebug) write(6,'(a)') ' Mapped elements'; flush(6)
 
 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes),source=0_pInt)
 call mesh_marc_map_nodes(mesh_Nnodes,FILEUNIT)                                                                  !ToDo: don't work on global variables
 if (myDebug) write(6,'(a)') ' Mapped nodes'; flush(6)
 
 call mesh_marc_build_nodes(FILEUNIT)                                                               !ToDo: don't work on global variables
 mesh_node = mesh_node0
 if (myDebug) write(6,'(a)') ' Built nodes'; flush(6)
 
 elemType = mesh_marc_count_cpSizes(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Counted CP sizes'; flush(6)
 
 call theMesh%init(elemType,mesh_node0)
 call theMesh%setNelems(mesh_NcpElems)
 
 call mesh_marc_build_elements(FILEUNIT)
 if (myDebug) write(6,'(a)') ' Built elements'; flush(6)
 call mesh_get_damaskOptions(mesh_periodicSurface,FILEUNIT)
 if (myDebug) write(6,'(a)') ' Got DAMASK options'; flush(6)
 close (FILEUNIT)
  
 
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

 if (usePingPong .and. (mesh_Nelems /= theMesh%nElems)) &
   call IO_error(600_pInt)                                                                          ! ping-pong must be disabled when having non-DAMASK elements
 if (debug_e < 1 .or. debug_e > theMesh%nElems) &
   call IO_error(602_pInt,ext_msg='element')                                                        ! selected element does not exist
 if (debug_i < 1 .or. debug_i > theMesh%elem%nIPs) &
   call IO_error(602_pInt,ext_msg='IP')                                                             ! selected element does not have requested IP

 FEsolving_execElem = [ 1_pInt,theMesh%nElems ]                                                      ! parallel loop bounds set to comprise all DAMASK elements
 allocate(FEsolving_execIP(2_pInt,theMesh%nElems), source=1_pInt)                                    ! parallel loop bounds set to comprise from first IP...
 FEsolving_execIP(2,:) = theMesh%elem%nIPs

 allocate(calcMode(theMesh%elem%nIPs,theMesh%nElems))
 calcMode = .false.                                                                                 ! pretend to have collected what first call is asking (F = I)
 calcMode(ip,mesh_FEasCP('elem',el)) = .true.                                                       ! first ip,el needs to be already pingponged to "calc"

 theMesh%homogenizationAt  = mesh_element(3,:)
 theMesh%microstructureAt  = mesh_element(4,:)

end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Figures out version of Marc input file format
!--------------------------------------------------------------------------------------------------
integer(pInt) function mesh_marc_get_fileFormat(fileUnit)
 use IO, only: &
   IO_lc, &
   IO_intValue, &
   IO_stringValue, &
   IO_stringPos

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line


 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)
   
   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'version') then
     mesh_marc_get_fileFormat = IO_intValue(line,chunkPos,2_pInt)
     exit
   endif
 enddo

620 end function mesh_marc_get_fileFormat


!--------------------------------------------------------------------------------------------------
!> @brief Figures out table styles for initial cond and hypoelastic
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_get_tableStyles(initialcond, hypoelastic,fileUnit)
 use IO, only: &
   IO_lc, &
   IO_intValue, &
   IO_stringValue, &
   IO_stringPos

 implicit none
 integer(pInt), intent(out) :: initialcond, hypoelastic
 integer(pInt), intent(in)  :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line

 initialcond = 0_pInt
 hypoelastic = 0_pInt

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)

   if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'table' .and. chunkPos(1_pInt) > 5) then
     initialcond = IO_intValue(line,chunkPos,4_pInt)
     hypoelastic = IO_intValue(line,chunkPos,5_pInt)
     exit
   endif
 enddo

620 end subroutine mesh_marc_get_tableStyles


!--------------------------------------------------------------------------------------------------
!> @brief Figures out material number of hypoelastic material
!--------------------------------------------------------------------------------------------------
function mesh_marc_get_matNumber(fileUnit,tableStyle)
 use IO, only: &
   IO_lc, &
   IO_intValue, &
   IO_stringValue, &
   IO_stringPos

 implicit none
 integer(pInt), intent(in) :: fileUnit, tableStyle
 integer(pInt), dimension(:), allocatable :: mesh_marc_get_matNumber

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: i, j, data_blocks
 character(len=300) line


 rewind(fileUnit)

 data_blocks = 1_pInt
 do
  read (fileUnit,'(A300)',END=620) line
  chunkPos = IO_stringPos(line)
  
  if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'hypoelastic') then
    read (fileUnit,'(A300)',END=620) line
    if (len(trim(line))/=0_pInt) then
      chunkPos = IO_stringPos(line)
      data_blocks = IO_intValue(line,chunkPos,1_pInt)
    endif
    allocate(mesh_marc_get_matNumber(data_blocks), source = 0_pInt)
    do i=1_pInt,data_blocks                                                                         ! read all data blocks
      read (fileUnit,'(A300)',END=620) line
      chunkPos = IO_stringPos(line)
      mesh_marc_get_matNumber(i) = IO_intValue(line,chunkPos,1_pInt)
      do j=1_pint,2_pInt + tableStyle                                                               ! read 2 or 3 remaining lines of data block
        read (fileUnit,'(A300)') line
      enddo
    enddo
    exit
  endif
 enddo

620 end function mesh_marc_get_matNumber


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of nodes and elements
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_count_nodesAndElements(nNodes, nElems, fileUnit)
 use IO, only: &
   IO_lc, &
   IO_stringValue, &
   IO_stringPos, &
   IO_IntValue

 implicit none
 integer(pInt), intent(in)  :: fileUnit
 integer(pInt), intent(out) :: nNodes, nElems
 
 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line

 nNodes = 0_pInt
 nElems = 0_pInt

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)

   if ( IO_lc(IO_StringValue(line,chunkPos,1_pInt)) == 'sizing') &
     nElems = IO_IntValue (line,chunkPos,3_pInt)
   if ( IO_lc(IO_StringValue(line,chunkPos,1_pInt)) == 'coordinates') then
     read (fileUnit,'(A300)') line
     chunkPos = IO_stringPos(line)
     nNodes = IO_IntValue (line,chunkPos,2_pInt)
     exit                                                                                           ! assumes that "coordinates" comes later in file
   endif
 enddo

620 end subroutine mesh_marc_count_nodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of element sets in mesh.
!--------------------------------------------------------------------------------------------------
 subroutine mesh_marc_count_elementSets(nElemSets,maxNelemInSet,fileUnit)
 use IO, only: &
   IO_lc, &
   IO_stringValue, &
   IO_stringPos, &
   IO_countContinuousIntValues

 implicit none
 integer(pInt), intent(in) :: fileUnit
 integer(pInt), intent(out) :: nElemSets, maxNelemInSet

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line

 nElemSets     = 0_pInt
 maxNelemInSet = 0_pInt

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)

   if ( IO_lc(IO_StringValue(line,chunkPos,1_pInt)) == 'define' .and. &
        IO_lc(IO_StringValue(line,chunkPos,2_pInt)) == 'element' ) then
     nElemSets = nElemSets + 1_pInt
     maxNelemInSet = max(maxNelemInSet, IO_countContinuousIntValues(fileUnit))
   endif
 enddo

620 end subroutine mesh_marc_count_elementSets


!--------------------------------------------------------------------------------------------------
!> @brief map element sets
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_elementSets(nameElemSet,mapElemSet,fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_continuousIntValues

 implicit none
 integer(pInt), intent(in) :: fileUnit
 character(len=64), dimension(:), intent(out) :: &
   nameElemSet
 integer(pInt), dimension(:,:), intent(out) :: &
   mapElemSet

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) :: elemSet
 
 elemSet = 0_pInt

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=640) line
   chunkPos = IO_stringPos(line)
   if( (IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'define' ) .and. &
       (IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'element' ) ) then
      elemSet = elemSet+1_pInt
      nameElemSet(elemSet)  = trim(IO_stringValue(line,chunkPos,4_pInt))
      mapElemSet(:,elemSet) = IO_continuousIntValues(fileUnit,size(mapElemSet,1)-1,nameElemSet,mapElemSet,size(nameElemSet))
   endif
 enddo

640 end subroutine mesh_marc_map_elementSets


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_elements(tableStyle,nameElemSet,mapElemSet,nElems,fileUnit)

 use math, only: math_sort
 use IO,   only: IO_lc, &
                 IO_intValue, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_continuousIntValues

 implicit none
 integer(pInt), intent(in) :: fileUnit,tableStyle,nElems
 character(len=64), intent(in), dimension(:) :: nameElemSet
 integer(pInt), dimension(:,:), intent(in) :: &
   mapElemSet

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line, &
                       tmp

 integer(pInt), dimension (1_pInt+nElems) :: contInts
 integer(pInt) :: i,cpElem

 cpElem = 0_pInt
 contInts = 0_pInt
 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=660) line
   chunkPos = IO_stringPos(line)
   if (MarcVersion < 13) then                                                                       ! Marc 2016 or earlier
     if( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'hypoelastic' ) then
       do i=1_pInt,3_pInt+TableStyle                                                                ! skip three (or four if new table style!) lines
         read (fileUnit,'(A300)') line
       enddo
       contInts = IO_continuousIntValues(fileUnit,nElems,nameElemSet,&
                                              mapElemSet,size(nameElemSet))
       exit
     endif  
   else                                                                                             ! Marc2017 and later
     if ( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'connectivity') then
       read (fileUnit,'(A300)',END=660) line
       chunkPos = IO_stringPos(line)
       if(any(Marc_matNumber==IO_intValue(line,chunkPos,6_pInt))) then
         do 
           read (fileUnit,'(A300)',END=660) line
           chunkPos = IO_stringPos(line)
           tmp = IO_lc(IO_stringValue(line,chunkPos,1_pInt))
           if (verify(trim(tmp),"0123456789")/=0) then                                              ! found keyword
             exit
           else
             contInts(1) = contInts(1) + 1_pInt  
             read (tmp,*) contInts(contInts(1)+1)     
           endif
         enddo
       endif  
     endif
   endif    
 enddo    
660 do i = 1_pInt,contInts(1)
      cpElem = cpElem+1_pInt
      mesh_mapFEtoCPelem(1,cpElem) = contInts(1_pInt+i)
      mesh_mapFEtoCPelem(2,cpElem) = cpElem
    enddo
 
call math_sort(mesh_mapFEtoCPelem,1_pInt,int(size(mesh_mapFEtoCPelem,2_pInt),pInt))                ! should be mesh_NcpElems

end subroutine mesh_marc_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_nodes(nNodes,fileUnit)

 use math, only: math_sort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_fixedIntValue

 implicit none
 integer(pInt), intent(in) :: fileUnit, nNodes

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line

 integer(pInt), dimension (nNodes) :: node_count
 integer(pInt) :: i

 node_count = 0_pInt

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=650) line
   chunkPos = IO_stringPos(line)
   if( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'coordinates' ) then
     read (fileUnit,'(A300)') line                                                                  ! skip crap line
     do i = 1_pInt,nNodes
       read (fileUnit,'(A300)') line
       mesh_mapFEtoCPnode(1_pInt,i) = IO_fixedIntValue (line,[0_pInt,10_pInt],1_pInt)
       mesh_mapFEtoCPnode(2_pInt,i) = i
     enddo
     exit
   endif
 enddo

650 call math_sort(mesh_mapFEtoCPnode,1_pInt,int(size(mesh_mapFEtoCPnode,2_pInt),pInt))

end subroutine mesh_marc_map_nodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_build_nodes(fileUnit)

 use IO, only: &
   IO_lc, &
   IO_stringValue, &
   IO_stringPos, &
   IO_fixedIntValue, &
   IO_fixedNoEFloatValue

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), dimension(5), parameter :: node_ends = int([0,10,30,50,70],pInt)
 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) :: i,j,m

 allocate ( mesh_node0 (3,mesh_Nnodes), source=0.0_pReal)

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=670) line
   chunkPos = IO_stringPos(line)
   if( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'coordinates' ) then
     read (fileUnit,'(A300)') line                                                                  ! skip crap line
     do i=1_pInt,mesh_Nnodes
       read (fileUnit,'(A300)') line
       m = mesh_FEasCP('node',IO_fixedIntValue(line,node_ends,1_pInt))
       do j = 1_pInt,3_pInt
         mesh_node0(j,m) = mesh_unitlength * IO_fixedNoEFloatValue(line,node_ends,j+1_pInt)
       enddo
     enddo
     exit
   endif
 enddo

670 mesh_node = mesh_node0

end subroutine mesh_marc_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets maximum count of nodes, IPs, IP neighbors, and cellnodes among cpElements.
!! Sets global values 'mesh_maxNnodes', 'mesh_maxNips', 'mesh_maxNipNeighbors',
!! and 'mesh_maxNcellnodes'
!--------------------------------------------------------------------------------------------------
integer(pInt) function mesh_marc_count_cpSizes(fileUnit)

 use IO,   only: IO_lc, &
                 IO_error, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_skipChunks
 use element

 implicit none
 integer(pInt), intent(in) :: fileUnit

 type(tElement) :: tempEl
 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer(pInt) :: i,t,g,e,c

 t = -1_pInt

 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=630) line
   chunkPos = IO_stringPos(line)
   if( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'connectivity' ) then
     read (fileUnit,'(A300)') line                                                                  ! Garbage line
     do i=1_pInt,mesh_Nelems                                                                        ! read all elements
       read (fileUnit,'(A300)') line
       chunkPos = IO_stringPos(line)                                                                ! limit to id and type
       if (t == -1_pInt) then
         t = FE_mapElemtype(IO_stringValue(line,chunkPos,2_pInt))
         call tempEl%init(t)
         mesh_marc_count_cpSizes = t
       else
         if (t /= FE_mapElemtype(IO_stringValue(line,chunkPos,2_pInt))) call IO_error(0_pInt)       !ToDo: error message
       endif
       call IO_skipChunks(fileUnit,tempEl%nNodes-(chunkPos(1_pInt)-2_pInt))
     enddo
     exit
   endif
 enddo

630 end function mesh_marc_count_cpSizes


!--------------------------------------------------------------------------------------------------
!> @brief Store FEid, type, mat, tex, and node list per element.
!! Allocates global array 'mesh_element'
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_build_elements(fileUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_fixedNoEFloatValue, &
                 IO_skipChunks, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_continuousIntValues, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) line

 integer(pInt), dimension(1_pInt+theMesh%nElems) :: contInts
 integer(pInt) :: i,j,t,sv,myVal,e,nNodesAlreadyRead

 allocate(mesh_element(4_pInt+theMesh%elem%nNodes,theMesh%nElems), source=0_pInt)
 mesh_elemType = -1_pInt


 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)
   if( IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'connectivity' ) then
     read (fileUnit,'(A300)',END=620) line                                                               ! garbage line
     do i = 1_pInt,mesh_Nelems
       read (fileUnit,'(A300)',END=620) line
       chunkPos = IO_stringPos(line)
       e = mesh_FEasCP('elem',IO_intValue(line,chunkPos,1_pInt))
       if (e /= 0_pInt) then                                                                        ! disregard non CP elems
         mesh_element(1,e) = -1_pInt                                                                ! DEPRECATED
         t = FE_mapElemtype(IO_StringValue(line,chunkPos,2_pInt))                                   ! elem type
         if (mesh_elemType /= t .and. mesh_elemType /= -1_pInt) &
           call IO_error(191,el=t,ip=mesh_elemType)
         mesh_elemType = t
         mesh_element(2,e) = t
         nNodesAlreadyRead = 0_pInt
         do j = 1_pInt,chunkPos(1)-2_pInt
           mesh_element(4_pInt+j,e) = mesh_FEasCP('node',IO_IntValue(line,chunkPos,j+2_pInt))          ! CP ids of nodes
         enddo
         nNodesAlreadyRead = chunkPos(1) - 2_pInt
         do while(nNodesAlreadyRead < theMesh%elem%nNodes)                                                 ! read on if not all nodes in one line
           read (fileUnit,'(A300)',END=620) line
           chunkPos = IO_stringPos(line)
           do j = 1_pInt,chunkPos(1)
             mesh_element(4_pInt+nNodesAlreadyRead+j,e) &
               = mesh_FEasCP('node',IO_IntValue(line,chunkPos,j))                                      ! CP ids of nodes
           enddo
           nNodesAlreadyRead = nNodesAlreadyRead + chunkPos(1)
         enddo
       endif
     enddo
     exit
   endif
 enddo

620 rewind(fileUnit)                                                                                ! just in case "initial state" appears before "connectivity"
 read (fileUnit,'(A300)',END=630) line
 do
   chunkPos = IO_stringPos(line)
   if( (IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == 'initial') .and. &
       (IO_lc(IO_stringValue(line,chunkPos,2_pInt)) == 'state') ) then
     if (initialcondTableStyle == 2_pInt) read (fileUnit,'(A300)',END=630) line                          ! read extra line for new style
     read (fileUnit,'(A300)',END=630) line                                                               ! read line with index of state var
     chunkPos = IO_stringPos(line)
     sv = IO_IntValue(line,chunkPos,1_pInt)                                                            ! figure state variable index
     if( (sv == 2_pInt).or.(sv == 3_pInt) ) then                                                    ! only state vars 2 and 3 of interest
       read (fileUnit,'(A300)',END=630) line                                                             ! read line with value of state var
       chunkPos = IO_stringPos(line)
       do while (scan(IO_stringValue(line,chunkPos,1_pInt),'+-',back=.true.)>1)                        ! is noEfloat value?
         myVal = nint(IO_fixedNoEFloatValue(line,[0_pInt,20_pInt],1_pInt),pInt)                     ! state var's value
         if (initialcondTableStyle == 2_pInt) then
           read (fileUnit,'(A300)',END=630) line                                                         ! read extra line
           read (fileUnit,'(A300)',END=630) line                                                         ! read extra line
         endif
         contInts = IO_continuousIntValues&                                                         ! get affected elements
                   (fileUnit,theMesh%nElems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
         do i = 1_pInt,contInts(1)
           e = mesh_FEasCP('elem',contInts(1_pInt+i))
           mesh_element(1_pInt+sv,e) = myVal
         enddo
         if (initialcondTableStyle == 0_pInt) read (fileUnit,'(A300)',END=630) line                      ! ignore IP range for old table style
         read (fileUnit,'(A300)',END=630) line
         chunkPos = IO_stringPos(line)
       enddo
     endif
   else
     read (fileUnit,'(A300)',END=630) line
   endif
 enddo

630 end subroutine mesh_marc_build_elements


!--------------------------------------------------------------------------------------------------
!> @brief get any additional damask options from input file, sets mesh_periodicSurface
!--------------------------------------------------------------------------------------------------
subroutine mesh_get_damaskOptions(periodic_surface,fileUnit)

use IO, only: &
  IO_lc, &
  IO_stringValue, &
  IO_stringPos

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 character(len=300) :: line
 integer :: myStat
 integer(pInt) :: chunk, Nchunks
 character(len=300) ::  v
 logical, dimension(3) :: periodic_surface
 

 periodic_surface = .false.
 myStat = 0
 rewind(fileUnit)
 do while(myStat == 0)
   read (fileUnit,'(a300)',iostat=myStat) line
   chunkPos = IO_stringPos(line)
   Nchunks = chunkPos(1)
   if (IO_lc(IO_stringValue(line,chunkPos,1_pInt)) == '$damask' .and. Nchunks > 1_pInt) then        ! found keyword for damask option and there is at least one more chunk to read
     select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
       case('periodic')                                                                             ! damask Option that allows to specify periodic fluxes
         do chunk = 3_pInt,Nchunks                                                                  ! loop through chunks (skipping the keyword)
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

 implicit none
 integer(pInt), dimension(:), allocatable :: &
   matchingNode2cellnode
 integer(pInt), dimension(:,:), allocatable :: &
   cellnodeParent
 integer(pInt), dimension(theMesh%elem%Ncellnodes) :: &
   localCellnode2globalCellnode
 integer(pInt) :: &
   e,n,i, &
   matchingNodeID, &
   localCellnodeID

 allocate(mesh_cell(FE_maxNcellnodesPerCell,theMesh%elem%nIPs,theMesh%nElems), source=0_pInt)
 allocate(matchingNode2cellnode(theMesh%nNodes),                            source=0_pInt)
 allocate(cellnodeParent(2_pInt,theMesh%elem%Ncellnodes*theMesh%nElems),       source=0_pInt)
 
 mesh_Ncells = theMesh%nElems*theMesh%elem%nIPs
!--------------------------------------------------------------------------------------------------
! Count cell nodes (including duplicates) and generate cell connectivity list
 mesh_Ncellnodes = 0_pInt

 do e = 1_pInt,theMesh%nElems
   localCellnode2globalCellnode = 0_pInt
   do i = 1_pInt,theMesh%elem%nIPs
     do n = 1_pInt,theMesh%elem%NcellnodesPerCell
       localCellnodeID = theMesh%elem%cell(n,i)
       if (localCellnodeID <= FE_NmatchingNodes(theMesh%elem%geomType)) then                                            ! this cell node is a matching node
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
   e,n,m, &
   localCellnodeID
 real(pReal), dimension(3) :: &
   myCoords

 mesh_build_cellnodes = 0.0_pReal
!$OMP PARALLEL DO PRIVATE(e,localCellnodeID,myCoords)
 do n = 1_pInt,Ncellnodes                                                                           ! loop over cell nodes
   e = mesh_cellnodeParent(1,n)
   localCellnodeID = mesh_cellnodeParent(2,n)
   myCoords = 0.0_pReal
   do m = 1_pInt,theMesh%elem%nNodes
     myCoords = myCoords + nodes(1:3,mesh_element(4_pInt+m,e)) &
                         * theMesh%elem%cellNodeParentNodeWeights(m,localCellnodeID)
   enddo
   mesh_build_cellnodes(1:3,n) = myCoords / sum(theMesh%elem%cellNodeParentNodeWeights(:,localCellnodeID))
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

  allocate(mesh_ipVolume(theMesh%elem%nIPs,theMesh%nElems),source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,m,subvolume)
   do e = 1_pInt,theMesh%nElems                                                                      ! loop over cpElems
     t = mesh_element(2_pInt,e)                                                                     ! get element type
     g = theMesh%elem%geomType
     c = theMesh%elem%cellType
     select case (c)

       case (1_pInt)                                                                                ! 2D 3node
         forall (i = 1_pInt:theMesh%elem%nIPs) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_areaTriangle(mesh_cellnode(1:3,mesh_cell(1,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(3,i,e)))

       case (2_pInt)                                                                                ! 2D 4node
         forall (i = 1_pInt:theMesh%elem%nIPs) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_areaTriangle(mesh_cellnode(1:3,mesh_cell(1,i,e)), &            ! here we assume a planar shape, so division in two triangles suffices
                                                  mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(3,i,e))) &
                              + math_areaTriangle(mesh_cellnode(1:3,mesh_cell(3,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(4,i,e)), &
                                                  mesh_cellnode(1:3,mesh_cell(1,i,e)))

       case (3_pInt)                                                                                ! 3D 4node
         forall (i = 1_pInt:theMesh%elem%nIPs) &                                                           ! loop over ips=cells in this element
           mesh_ipVolume(i,e) = math_volTetrahedron(mesh_cellnode(1:3,mesh_cell(1,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(2,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(3,i,e)), &
                                                    mesh_cellnode(1:3,mesh_cell(4,i,e)))

       case (4_pInt)                                                                                ! 3D 8node
         m = FE_NcellnodesPerCellface(c)
         do i = 1_pInt,theMesh%elem%nIPs                                                                   ! loop over ips=cells in this element
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
   allocate(mesh_ipCoordinates(3,theMesh%elem%nIPs,theMesh%nElems),source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,myCoords)
 do e = 1_pInt,theMesh%nElems                                                                        ! loop over cpElems
   t = mesh_element(2_pInt,e)                                                                       ! get element type
   g = theMesh%elem%geomType
   c = theMesh%elem%cellType
   do i = 1_pInt,theMesh%elem%nIPs
     myCoords = 0.0_pReal
     do n = 1_pInt,theMesh%elem%nCellnodesPerCell
       myCoords = myCoords + mesh_cellnode(1:3,mesh_cell(n,i,e))
     enddo
     mesh_ipCoordinates(1:3,i,e) = myCoords / real(theMesh%elem%nCellnodesPerCell,pReal)
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
 g = theMesh%elem%geomType
 c = theMesh%elem%cellType
 mesh_cellCenterCoordinates = 0.0_pReal
 do n = 1_pInt,theMesh%elem%nCellnodesPerCell
   mesh_cellCenterCoordinates = mesh_cellCenterCoordinates + mesh_cellnode(1:3,mesh_cell(n,ip,el))
 enddo
 mesh_cellCenterCoordinates = mesh_cellCenterCoordinates / real(theMesh%elem%nCellnodesPerCell,pReal)

 end function mesh_cellCenterCoordinates





!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipAreas
 use math, only: &
   math_cross

 implicit none
 integer(pInt) :: e,t,g,c,i,f,n,m
 real(pReal), dimension (3,FE_maxNcellnodesPerCellface) :: nodePos, normals
 real(pReal), dimension(3) :: normal

 allocate(mesh_ipArea(theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems), source=0.0_pReal)
 allocate(mesh_ipAreaNormal(3_pInt,theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems), source=0.0_pReal)

 !$OMP PARALLEL DO PRIVATE(t,g,c,nodePos,normal,normals)
   do e = 1_pInt,theMesh%nElems                                                                      ! loop over cpElems
     t = mesh_element(2_pInt,e)                                                                     ! get element type
     g = theMesh%elem%geomType
     c = theMesh%elem%cellType
     select case (c)

       case (1_pInt,2_pInt)                                                                         ! 2D 3 or 4 node
         do i = 1_pInt,theMesh%elem%nIPs
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
         do i = 1_pInt,theMesh%elem%nIPs
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             normal = math_cross(nodePos(1:3,2) - nodePos(1:3,1), &
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
         do i = 1_pInt,theMesh%elem%nIPs
           do f = 1_pInt,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(FE_cellface(n,f,c),i,e))
             forall(n = 1_pInt:FE_NcellnodesPerCellface(c)) &
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

 do e = 1_pInt,theMesh%nElems
   g = theMesh%elem%geomType
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

 do e = 1_pInt,theMesh%nElems
   g = theMesh%elem%geomType
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

 allocate(mesh_ipNeighborhood(3,theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems))
 mesh_ipNeighborhood = 0_pInt


 do myElem = 1_pInt,theMesh%nElems                                                                    ! loop over cpElems
   myType = theMesh%elem%geomType
   do myIP = 1_pInt,theMesh%elem%nIPs

     do neighbor = 1_pInt,FE_NipNeighbors(theMesh%elem%cellType)                                       ! loop over neighbors of IP
       neighboringIPkey = theMesh%elem%IPneighbor(neighbor,myIP)

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
           neighboringType = theMesh%elem%geomType

           !*** trivial solution if neighbor has only one IP

           if (theMesh%elem%nIPs == 1_pInt) then
             mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
             mesh_ipNeighborhood(2,neighbor,myIP,myElem) = 1_pInt
             cycle
           endif

           !*** find those nodes which build the link to the neighbor

           NlinkedNodes = 0_pInt
           linkedNodes = 0_pInt
           do a = 1_pInt,theMesh%elem%maxNnodeAtIP
             anchor = theMesh%elem%NnodeAtIP(a,myIP)
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

 checkCandidateIP: do candidateIP = 1_pInt,theMesh%elem%nIPs
             NmatchingNodes = 0_pInt
             matchingNodes = 0_pInt
             do a = 1_pInt,theMesh%elem%maxNnodeAtIP
               anchor = theMesh%elem%NnodeAtIP(a,candidateIP)
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
 do myElem = 1_pInt,theMesh%nElems                                                                    ! loop over cpElems
   myType = theMesh%elem%geomType
   do myIP = 1_pInt,theMesh%elem%nIPs
     do neighbor = 1_pInt,FE_NipNeighbors(theMesh%elem%cellType)                                       ! loop over neighbors of IP
       neighboringElem = mesh_ipNeighborhood(1,neighbor,myIP,myElem)
       neighboringIP   = mesh_ipNeighborhood(2,neighbor,myIP,myElem)
       if (neighboringElem > 0_pInt .and. neighboringIP > 0_pInt) then                               ! if neighbor exists ...
         neighboringType = theMesh%elem%geomType
         do pointingToMe = 1_pInt,FE_NipNeighbors(theMesh%elem%cellType)                      ! find neighboring index that points from my neighbor to myself
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
 
 contains
 
 !--------------------------------------------------------------------------------------------------
!> @brief find face-matching element of same type
!--------------------------------------------------------------------------------------------------
subroutine mesh_faceMatch(elem, face ,matchingElem, matchingFace)

implicit none
integer(pInt), intent(out) ::     matchingElem, &                                                   ! matching CP element ID
                                  matchingFace                                                      ! matching face ID
integer(pInt), intent(in) ::      face, &                                                           ! face ID
                                  elem                                                              ! CP elem ID
integer(pInt), dimension(FE_NmatchingNodesPerFace(face,theMesh%elem%geomType)) :: &
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
myType =theMesh%elem%geomType

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
    candidateType = theMesh%elem%geomType
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

end subroutine mesh_build_ipNeighborhood


!--------------------------------------------------------------------------------------------------
!> @brief mapping of FE element types to internal representation
!--------------------------------------------------------------------------------------------------
integer(pInt) function FE_mapElemtype(what)
 use IO, only: IO_lc, IO_error

 implicit none
 character(len=*), intent(in) :: what

 select case (IO_lc(what))
    case (   '6')
      FE_mapElemtype = 1_pInt            ! Two-dimensional Plane Strain Triangle
    case ( '155', &
           '125', &
           '128')
      FE_mapElemtype = 2_pInt            ! Two-dimensional Plane Strain triangle (155: cubic shape function, 125/128: second order isoparametric)
    case ( '11')
      FE_mapElemtype = 3_pInt            ! Arbitrary Quadrilateral Plane-strain
    case ( '27')
      FE_mapElemtype = 4_pInt            ! Plane Strain, Eight-node Distorted Quadrilateral
    case ( '54')
      FE_mapElemtype = 5_pInt            ! Plane Strain, Eight-node Distorted Quadrilateral with reduced integration
    case ( '134')
      FE_mapElemtype = 6_pInt            ! Three-dimensional Four-node Tetrahedron
    case ( '157')
      FE_mapElemtype = 7_pInt            ! Three-dimensional, Low-order, Tetrahedron, Herrmann Formulations
    case ( '127')
      FE_mapElemtype = 8_pInt            ! Three-dimensional Ten-node Tetrahedron
    case ( '136')
      FE_mapElemtype = 9_pInt            ! Three-dimensional Arbitrarily Distorted Pentahedral
    case ( '117', &
           '123')
      FE_mapElemtype = 10_pInt           ! Three-dimensional Arbitrarily Distorted linear hexahedral with reduced integration
    case ( '7')
      FE_mapElemtype = 11_pInt           ! Three-dimensional Arbitrarily Distorted Brick
    case ( '57')
      FE_mapElemtype = 12_pInt           ! Three-dimensional Arbitrarily Distorted quad hexahedral with reduced integration
    case ( '21')
      FE_mapElemtype = 13_pInt           ! Three-dimensional Arbitrarily Distorted quadratic hexahedral
    case default
      call IO_error(error_ID=190_pInt,ext_msg=IO_lc(what))
 end select

end function FE_mapElemtype


!--------------------------------------------------------------------------------------------------
!> @brief get properties of different types of finite elements
!> @details assign globals       FE_cellface
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_FEdata

 implicit none
 integer(pInt) :: me
 allocate(FE_cellface(FE_maxNcellnodesPerCellface,FE_maxNcellfaces,FE_Ncelltypes),  source=0_pInt)

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

end module mesh
