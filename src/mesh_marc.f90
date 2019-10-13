!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solver MSC.Marc
!--------------------------------------------------------------------------------------------------
module mesh
  use IO
  use prec
  use math
  use mesh_base
  use DAMASK_interface
  use IO
  use debug
  use numerics
  use FEsolving
  use element
  use discretization
  use geometry_plastic_nonlocal
  use HDF5_utilities
  use results

  implicit none
  private
    
  real(pReal), public, protected :: &
    mesh_unitlength                                                                                  !< physical length of one unit in mesh

!-------------------------------------------------------------------------------------------------- 
! public variables (DEPRECATED)
   
 real(pReal), dimension(:,:,:), allocatable, public :: &
   mesh_ipCoordinates                                                                               !< IP x,y,z coordinates (after deformation!)

 real(pReal), dimension(:,:), allocatable, public :: &
   mesh_cellnode                                                                                    !< cell node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)
!--------------------------------------------------------------------------------------------------

 integer, dimension(:,:), allocatable :: &
   mesh_FEnodes
   
 real(pReal), dimension(:,:), allocatable :: &
   mesh_node, &                                                                                     !< node x,y,z coordinates (after deformation! ONLY FOR MARC!!!
   mesh_node0                                                                                       !< node x,y,z coordinates (initially!)

! -------------------------------------------------------------------------------------------------- 

 type(tMesh) :: theMesh
 

 integer:: &
   mesh_Ncellnodes, &                                                                                  !< total number of cell nodes in mesh (including duplicates)
   mesh_elemType, &                                                                                 !< Element type of the mesh (only support homogeneous meshes)
   mesh_Nnodes, &                                                                                   !< total number of nodes in mesh
   mesh_Ncells, &                                                                                   !< total number of cells in mesh
   mesh_maxNsharedElems                                                                             !< max number of CP elements sharing a node


integer, dimension(:,:), allocatable :: &
   mesh_cellnodeParent                                                                              !< cellnode's parent element ID, cellnode's intra-element ID

 integer,dimension(:,:,:), allocatable :: &
   mesh_cell2, &                                                                                        !< cell connectivity for each element,ip/cell
   mesh_cell                                                                                        !< cell connectivity for each element,ip/cell


 integer, dimension(4), parameter :: FE_NipNeighbors = &                                            !< number of ip neighbors / cell faces in a specific cell type
 int([&
      3, & ! (2D 3node)
      4, & ! (2D 4node)
      4, & ! (3D 4node)
      6  & ! (3D 8node)
  ],pInt)

 integer,     dimension(:),   allocatable :: &
   microstructureAt, &
   homogenizationAt
    
 integer :: &
   mesh_NelemSets
   
 character(len=64), dimension(:), allocatable :: &
   mesh_nameElemSet
 integer, dimension(:,:), allocatable :: &
   mesh_mapElemSet                                                                                  !< list of elements in elementSet
 integer, dimension(:,:), allocatable, target :: &
   mesh_mapFEtoCPelem, &                                                                            !< [sorted FEid, corresponding CPid]
   mesh_mapFEtoCPnode                                                                               !< [sorted FEid, corresponding CPid]

 integer, dimension(:,:,:,:), allocatable :: &
   mesh_ipNeighborhood2                                                                             !< 6 or less neighboring IPs as [element_num, IP_index, neighbor_index that points to me]


 public :: &
   mesh_init, &
   mesh_build_cellnodes, &
   mesh_build_ipCoordinates, &
   mesh_FEasCP
 
 
contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

  integer, intent(in) :: el, ip
   
  integer, parameter  :: FILEUNIT = 222
  character(len=pStringLen), dimension(:), allocatable :: inputFile                                 !< file content, separated per lines

  integer :: j, fileFormatVersion, elemType, &
    mesh_maxNelemInSet, &
    mesh_nElems, &
    hypoelasticTableStyle, &
    initialcondTableStyle
  integer, dimension(:), allocatable :: &
   marc_matNumber                                                                                   !< array of material numbers for hypoelastic material (Marc only)
  logical :: myDebug
   real(pReal), dimension(:,:,:), allocatable:: &
   mesh_ipArea                                                                                      !< area of interface to neighboring IP (initially!)
 real(pReal),dimension(:,:,:,:), allocatable :: &
   mesh_ipAreaNormal                                                                                !< area normal of interface to neighboring IP (initially!)
   
 real(pReal), dimension(:,:), allocatable :: &
   ip_reshaped
 
  write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'
 
  mesh_unitlength = numerics_unitlength                                                              ! set physical extent of a length unit in mesh
  inputFile = IO_read_ASCII(trim(modelName)//trim(InputFileExtension))
  myDebug = (iand(debug_level(debug_mesh),debug_levelBasic) /= 0)
 
  ! parsing Marc input file
  call IO_open_inputFile(FILEUNIT,modelName)
  fileFormatVersion = mesh_marc_get_fileFormat(inputFile)
  call mesh_marc_get_tableStyles(initialcondTableStyle,hypoelasticTableStyle,inputFile)
  if (fileFormatVersion > 12) &
    marc_matNumber = mesh_marc_get_matNumber(FILEUNIT,hypoelasticTableStyle)
  call mesh_marc_count_nodesAndElements(mesh_nNodes, mesh_nElems, FILEUNIT)
  call mesh_marc_count_elementSets(mesh_NelemSets,mesh_maxNelemInSet,FILEUNIT) 
  allocate(mesh_nameElemSet(mesh_NelemSets)); mesh_nameElemSet = 'n/a'
  allocate(mesh_mapElemSet(1+mesh_maxNelemInSet,mesh_NelemSets),source=0)
  call mesh_marc_map_elementSets(mesh_nameElemSet,mesh_mapElemSet,FILEUNIT)
  allocate (mesh_mapFEtoCPelem(2,mesh_nElems), source = 0)
  call mesh_marc_map_elements(hypoelasticTableStyle,mesh_nameElemSet,mesh_mapElemSet,&
                              mesh_nElems,fileFormatVersion,marc_matNumber,FILEUNIT)  
  allocate (mesh_mapFEtoCPnode(2,mesh_Nnodes),source=0)
  call mesh_marc_map_nodes(mesh_Nnodes,inputFile)                                                                  !ToDo: don't work on global variables
  
  mesh_node0 = mesh_marc_build_nodes(mesh_Nnodes,inputFile)
  mesh_node  = mesh_node0
  
  elemType = mesh_marc_getElemType(mesh_nElems,FILEUNIT)
  if (myDebug) write(6,'(a)') ' Counted CP sizes'; flush(6)
  
  call theMesh%init('mesh',elemType,mesh_node0)
  call theMesh%setNelems(mesh_nElems)
  
  allocate(microstructureAt(theMesh%nElems), source=0)
  allocate(homogenizationAt(theMesh%nElems), source=0)
 
  mesh_FEnodes = mesh_marc_buildElements(theMesh%nElems,theMesh%elem%nNodes,FILEUNIT)
  call mesh_marc_buildElements2(microstructureAt,homogenizationAt, &
                               mesh_nElems,theMesh%elem%nNodes,initialcondTableStyle,FILEUNIT)
  if (myDebug) write(6,'(a)') ' Built elements'; flush(6)
  close (FILEUNIT)


#if defined(DAMASK_HDF5)
  call results_openJobFile
  call HDF5_closeGroup(results_addGroup('geometry'))
  call results_writeDataset('geometry',mesh_FEnodes,'C',&
                            'connectivity of the elements','-')
  call results_closeJobFile
#endif
 
  call buildCells(theMesh,theMesh%elem,mesh_FEnodes)
   
  call mesh_build_cellconnectivity
  if (myDebug) write(6,'(a)') ' Built cell connectivity'; flush(6)
  mesh_cellnode = mesh_build_cellnodes()
  if (myDebug) write(6,'(a)') ' Built cell nodes'; flush(6)
 
  allocate(mesh_ipCoordinates(3,theMesh%elem%nIPs,theMesh%nElems),source=0.0_pReal)
  call mesh_build_ipCoordinates
  if (myDebug) write(6,'(a)') ' Built IP coordinates'; flush(6)

  allocate(mesh_ipArea(theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems))
  allocate(mesh_ipAreaNormal(3,theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems))
  call mesh_build_ipAreas(mesh_ipArea,mesh_ipAreaNormal)
  if (myDebug) write(6,'(a)') ' Built IP areas'; flush(6)
 
  call IP_neighborhood2
  if (myDebug) write(6,'(a)') ' Built IP neighborhood'; flush(6)
 
  if (usePingPong .and. (mesh_Nelems /= theMesh%nElems)) &
    call IO_error(600)                                                                          ! ping-pong must be disabled when having non-DAMASK elements
  if (debug_e < 1 .or. debug_e > theMesh%nElems) &
    call IO_error(602,ext_msg='element')                                                        ! selected element does not exist
  if (debug_i < 1 .or. debug_i > theMesh%elem%nIPs) &
    call IO_error(602,ext_msg='IP')                                                             ! selected element does not have requested IP
 
  FEsolving_execElem = [ 1,theMesh%nElems ]                                                      ! parallel loop bounds set to comprise all DAMASK elements
  allocate(FEsolving_execIP(2,theMesh%nElems), source=1)                                    ! parallel loop bounds set to comprise from first IP...
  FEsolving_execIP(2,:) = theMesh%elem%nIPs
 
  allocate(calcMode(theMesh%elem%nIPs,theMesh%nElems))
  calcMode = .false.                                                                                 ! pretend to have collected what first call is asking (F = I)
  calcMode(ip,mesh_FEasCP('elem',el)) = .true.                                                       ! first ip,el needs to be already pingponged to "calc"
 
  ip_reshaped = reshape(mesh_ipCoordinates,[3,theMesh%elem%nIPs*theMesh%nElems])
  call discretization_init(microstructureAt,homogenizationAt,&
                           ip_reshaped,&
                           mesh_node0)
#if defined(DAMASK_HDF5)
  call results_openJobFile
  call results_writeDataset('geometry',ip_reshaped,'x_c', &
                            'cell center coordinates','m')
  call results_writeDataset('geometry',mesh_node0,'x_n', &
                            'nodal coordinates','m')
  call results_closeJobFile()
#endif
  call geometry_plastic_nonlocal_setIPvolume(IPvolume())
  call geometry_plastic_nonlocal_setIPneighborhood(mesh_ipNeighborhood2)
  call geometry_plastic_nonlocal_setIParea(mesh_ipArea)
  call geometry_plastic_nonlocal_setIPareaNormal(mesh_ipAreaNormal)

end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Figures out version of Marc input file format
!--------------------------------------------------------------------------------------------------
integer function mesh_marc_get_fileFormat(fileContent)
 
  character(len=pStringLen), dimension(:), intent(in) :: fileContent                                !< file content, separated per lines
  
  integer, allocatable, dimension(:) :: chunkPos
  integer :: l
  
  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if ( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'version') then
      mesh_marc_get_fileFormat = IO_intValue(fileContent(l),chunkPos,2)
      exit
    endif
  enddo

end function mesh_marc_get_fileFormat


!--------------------------------------------------------------------------------------------------
!> @brief Figures out table styles for initial cond and hypoelastic
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_get_tableStyles(initialcond,hypoelastic,fileContent)
 
  integer, intent(out) :: initialcond, hypoelastic
  character(len=pStringLen), dimension(:), intent(in) :: fileContent                                !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: l
  
  initialcond = 0
  hypoelastic = 0

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if ( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'table' .and. chunkPos(1) > 5) then
      initialcond = IO_intValue(fileContent(l),chunkPos,4)
      hypoelastic = IO_intValue(fileContent(l),chunkPos,5)
      exit
    endif
  enddo

end subroutine mesh_marc_get_tableStyles


!--------------------------------------------------------------------------------------------------
!> @brief Figures out material number of hypoelastic material
!--------------------------------------------------------------------------------------------------
function mesh_marc_get_matNumber(fileUnit,tableStyle)
 
  integer, intent(in)                :: fileUnit, tableStyle
  integer, dimension(:), allocatable :: mesh_marc_get_matNumber

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i, j, data_blocks
  character(len=300) :: line

  data_blocks = 1
  rewind(fileUnit)
  do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)
   
   if ( IO_lc(IO_stringValue(line,chunkPos,1)) == 'hypoelastic') then
     read (fileUnit,'(A300)',END=620) line
     if (len(trim(line))/=0) then
       chunkPos = IO_stringPos(line)
       data_blocks = IO_intValue(line,chunkPos,1)
     endif
     allocate(mesh_marc_get_matNumber(data_blocks), source = 0)
     do i=1,data_blocks                                                                             ! read all data blocks
       read (fileUnit,'(A300)',END=620) line
       chunkPos = IO_stringPos(line)
       mesh_marc_get_matNumber(i) = IO_intValue(line,chunkPos,1)
       do j=1,2 + tableStyle                                                                        ! read 2 or 3 remaining lines of data block
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
 
  integer, intent(in)  :: fileUnit
  integer, intent(out) :: nNodes, nElems
  
  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) :: line

  nNodes = 0
  nElems = 0

  rewind(fileUnit)
  do
    read (fileUnit,'(A300)',END=620) line
    chunkPos = IO_stringPos(line)

    if ( IO_lc(IO_StringValue(line,chunkPos,1)) == 'sizing') &
      nElems = IO_IntValue (line,chunkPos,3)
    if ( IO_lc(IO_StringValue(line,chunkPos,1)) == 'coordinates') then
      read (fileUnit,'(A300)') line
      chunkPos = IO_stringPos(line)
      nNodes = IO_IntValue (line,chunkPos,2)
      exit                                                                                          ! assumes that "coordinates" comes later in file
    endif
  enddo

620 end subroutine mesh_marc_count_nodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of element sets in mesh.
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_count_elementSets(nElemSets,maxNelemInSet,fileUnit)
 
  integer, intent(out) :: nElemSets, maxNelemInSet
  integer, intent(in)  :: fileUnit

  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) :: line

  nElemSets     = 0
  maxNelemInSet = 0

  rewind(fileUnit)
  do
    read (fileUnit,'(A300)',END=620) line
    chunkPos = IO_stringPos(line)

    if ( IO_lc(IO_StringValue(line,chunkPos,1)) == 'define' .and. &
         IO_lc(IO_StringValue(line,chunkPos,2)) == 'element' ) then
      nElemSets = nElemSets + 1
      maxNelemInSet = max(maxNelemInSet, IO_countContinuousIntValues(fileUnit))
    endif
  enddo

620 end subroutine mesh_marc_count_elementSets


!--------------------------------------------------------------------------------------------------
!> @brief map element sets
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_elementSets(nameElemSet,mapElemSet,fileUnit)
 
  character(len=64), dimension(:),   intent(out) :: nameElemSet
  integer,           dimension(:,:), intent(out) :: mapElemSet
  integer,                           intent(in)  :: fileUnit

  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) :: line
  integer :: elemSet
  
  elemSet = 0

  rewind(fileUnit)
  do
    read (fileUnit,'(A300)',END=620) line
    chunkPos = IO_stringPos(line)
    if( (IO_lc(IO_stringValue(line,chunkPos,1)) == 'define' ) .and. &
        (IO_lc(IO_stringValue(line,chunkPos,2)) == 'element' ) ) then
       elemSet = elemSet+1
       nameElemSet(elemSet)  = trim(IO_stringValue(line,chunkPos,4))
       mapElemSet(:,elemSet) = IO_continuousIntValues(fileUnit,size(mapElemSet,1)-1,nameElemSet,mapElemSet,size(nameElemSet))
    endif
  enddo

620 end subroutine mesh_marc_map_elementSets



!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_elements(tableStyle,nameElemSet,mapElemSet,nElems,fileFormatVersion,matNumber,fileUnit)
 
 integer, intent(in) :: fileUnit,tableStyle,nElems,fileFormatVersion
 integer, dimension(:), intent(in) :: matNumber
 character(len=64), intent(in), dimension(:) :: nameElemSet
 integer, dimension(:,:), intent(in) :: &
   mapElemSet

 integer, allocatable, dimension(:) :: chunkPos
 character(len=300) :: line, &
                       tmp

 integer, dimension (1+nElems) :: contInts
 integer :: i,cpElem

 cpElem = 0
 contInts = 0
 rewind(fileUnit)
 do
   read (fileUnit,'(A300)',END=620) line
   chunkPos = IO_stringPos(line)
   if (fileFormatVersion < 13) then                                                                       ! Marc 2016 or earlier
     if( IO_lc(IO_stringValue(line,chunkPos,1)) == 'hypoelastic' ) then
       do i=1,3+TableStyle                                                                          ! skip three (or four if new table style!) lines
         read (fileUnit,'(A300)') line
       enddo
       contInts = IO_continuousIntValues(fileUnit,nElems,nameElemSet,&
                                              mapElemSet,size(nameElemSet))
       exit
     endif  
   else                                                                                             ! Marc2017 and later
     if ( IO_lc(IO_stringValue(line,chunkPos,1)) == 'connectivity') then
       read (fileUnit,'(A300)',END=620) line
       chunkPos = IO_stringPos(line)
       if(any(matNumber==IO_intValue(line,chunkPos,6))) then
         do 
           read (fileUnit,'(A300)',END=620) line
           chunkPos = IO_stringPos(line)
           tmp = IO_lc(IO_stringValue(line,chunkPos,1))
           if (verify(trim(tmp),"0123456789")/=0) then                                              ! found keyword
             exit
           else
             contInts(1) = contInts(1) + 1  
             read (tmp,*) contInts(contInts(1)+1)     
           endif
         enddo
       endif  
     endif
   endif    
 enddo    
620 do i = 1,contInts(1)
      cpElem = cpElem+1
      mesh_mapFEtoCPelem(1,cpElem) = contInts(1+i)
      mesh_mapFEtoCPelem(2,cpElem) = cpElem
    enddo
 
call math_sort(mesh_mapFEtoCPelem)

end subroutine mesh_marc_map_elements


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_map_nodes(nNodes,fileContent)

  integer, intent(in)                                 :: nNodes
  character(len=pStringLen), dimension(:), intent(in) :: fileContent                                !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i, l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'coordinates' ) then
      do i = 1,nNodes
        mesh_mapFEtoCPnode(1,i) = IO_fixedIntValue (fileContent(l+1+i),[0,10],1)
        mesh_mapFEtoCPnode(2,i) = i
      enddo
      exit
    endif
  enddo

  call math_sort(mesh_mapFEtoCPnode)

end subroutine mesh_marc_map_nodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!--------------------------------------------------------------------------------------------------
function mesh_marc_build_nodes(nNode,fileContent) result(nodes)

  integer,                                 intent(in) :: nNode
  character(len=pStringLen), dimension(:), intent(in) :: fileContent                                !< file content, separated per lines
  
  real(pReal), dimension(3,nNode)    :: nodes 
  integer, dimension(5), parameter   :: node_ends = [0,10,30,50,70]
  integer, allocatable, dimension(:) :: chunkPos
  integer :: i,j,m,l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'coordinates' ) then
      do i=1,nNode
        m = mesh_FEasCP('node',IO_fixedIntValue(fileContent(l+1+i),node_ends,1))
        do j = 1,3
          nodes(j,m) = mesh_unitlength * IO_fixedNoEFloatValue(fileContent(l+1+i),node_ends,j+1)
        enddo
      enddo
      exit
    endif
  enddo

end function mesh_marc_build_nodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets element type (and checks if the whole mesh comprises of only one type)
!--------------------------------------------------------------------------------------------------
integer function mesh_marc_getElemType(nElem,fileUnit)

  integer, intent(in) :: &
    nElem, &
    fileUnit

  type(tElement) :: tempEl
  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) :: line
  integer :: i,t

  t = -1

  rewind(fileUnit)
  do
    read (fileUnit,'(A300)',END=620) line
    chunkPos = IO_stringPos(line)
    if( IO_lc(IO_stringValue(line,chunkPos,1)) == 'connectivity' ) then
      read (fileUnit,'(A300)') line                                                                 ! Garbage line
      do i=1,nElem                                                                                  ! read all elements
        read (fileUnit,'(A300)') line
        chunkPos = IO_stringPos(line)
        if (t == -1) then
          t = mapElemtype(IO_stringValue(line,chunkPos,2))
          call tempEl%init(t)
          mesh_marc_getElemType = t
        else
          if (t /= mapElemtype(IO_stringValue(line,chunkPos,2))) call IO_error(191,el=t,ip=i)
        endif
        call IO_skipChunks(fileUnit,tempEl%nNodes-(chunkPos(1)-2))
      enddo
      exit
    endif
  enddo

  contains 

  !--------------------------------------------------------------------------------------------------
  !> @brief mapping of Marc element types to internal representation
  !--------------------------------------------------------------------------------------------------
  integer function mapElemtype(what)
  
   character(len=*), intent(in) :: what
  
   select case (IO_lc(what))
      case (   '6')
        mapElemtype = 1            ! Two-dimensional Plane Strain Triangle
      case ( '155', &
             '125', &
             '128')
        mapElemtype = 2            ! Two-dimensional Plane Strain triangle (155: cubic shape function, 125/128: second order isoparametric)
      case ( '11')
        mapElemtype = 3            ! Arbitrary Quadrilateral Plane-strain
      case ( '27')
        mapElemtype = 4            ! Plane Strain, Eight-node Distorted Quadrilateral
      case ( '54')
        mapElemtype = 5            ! Plane Strain, Eight-node Distorted Quadrilateral with reduced integration
      case ( '134')
        mapElemtype = 6            ! Three-dimensional Four-node Tetrahedron
      case ( '157')
        mapElemtype = 7            ! Three-dimensional, Low-order, Tetrahedron, Herrmann Formulations
      case ( '127')
        mapElemtype = 8            ! Three-dimensional Ten-node Tetrahedron
      case ( '136')
        mapElemtype = 9            ! Three-dimensional Arbitrarily Distorted Pentahedral
      case ( '117', &
             '123')
        mapElemtype = 10           ! Three-dimensional Arbitrarily Distorted linear hexahedral with reduced integration
      case ( '7')
        mapElemtype = 11           ! Three-dimensional Arbitrarily Distorted Brick
      case ( '57')
        mapElemtype = 12           ! Three-dimensional Arbitrarily Distorted quad hexahedral with reduced integration
      case ( '21')
        mapElemtype = 13           ! Three-dimensional Arbitrarily Distorted quadratic hexahedral
      case default
        call IO_error(error_ID=190,ext_msg=IO_lc(what))
   end select

end function mapElemtype


620 end function mesh_marc_getElemType


!--------------------------------------------------------------------------------------------------
!> @brief Stores node IDs
!--------------------------------------------------------------------------------------------------
function mesh_marc_buildElements(nElem,nNodes,fileUnit)
 
  integer, intent(in) :: &
    nElem, &
    nNodes, &                                                                                       !< number of nodes per element
    fileUnit
    
  integer, dimension(nElem,nNodes) :: &
    mesh_marc_buildElements
 
  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) line
 
  integer, dimension(1+nElem) :: contInts
  integer :: i,j,t,sv,myVal,e,nNodesAlreadyRead
 
  rewind(fileUnit)
  do
    read (fileUnit,'(A300)',END=620) line
    chunkPos = IO_stringPos(line)
    if( IO_lc(IO_stringValue(line,chunkPos,1)) == 'connectivity' ) then
      read (fileUnit,'(A300)',END=620) line                                                         ! garbage line
      do i = 1,nElem
        read (fileUnit,'(A300)',END=620) line
        chunkPos = IO_stringPos(line)
        e = mesh_FEasCP('elem',IO_intValue(line,chunkPos,1))
        if (e /= 0) then                                                                            ! disregard non CP elems
          nNodesAlreadyRead = 0
          do j = 1,chunkPos(1)-2
            mesh_marc_buildElements(j,e) = mesh_FEasCP('node',IO_IntValue(line,chunkPos,j+2))                  ! CP ids of nodes
          enddo
          nNodesAlreadyRead = chunkPos(1) - 2
          do while(nNodesAlreadyRead < nNodes)                                                      ! read on if not all nodes in one line
            read (fileUnit,'(A300)',END=620) line
            chunkPos = IO_stringPos(line)
            do j = 1,chunkPos(1)
              mesh_marc_buildElements(nNodesAlreadyRead+j,e) = mesh_FEasCP('node',IO_IntValue(line,chunkPos,j)) ! CP ids of nodes
            enddo
            nNodesAlreadyRead = nNodesAlreadyRead + chunkPos(1)
          enddo
        endif
      enddo
      exit
    endif
  enddo
620 rewind(fileUnit)
 
end function mesh_marc_buildElements


!--------------------------------------------------------------------------------------------------
!> @brief Stores homogenization and microstructure ID
!--------------------------------------------------------------------------------------------------
subroutine mesh_marc_buildElements2(microstructureAt,homogenizationAt, &
                                   nElem,nNodes,initialcondTableStyle,fileUnit)

  integer, dimension(:), intent(out) :: &
    microstructureAt, &
    homogenizationAt
  integer, intent(in) :: &
    nElem, &
    nNodes, &                                                                                       !< number of nodes per element
    initialcondTableStyle, &
    fileUnit

  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) line
 
  integer, dimension(1+nElem) :: contInts
  integer :: i,j,t,sv,myVal,e,nNodesAlreadyRead

  read (fileUnit,'(A300)',END=630) line
  do
    chunkPos = IO_stringPos(line)
    if( (IO_lc(IO_stringValue(line,chunkPos,1)) == 'initial') .and. &
        (IO_lc(IO_stringValue(line,chunkPos,2)) == 'state') ) then
      if (initialcondTableStyle == 2) read (fileUnit,'(A300)',END=630) line                         ! read extra line for new style
      read (fileUnit,'(A300)',END=630) line                                                         ! read line with index of state var
      chunkPos = IO_stringPos(line)
      sv = IO_IntValue(line,chunkPos,1)                                                             ! figure state variable index
      if( (sv == 2).or.(sv == 3) ) then                                                             ! only state vars 2 and 3 of interest
        read (fileUnit,'(A300)',END=630) line                                                       ! read line with value of state var
        chunkPos = IO_stringPos(line)
        do while (scan(IO_stringValue(line,chunkPos,1),'+-',back=.true.)>1)                         ! is noEfloat value?
          myVal = nint(IO_fixedNoEFloatValue(line,[0,20],1),pInt)                                   ! state var's value
          if (initialcondTableStyle == 2) then
            read (fileUnit,'(A300)',END=630) line                                                   ! read extra line
            read (fileUnit,'(A300)',END=630) line                                                   ! read extra line
          endif
          contInts = IO_continuousIntValues&                                                        ! get affected elements
                    (fileUnit,theMesh%nElems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
          do i = 1,contInts(1)
            e = mesh_FEasCP('elem',contInts(1+i))
            if (sv == 2) microstructureAt(e) = myVal
            if (sv == 3) homogenizationAt(e) = myVal
          enddo
          if (initialcondTableStyle == 0) read (fileUnit,'(A300)',END=630) line                     ! ignore IP range for old table style
          read (fileUnit,'(A300)',END=630) line
          chunkPos = IO_stringPos(line)
        enddo
      endif
    else
      read (fileUnit,'(A300)',END=630) line
    endif
  enddo
 
630 end subroutine mesh_marc_buildElements2


subroutine buildCells(thisMesh,elem,connectivity_elem)

  class(tMesh)                       :: thisMesh
  type(tElement),         intent(in) :: elem
  integer,dimension(:,:), intent(in) :: connectivity_elem
  
  integer,dimension(:),     allocatable :: candidates_local
  integer,dimension(:,:),   allocatable :: parentsAndWeights,candidates_global, connectivity_cell_reshape
  integer,dimension(:,:,:), allocatable :: connectivity_cell
  
  real(pReal), dimension(:,:), allocatable :: nodes_new,nodes
  integer :: e, n, c, p, s,i,m,j,nParentNodes,nCellNode,Nelem,candidateID
  
  Nelem = thisMesh%Nelems

!---------------------------------------------------------------------------------------------------
! initialize global connectivity to negative local connectivity
  allocate(connectivity_cell(elem%NcellNodesPerCell,elem%nIPs,Nelem))
  connectivity_cell = -spread(elem%cell,3,Nelem)                                                    ! local cell node ID
  
!---------------------------------------------------------------------------------------------------
! set connectivity of cell nodes that coincide with FE nodes (defined by 1 parent node)
! and renumber local (negative) to global (positive) node ID
  do e = 1, Nelem
    do c = 1, elem%NcellNodes
      realNode: if (count(elem%cellNodeParentNodeWeights(:,c) /= 0) == 1) then
        where(connectivity_cell(:,:,e) == -c)
          connectivity_cell(:,:,e) = connectivity_elem(c,e)
        end where
      endif realNode
    enddo
  enddo
  nCellNode = thisMesh%nNodes
  

!---------------------------------------------------------------------------------------------------
! set connectivity of cell nodes that are defined by 2,...,nNodes real nodes
  do nParentNodes = 2, elem%nNodes
  
    ! get IDs of local cell nodes that are defined by the current number of parent nodes
    candidates_local = [integer::]
    do c = 1, elem%NcellNodes
      if (count(elem%cellNodeParentNodeWeights(:,c) /= 0) == nParentNodes) &
        candidates_local = [candidates_local,c]
    enddo
    s = size(candidates_local)
    
    if (allocated(candidates_global)) deallocate(candidates_global)
    allocate(candidates_global(nParentNodes*2+2,s*Nelem))                                           ! stores parent node ID + weight together with element ID and cellnode id (local)
    parentsAndWeights = reshape([(0, i = 1,2*nParentNodes)],[nParentNodes,2])                       ! (re)allocate
    
    do e = 1, Nelem
      do i = 1, size(candidates_local)
        candidateID = (e-1)*size(candidates_local)+i                                                ! including duplicates, runs to (Nelem*size(candidates_local))
        c = candidates_local(i)                                                                     ! c is local cellnode ID for connectivity
        p = 0
        do j = 1, size(elem%cellNodeParentNodeWeights(:,c))
          if (elem%cellNodeParentNodeWeights(j,c) /= 0) then                                        ! real node 'j' partly defines cell node 'c'
            p = p + 1                                                               
            parentsAndWeights(p,1:2) = [connectivity_elem(j,e),elem%cellNodeParentNodeWeights(j,c)]
          endif
        enddo
        ! store (and order) real node IDs and their weights together with the element number and local ID
        do p = 1, nParentNodes
          m = maxloc(parentsAndWeights(:,1),1)
          
          candidates_global(p,                                candidateID) = parentsAndWeights(m,1)
          candidates_global(p+nParentNodes,                   candidateID) = parentsAndWeights(m,2)
          candidates_global(nParentNodes*2+1:nParentNodes*2+2,candidateID) = [e,c]
          
          parentsAndWeights(m,1) = -huge(parentsAndWeights(m,1))                                    ! out of the competition
        enddo
      enddo
    enddo

    ! sort according to real node IDs + weight (from left to right)
    call math_sort(candidates_global,sortDim=1)                                                     ! sort according to first column
    
    do p = 2, nParentNodes*2
      n = 1
      do while(n <= size(candidates_local)*Nelem)
        j=0
        do while (n+j<= size(candidates_local)*Nelem)
          if (candidates_global(p-1,n+j)/=candidates_global(p-1,n)) exit
          j = j + 1
        enddo
        e = n+j-1
        if (any(candidates_global(p,n:e)/=candidates_global(p,n))) &
          call math_sort(candidates_global(:,n:e),sortDim=p)
        n = e+1
      enddo
    enddo
    
    i = uniqueRows(candidates_global(1:2*nParentNodes,:))
      
    
    ! calculate coordinates of cell nodes and insert their ID into the cell conectivity
    nodes_new = reshape([(0.0_pReal,j = 1, 3*i)], [3,i])
    
    i = 1
    n = 1
    do while(n <= size(candidates_local)*Nelem)
     j=0
     parentsAndWeights(:,1) = candidates_global(1:nParentNodes,n+j)
     parentsAndWeights(:,2) = candidates_global(nParentNodes+1:nParentNodes*2,n+j)
     e = candidates_global(nParentNodes*2+1,n+j)
     c = candidates_global(nParentNodes*2+2,n+j)
     do m = 1, nParentNodes
       nodes_new(:,i) = nodes_new(:,i) &
                      + thisMesh%node_0(:,parentsAndWeights(m,1)) * real(parentsAndWeights(m,2),pReal)
     enddo
     nodes_new(:,i) =  nodes_new(:,i)/real(sum(parentsAndWeights(:,2)),pReal)

     do while (n+j<= size(candidates_local)*Nelem)
       if (any(candidates_global(1:2*nParentNodes,n+j)/=candidates_global(1:2*nParentNodes,n))) exit
       where (connectivity_cell(:,:,candidates_global(nParentNodes*2+1,n+j)) == -candidates_global(nParentNodes*2+2,n+j)) ! still locally defined
         connectivity_cell(:,:,candidates_global(nParentNodes*2+1,n+j)) = nCellNode + i                                   ! gets current new cell node id
       end where
       
        j = j + 1
     enddo
     i=i+1
     n = n+j

    enddo
    nCellNode = nCellNode + i
    if (i/=0) nodes = reshape([nodes,nodes_new],[3,nCellNode])
  enddo
  thisMesh%node_0 = nodes
  mesh_cell2 = connectivity_cell
 
#if defined(DAMASK_HDF5) 
  connectivity_cell_reshape = reshape(connectivity_cell,[elem%NcellNodesPerCell,elem%nIPs*Nelem])
  call results_openJobFile
  call results_writeDataset('geometry',connectivity_cell_reshape,'c',&
                            'connectivity of the cells','-')
  call results_closeJobFile
#endif  

  contains

  !------------------------------------------------------------------------------------------------
  !> @brief count unique rows (same rows need to be stored consecutively)
  !------------------------------------------------------------------------------------------------
  pure function uniqueRows(A) result(u)
    
    integer, dimension(:,:), intent(in) :: A                                                        !< array, rows need to be sorted

    integer :: &
      u, &                                                                                          !< # of unique rows
      r, &                                                                                          !< row counter
      d                                                                                             !< duplicate counter

    u = 0
    r = 1
    do while(r <= size(A,2))
      d = 0
      do while (r+d<= size(A,2))
        if (any(A(:,r)/=A(:,r+d))) exit
        d = d+1
      enddo
      u = u+1
      r = r+d
    enddo
      
  end function uniqueRows

end subroutine buildCells


!--------------------------------------------------------------------------------------------------
!> @brief Split CP elements into cells.
!> @details Build a mapping between cells and the corresponding cell nodes ('mesh_cell').
!> Cell nodes that are also matching nodes are unique in the list of cell nodes,
!> all others (currently) might be stored more than once.
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_cellconnectivity

 integer, dimension(:), allocatable :: &
   matchingNode2cellnode
   
 integer, dimension(10), parameter :: FE_NmatchingNodes = &               !< number of nodes that are needed for face matching in a specific type of element geometry
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

 integer, dimension(:,:), allocatable :: &
   cellnodeParent
 integer, dimension(theMesh%elem%Ncellnodes) :: &
   localCellnode2globalCellnode
 integer :: &
   e,n,i, &
   matchingNodeID, &
   localCellnodeID

 allocate(mesh_cell(theMesh%elem%NcellNodesPerCell,theMesh%elem%nIPs,theMesh%nElems), source=0)
 allocate(matchingNode2cellnode(theMesh%nNodes),                            source=0)
 allocate(cellnodeParent(2,theMesh%elem%Ncellnodes*theMesh%nElems),       source=0)
 
 mesh_Ncells = theMesh%nElems*theMesh%elem%nIPs
!--------------------------------------------------------------------------------------------------
! Count cell nodes (including duplicates) and generate cell connectivity list
 mesh_Ncellnodes = 0

 do e = 1,theMesh%nElems
   localCellnode2globalCellnode = 0
   do i = 1,theMesh%elem%nIPs
     do n = 1,theMesh%elem%NcellnodesPerCell
       localCellnodeID = theMesh%elem%cell(n,i)
       if (localCellnodeID <= FE_NmatchingNodes(theMesh%elem%geomType)) then                                            ! this cell node is a matching node
         matchingNodeID = mesh_FEnodes(localCellnodeID,e)
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
function mesh_build_cellnodes()

 
 real(pReal), dimension(3,mesh_Ncellnodes) :: mesh_build_cellnodes

 integer :: &
   e,n,m, &
   localCellnodeID
 real(pReal), dimension(3) :: &
   myCoords

 mesh_build_cellnodes = 0.0_pReal
 do n = 1,mesh_Ncellnodes                                                                           ! loop over cell nodes
   e = mesh_cellnodeParent(1,n)
   localCellnodeID = mesh_cellnodeParent(2,n)
   myCoords = 0.0_pReal
   do m = 1,theMesh%elem%nNodes
     myCoords = myCoords + mesh_node(1:3,mesh_FEnodes(m,e)) &
                         * theMesh%elem%cellNodeParentNodeWeights(m,localCellnodeID)
   enddo
   mesh_build_cellnodes(1:3,n) = myCoords / sum(theMesh%elem%cellNodeParentNodeWeights(:,localCellnodeID))
 enddo

end function mesh_build_cellnodes


!---------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume.
!> @details The IP volume is calculated differently depending on the cell type.
!> 2D cells assume an element depth of one in order to calculate the volume.
!> For the hexahedral cell we subdivide the cell into subvolumes of pyramidal
!> shape with a cell face as basis and the central ip at the tip. This subvolume is
!> calculated as an average of four tetrahedals with three corners on the cell face
!> and one corner at the central ip.
!---------------------------------------------------------------------------------------------------
function IPvolume()
 
  real(pReal), dimension(theMesh%elem%nIPs,theMesh%nElems) :: IPvolume

  integer :: e,i,c,m,f,n
  real(pReal), dimension(size(theMesh%elem%cellFace,1),size(theMesh%elem%cellFace,2)) :: subvolume

  c = theMesh%elem%cellType
  m = size(theMesh%elem%cellFace,1)

  do e = 1,theMesh%nElems
    select case (c)

      case (1)                                                                                ! 2D 3node
        forall (i = 1:theMesh%elem%nIPs) &                                                           ! loop over ips=cells in this element
          IPvolume(i,e) = math_areaTriangle(theMesh%node_0(1:3,mesh_cell2(1,i,e)), &
                                            theMesh%node_0(1:3,mesh_cell2(2,i,e)), &
                                            theMesh%node_0(1:3,mesh_cell2(3,i,e)))

      case (2)                                                                                ! 2D 4node
        forall (i = 1:theMesh%elem%nIPs) &                                                           ! loop over ips=cells in this element
          IPvolume(i,e) = math_areaTriangle(theMesh%node_0(1:3,mesh_cell2(1,i,e)), &            ! here we assume a planar shape, so division in two triangles suffices
                                            theMesh%node_0(1:3,mesh_cell2(2,i,e)), &
                                            theMesh%node_0(1:3,mesh_cell2(3,i,e))) &
                        + math_areaTriangle(theMesh%node_0(1:3,mesh_cell2(3,i,e)), &
                                            theMesh%node_0(1:3,mesh_cell2(4,i,e)), &
                                            theMesh%node_0(1:3,mesh_cell2(1,i,e)))

      case (3)                                                                                ! 3D 4node
        forall (i = 1:theMesh%elem%nIPs) &                                                           ! loop over ips=cells in this element
          IPvolume(i,e) = math_volTetrahedron(theMesh%node_0(1:3,mesh_cell2(1,i,e)), &
                                              theMesh%node_0(1:3,mesh_cell2(2,i,e)), &
                                              theMesh%node_0(1:3,mesh_cell2(3,i,e)), &
                                              theMesh%node_0(1:3,mesh_cell2(4,i,e)))

      case (4)                                                                                ! 3D 8node
        do i = 1,theMesh%elem%nIPs                                                                   ! loop over ips=cells in this element
          subvolume = 0.0_pReal
          forall(f = 1:FE_NipNeighbors(c), n = 1:m) &
            subvolume(n,f) = math_volTetrahedron(&
                               mesh_cellnode(1:3,mesh_cell(theMesh%elem%cellface(      n     ,f),i,e)), &
                               mesh_cellnode(1:3,mesh_cell(theMesh%elem%cellface(1+mod(n  ,m),f),i,e)), &
                               mesh_cellnode(1:3,mesh_cell(theMesh%elem%cellface(1+mod(n+1,m),f),i,e)), &
                               mesh_ipCoordinates(1:3,i,e))
          IPvolume(i,e) = 0.5_pReal * sum(subvolume)                                         ! each subvolume is based on four tetrahedrons, altough the face consists of only two triangles -> averaging factor two
        enddo

    end select
  enddo

end function IPvolume


!---------------------------------------------------------------------------------------------------
!> @brief cell neighborhood
!---------------------------------------------------------------------------------------------------
subroutine IP_neighborhood2

  integer, dimension(:,:), allocatable :: faces
  integer, dimension(:),   allocatable :: face
  integer :: e,i,f,c,m,n,j,k,l,p, current, next,i2,e2,n2,k2
  logical :: match

  allocate(faces(size(theMesh%elem%cellface,1)+3,size(theMesh%elem%cellface,2)*theMesh%elem%nIPs*theMesh%Nelems))

  ! store cell face definitions
  f = 0
  do e = 1,theMesh%nElems
    do i = 1,theMesh%elem%nIPs
      do n = 1, theMesh%elem%nIPneighbors
        f = f + 1
        face = mesh_cell2(theMesh%elem%cellFace(:,n),i,e)
        storeSorted: do j = 1, size(face)
          faces(j,f) = maxval(face)
          face(maxloc(face)) = -huge(1)
        enddo storeSorted
        faces(j:j+2,f) = [e,i,n]
      enddo
    enddo
  enddo
  
  ! sort ..
  call math_sort(faces,sortDim=1)
  do p = 2, size(faces,1)-2
     n = 1
     do while(n <= size(faces,2))
      j=0
      do while (n+j<= size(faces,2))
        if (faces(p-1,n+j)/=faces(p-1,n)) exit
        j = j + 1
      enddo
      e = n+j-1
      if (any(faces(p,n:e)/=faces(p,n))) call math_sort(faces(:,n:e),sortDim=p)
      n = e+1
     enddo
  enddo

 allocate(mesh_ipNeighborhood2(3,theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems),source=0)
 
 ! find IP neighbors
 f = 1
 do while(f <= size(faces,2))
   e = faces(size(theMesh%elem%cellface,1)+1,f)
   i = faces(size(theMesh%elem%cellface,1)+2,f)
   n = faces(size(theMesh%elem%cellface,1)+3,f)
   
   if (f < size(faces,2)) then
     match = all(faces(1:size(theMesh%elem%cellface,1),f) == faces(1:size(theMesh%elem%cellface,1),f+1))
     e2 = faces(size(theMesh%elem%cellface,1)+1,f+1)
     i2 = faces(size(theMesh%elem%cellface,1)+2,f+1)
     n2 = faces(size(theMesh%elem%cellface,1)+3,f+1)
   else
     match = .false.
   endif
   
   if (match) then
     if (e == e2) then ! same element. MD: I don't think that we need this (not even for other elements)
       k  = theMesh%elem%IPneighbor(n, i)
       k2 = theMesh%elem%IPneighbor(n2,i2)
     endif
     mesh_ipNeighborhood2(1:3,n, i, e)  = [e2,i2,n2]
     mesh_ipNeighborhood2(1:3,n2,i2,e2) = [e, i, n]
     f = f +1
   endif
   f = f +1 
 enddo

end subroutine IP_neighborhood2


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP Coordinates.
! Marc however only provides nodal displacements,
! so in this case the ip coordinates are always calculated on the basis of this subroutine.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR THE MOMENT THIS SUBROUTINE ACTUALLY CALCULATES THE CELL CENTER AND NOT THE IP COORDINATES,
! AS THE IP IS NOT (ALWAYS) LOCATED IN THE CENTER OF THE IP VOLUME.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipCoordinates

  integer :: e,i,n
  real(pReal), dimension(3) :: myCoords


  do e = 1,theMesh%nElems
    do i = 1,theMesh%elem%nIPs
      myCoords = 0.0_pReal
      do n = 1,theMesh%elem%nCellnodesPerCell
        myCoords = myCoords + mesh_cellnode(1:3,mesh_cell2(n,i,e))
      enddo
      mesh_ipCoordinates(1:3,i,e) = myCoords / real(theMesh%elem%nCellnodesPerCell,pReal)
    enddo
  enddo

end subroutine mesh_build_ipCoordinates


!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas, allocate globals '_ipArea', and '_ipAreaNormal'
!--------------------------------------------------------------------------------------------------
subroutine mesh_build_ipAreas(ipArea,ipAreaNormal)

 integer :: e,c,i,f,n,m
 real(pReal), dimension(theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems), intent(out) ::  ipArea
 real(pReal), dimension(3,theMesh%elem%nIPneighbors,theMesh%elem%nIPs,theMesh%nElems), intent(out) :: ipAreaNormal
 real(pReal), dimension (3,size(theMesh%elem%cellFace,2)) :: nodePos, normals
 real(pReal), dimension(3) :: normal
 
 c = theMesh%elem%cellType


   do e = 1,theMesh%nElems                                                                      ! loop over cpElems
     select case (c)

       case (1,2)                                                                         ! 2D 3 or 4 node
         do i = 1,theMesh%elem%nIPs
           do f = 1,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1: size(theMesh%elem%cellface,1)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(theMesh%elem%cellface(n,f),i,e))
             normal(1) =   nodePos(2,2) - nodePos(2,1)                                              ! x_normal =  y_connectingVector
             normal(2) = -(nodePos(1,2) - nodePos(1,1))                                             ! y_normal = -x_connectingVector
             normal(3) = 0.0_pReal
             ipArea(f,i,e) = norm2(normal)
             ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (3)                                                                                ! 3D 4node
         do i = 1,theMesh%elem%nIPs
           do f = 1,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1: size(theMesh%elem%cellface,1)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(theMesh%elem%cellface(n,f),i,e))
             normal = math_cross(nodePos(1:3,2) - nodePos(1:3,1), &
                                 nodePos(1:3,3) - nodePos(1:3,1))
             ipArea(f,i,e) = norm2(normal)
             ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)                             ! ensure unit length of area normal
           enddo
         enddo

       case (4)                                                                                ! 3D 8node
         ! for this cell type we get the normal of the quadrilateral face as an average of
         ! four normals of triangular subfaces; since the face consists only of two triangles,
         ! the sum has to be divided by two; this whole prcedure tries to compensate for
         ! probable non-planar cell surfaces
         m = size(theMesh%elem%cellFace,1)
         do i = 1,theMesh%elem%nIPs
           do f = 1,FE_NipNeighbors(c)                                                         ! loop over cell faces
             forall(n = 1: size(theMesh%elem%cellface,1)) &
               nodePos(1:3,n) = mesh_cellnode(1:3,mesh_cell(theMesh%elem%cellface(n,f),i,e))
             forall(n = 1: size(theMesh%elem%cellface,1)) &
               normals(1:3,n) = 0.5_pReal &
                              * math_cross(nodePos(1:3,1+mod(n  ,m)) - nodePos(1:3,n), &
                                           nodePos(1:3,1+mod(n+1,m)) - nodePos(1:3,n))
             normal = 0.5_pReal * sum(normals,2)
             ipArea(f,i,e) = norm2(normal)
             ipAreaNormal(1:3,f,i,e) = normal / norm2(normal)
           enddo
         enddo

     end select
   enddo
   
end subroutine mesh_build_ipAreas


!--------------------------------------------------------------------------------------------------
!> @brief Gives the FE to CP ID mapping by binary search through lookup array
!! valid questions (what) are 'elem', 'node'
!--------------------------------------------------------------------------------------------------
integer function mesh_FEasCP(what,myID)

  character(len=*), intent(in) :: what
  integer,          intent(in) :: myID
 
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
 
  if (lookupMap(1,lower) == myID) then                                                               ! check at bounds QUESTION is it valid to extend bounds by 1 and just do binary search w/o init check at bounds?
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
