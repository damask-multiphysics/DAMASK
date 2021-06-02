!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Sets up the mesh for the solver MSC.Marc
!--------------------------------------------------------------------------------------------------
module discretization_marc
  use IO
  use prec
  use math
  use DAMASK_interface
  use IO
  use config
  use element
  use discretization
  use geometry_plastic_nonlocal
  use results

  implicit none
  private

  real(pReal),                         public, protected :: &
    mesh_unitlength                                                                                 !< physical length of one unit in mesh MD: needs systematic_name

  integer,  dimension(:), allocatable, public, protected :: &
    discretization_Marc_FEM2DAMASK_elem, &                                                          !< DAMASK element ID for Marc element ID
    discretization_Marc_FEM2DAMASK_node                                                             !< DAMASK node ID for Marc node ID


  type tCellNodeDefinition
    integer, dimension(:,:), allocatable :: parents
    integer, dimension(:,:), allocatable :: weights
  end type tCellNodeDefinition

  type(tCellNodeDefinition), dimension(:), allocatable :: cellNodeDefinition

  integer, dimension(:,:,:), allocatable :: &
    connectivity_cell                                                                               !< cell connectivity for each element,ip/cell

  public :: &
    discretization_Marc_init, &
    discretization_Marc_updateNodeAndIpCoords, &
    discretization_Marc_FEM2DAMASK_cell

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine discretization_Marc_init

  real(pReal), dimension(:,:),     allocatable :: &
   node0_elem, &                                                                                    !< node x,y,z coordinates (initially!)
   node0_cell
  type(tElement) :: elem

  integer,     dimension(:),       allocatable :: &
    materialAt
  integer:: &
    Nelems, &                                                                                       !< total number of elements in the mesh
    debug_e, debug_i

  real(pReal), dimension(:,:),     allocatable :: &
    IP_reshaped
  integer,     dimension(:,:),     allocatable :: &
    connectivity_elem
  real(pReal), dimension(:,:,:,:), allocatable :: &
    unscaledNormals

  class(tNode), pointer :: &
    num_commercialFEM

  print'(/,a)', ' <<<+-  discretization_marc init  -+>>>'; flush(6)

!---------------------------------------------------------------------------------
! read debug parameters
  debug_e = config_debug%get_asInt('element',defaultVal=1)
  debug_i = config_debug%get_asInt('integrationpoint',defaultVal=1)

!--------------------------------------------------------------------------------
! read numerics parameter and do sanity check
  num_commercialFEM => config_numerics%get('commercialFEM',defaultVal = emptyDict)
  mesh_unitlength = num_commercialFEM%get_asFloat('unitlength',defaultVal=1.0_pReal)                ! set physical extent of a length unit in mesh
  if (mesh_unitlength <= 0.0_pReal) call IO_error(301,ext_msg='unitlength')

  call inputRead(elem,node0_elem,connectivity_elem,materialAt)
  nElems = size(connectivity_elem,2)

  if (debug_e < 1 .or. debug_e > nElems)    call IO_error(602,ext_msg='element')
  if (debug_i < 1 .or. debug_i > elem%nIPs) call IO_error(602,ext_msg='IP')

  allocate(cellNodeDefinition(elem%nNodes-1))
  allocate(connectivity_cell(elem%NcellNodesPerCell,elem%nIPs,nElems))

  call buildCells(connectivity_cell,cellNodeDefinition,&
                  elem,connectivity_elem)
  node0_cell  = buildCellNodes(node0_elem)

  IP_reshaped = buildIPcoordinates(node0_cell)

  call discretization_init(materialAt, IP_reshaped, node0_cell)

  call writeGeometry(elem,connectivity_elem,&
                     reshape(connectivity_cell,[elem%NcellNodesPerCell,elem%nIPs*nElems]),&
                     node0_cell,IP_reshaped)

!--------------------------------------------------------------------------------------------------
! geometry information required by the nonlocal CP model
  call geometry_plastic_nonlocal_setIPvolume(IPvolume(elem,node0_cell))
  unscaledNormals = IPareaNormal(elem,nElems,node0_cell)
  call geometry_plastic_nonlocal_setIParea(norm2(unscaledNormals,1))
  call geometry_plastic_nonlocal_setIPareaNormal(unscaledNormals/spread(norm2(unscaledNormals,1),1,3))
  call geometry_plastic_nonlocal_setIPneighborhood(IPneighborhood(elem))
  call geometry_plastic_nonlocal_results

end subroutine discretization_Marc_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate and set current nodal and IP positions (including cell nodes)
!--------------------------------------------------------------------------------------------------
subroutine discretization_Marc_updateNodeAndIpCoords(d_n)

  real(pReal), dimension(:,:), intent(in)  :: d_n

  real(pReal), dimension(:,:), allocatable :: node_cell


  node_cell = buildCellNodes(discretization_NodeCoords0(1:3,1:maxval(discretization_Marc_FEM2DAMASK_node)) + d_n)

  call discretization_setNodeCoords(node_cell)
  call discretization_setIPcoords(buildIPcoordinates(node_cell))

end subroutine discretization_Marc_updateNodeAndIpCoords


!--------------------------------------------------------------------------------------------------
!> @brief Calculate and set current nodal and IP positions (including cell nodes)
!--------------------------------------------------------------------------------------------------
function discretization_marc_FEM2DAMASK_cell(IP_FEM,elem_FEM) result(cell)

  integer, intent(in) :: IP_FEM, elem_FEM
  integer :: cell

  real(pReal), dimension(:,:), allocatable :: node_cell


  cell = (discretization_Marc_FEM2DAMASK_elem(elem_FEM)-1)*discretization_nIPs + IP_FEM


end function discretization_marc_FEM2DAMASK_cell


!--------------------------------------------------------------------------------------------------
!> @brief Write all information needed for the DADF5 geometry
!--------------------------------------------------------------------------------------------------
subroutine writeGeometry(elem, &
                         connectivity_elem,connectivity_cell_reshaped, &
                         coordinates_nodes,coordinates_points)

  type(tElement),              intent(in) :: &
    elem
  integer, dimension(:,:),     intent(in) :: &
    connectivity_elem, &
    connectivity_cell_reshaped
  real(pReal), dimension(:,:), intent(in) :: &
    coordinates_nodes, &
    coordinates_points


  call results_openJobFile
  call results_closeGroup(results_addGroup('geometry'))

  call results_writeDataset(connectivity_elem,'geometry','T_e',&
                            'connectivity of the elements','-')

  call results_writeDataset(connectivity_cell_reshaped,'geometry','T_c', &
                            'connectivity of the cells','-')
  call results_addAttribute('VTK_TYPE',elem%vtkType,'geometry/T_c')

  call results_writeDataset(coordinates_nodes,'geometry','x_n', &
                            'initial coordinates of the nodes','m')

  call results_writeDataset(coordinates_points,'geometry','x_p', &
                            'initial coordinates of the materialpoints (cell centers)','m')

  call results_closeJobFile

end subroutine writeGeometry


!--------------------------------------------------------------------------------------------------
!> @brief Read mesh from marc input file
!--------------------------------------------------------------------------------------------------
subroutine inputRead(elem,node0_elem,connectivity_elem,materialAt)

  type(tElement), intent(out) :: elem
  real(pReal), dimension(:,:), allocatable, intent(out) :: &
    node0_elem                                                                                      !< node x,y,z coordinates (initially!)
  integer, dimension(:,:),     allocatable, intent(out) :: &
    connectivity_elem
  integer,     dimension(:),   allocatable, intent(out) :: &
    materialAt

  integer :: &
    fileFormatVersion, &
    hypoelasticTableStyle, &
    initialcondTableStyle, &
    nNodes, &
    nElems
  integer, dimension(:), allocatable :: &
    matNumber                                                                                       !< material numbers for hypoelastic material
  character(len=pStringLen), dimension(:), allocatable :: &
    inputFile, &                                                                                    !< file content, separated per lines
    nameElemSet
  integer, dimension(:,:), allocatable :: &
    mapElemSet                                                                                      !< list of elements in elementSet


  inputFile = IO_readlines(trim(getSolverJobName())//trim(InputFileExtension))
  call inputRead_fileFormat(fileFormatVersion, &
                            inputFile)
  call inputRead_tableStyles(initialcondTableStyle,hypoelasticTableStyle, &
                             inputFile)
  if (fileFormatVersion > 12) &
    call inputRead_matNumber(matNumber, &
                             hypoelasticTableStyle,inputFile)
  call inputRead_NnodesAndElements(nNodes,nElems,&
                                   inputFile)


  call inputRead_mapElemSets(nameElemSet,mapElemSet,&
                             inputFile)

  call inputRead_elemType(elem, &
                          nElems,inputFile)

  call inputRead_mapElems(discretization_Marc_FEM2DAMASK_elem,&
                          nElems,elem%nNodes,inputFile)

  call inputRead_mapNodes(discretization_Marc_FEM2DAMASK_node,&
                          nNodes,inputFile)

  call inputRead_elemNodes(node0_elem, &
                           Nnodes,inputFile)

  connectivity_elem = inputRead_connectivityElem(nElems,elem%nNodes,inputFile)

  call inputRead_material(materialAt, &
                          nElems,elem%nNodes,nameElemSet,mapElemSet,&
                          initialcondTableStyle,inputFile)
end subroutine inputRead



!--------------------------------------------------------------------------------------------------
!> @brief Figures out version of Marc input file format
!--------------------------------------------------------------------------------------------------
subroutine inputRead_fileFormat(fileFormat,fileContent)

  integer,                        intent(out) :: fileFormat
  character(len=*), dimension(:), intent(in)  :: fileContent                                        !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 2) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'version') then
      fileFormat = IO_intValue(fileContent(l),chunkPos,2)
      exit
    endif
  enddo

end subroutine inputRead_fileFormat


!--------------------------------------------------------------------------------------------------
!> @brief Figures out table styles for initial cond and hypoelastic
!--------------------------------------------------------------------------------------------------
subroutine inputRead_tableStyles(initialcond,hypoelastic,fileContent)

  integer,                        intent(out) :: initialcond, hypoelastic
  character(len=*), dimension(:), intent(in)  :: fileContent                                        !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: l

  initialcond = 0
  hypoelastic = 0

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 6) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'table') then
      initialcond = IO_intValue(fileContent(l),chunkPos,4)
      hypoelastic = IO_intValue(fileContent(l),chunkPos,5)
      exit
    endif
  enddo

end subroutine inputRead_tableStyles


!--------------------------------------------------------------------------------------------------
!> @brief Figures out material number of hypoelastic material
!--------------------------------------------------------------------------------------------------
subroutine inputRead_matNumber(matNumber, &
                               tableStyle,fileContent)

  integer, allocatable, dimension(:), intent(out) :: matNumber
  integer,                            intent(in)  :: tableStyle
  character(len=*),     dimension(:), intent(in)  :: fileContent                                    !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i, j, data_blocks, l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'hypoelastic') then
      if (len_trim(fileContent(l+1))/=0) then
        chunkPos = IO_stringPos(fileContent(l+1))
        data_blocks = IO_intValue(fileContent(l+1),chunkPos,1)
      else
        data_blocks = 1
      endif
      allocate(matNumber(data_blocks), source = 0)
      do i = 0, data_blocks - 1
        j = i*(2+tableStyle) + 1
        chunkPos = IO_stringPos(fileContent(l+1+j))
        matNumber(i+1) = IO_intValue(fileContent(l+1+j),chunkPos,1)
      enddo
      exit
    endif
  enddo

end subroutine inputRead_matNumber


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of nodes and elements
!--------------------------------------------------------------------------------------------------
subroutine inputRead_NnodesAndElements(nNodes,nElems,&
                                       fileContent)

  integer,                        intent(out) :: nNodes, nElems
  character(len=*), dimension(:), intent(in)  :: fileContent                                        !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: l

  nNodes = 0
  nElems = 0

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if    (IO_lc(IO_StringValue(fileContent(l),chunkPos,1)) == 'sizing') then
      nElems = IO_IntValue (fileContent(l),chunkPos,3)
    elseif(IO_lc(IO_StringValue(fileContent(l),chunkPos,1)) == 'coordinates') then
      chunkPos = IO_stringPos(fileContent(l+1))
      nNodes = IO_IntValue (fileContent(l+1),chunkPos,2)
    endif
  enddo

end subroutine inputRead_NnodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of element sets in mesh.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_NelemSets(nElemSets,maxNelemInSet,&
                               fileContent)

  integer,                            intent(out) :: nElemSets, maxNelemInSet
  character(len=*),     dimension(:), intent(in)  :: fileContent                                    !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer                            :: i,l,elemInCurrentSet

  nElemSets     = 0
  maxNelemInSet = 0

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 2) cycle
    if(IO_lc(IO_StringValue(fileContent(l),chunkPos,1)) == 'define' .and. &
       IO_lc(IO_StringValue(fileContent(l),chunkPos,2)) == 'element') then
      nElemSets = nElemSets + 1

      chunkPos = IO_stringPos(fileContent(l+1))
      if(containsRange(fileContent(l+1),chunkPos)) then
        elemInCurrentSet = 1 + abs( IO_intValue(fileContent(l+1),chunkPos,3) &
                                   -IO_intValue(fileContent(l+1),chunkPos,1))
      else
        elemInCurrentSet = 0
        i = 0
        do while (.true.)
          i = i + 1
          chunkPos = IO_stringPos(fileContent(l+i))
          elemInCurrentSet = elemInCurrentSet + chunkPos(1) - 1                                     ! add line's count when assuming 'c'
          if(IO_lc(IO_stringValue(fileContent(l+i),chunkPos,chunkPos(1))) /= 'c') then              ! line finished, read last value
            elemInCurrentSet = elemInCurrentSet + 1                                                 ! data ended
            exit
          endif
        enddo
      endif
      maxNelemInSet = max(maxNelemInSet, elemInCurrentSet)
    endif
  enddo

end subroutine inputRead_NelemSets


!--------------------------------------------------------------------------------------------------
!> @brief map element sets
!--------------------------------------------------------------------------------------------------
subroutine inputRead_mapElemSets(nameElemSet,mapElemSet,&
                                 fileContent)

  character(len=pStringLen), dimension(:),   allocatable, intent(out) :: nameElemSet
  integer,                   dimension(:,:), allocatable, intent(out) :: mapElemSet
  character(len=*),          dimension(:),                intent(in)  :: fileContent                !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: elemSet, NelemSets, maxNelemInSet,l


  call inputRead_NelemSets(NelemSets,maxNelemInSet,fileContent)
  allocate(nameElemSet(NelemSets)); nameElemSet = 'n/a'
  allocate(mapElemSet(1+maxNelemInSet,NelemSets),source=0)
  elemSet = 0

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 2) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'define' .and. &
       IO_lc(IO_stringValue(fileContent(l),chunkPos,2)) == 'element') then
       elemSet = elemSet+1
       nameElemSet(elemSet)  = trim(IO_stringValue(fileContent(l),chunkPos,4))
       mapElemSet(:,elemSet) = continuousIntValues(fileContent(l+1:),size(mapElemSet,1)-1,nameElemSet,mapElemSet,size(nameElemSet))
    endif
  enddo

end subroutine inputRead_mapElemSets


!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_mapElems(FEM2DAMASK, &
                              nElems,nNodesPerElem,fileContent)

  integer, allocatable, dimension(:), intent(out) :: FEM2DAMASK

  integer,                            intent(in)  :: nElems, &                                      !< number of elements
                                                     nNodesPerElem                                  !< number of nodes per element
  character(len=*),     dimension(:), intent(in)  :: fileContent                                    !< file content, separated per lines

  integer, dimension(2,nElems)       :: map_unsorted
  integer, allocatable, dimension(:) :: chunkPos
  integer :: i,j,l,nNodesAlreadyRead

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'connectivity') then
      j = 0
      do i = 1,nElems
        chunkPos = IO_stringPos(fileContent(l+1+i+j))
        map_unsorted(:,i) = [IO_intValue(fileContent(l+1+i+j),chunkPos,1),i]
        nNodesAlreadyRead = chunkPos(1) - 2
        do while(nNodesAlreadyRead < nNodesPerElem)                                                 ! read on if not all nodes in one line
          j = j + 1
          chunkPos = IO_stringPos(fileContent(l+1+i+j))
          nNodesAlreadyRead = nNodesAlreadyRead + chunkPos(1)
        enddo
      enddo
      exit
    endif
  enddo

  call math_sort(map_unsorted)
  allocate(FEM2DAMASK(minval(map_unsorted(1,:)):maxval(map_unsorted(1,:))),source=-1)
  do i = 1, nElems
    FEM2DAMASK(map_unsorted(1,i)) = map_unsorted(2,i)
  enddo

end subroutine inputRead_mapElems


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_mapNodes(FEM2DAMASK, &
                              nNodes,fileContent)

  integer, allocatable, dimension(:), intent(out) :: FEM2DAMASK

  integer,                            intent(in)  :: nNodes                                         !< number of nodes
  character(len=*),     dimension(:), intent(in)  :: fileContent                                    !< file content, separated per lines

  integer, dimension(2,nNodes)       :: map_unsorted
  integer, allocatable, dimension(:) :: chunkPos
  integer :: i, l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'coordinates') then
      chunkPos = [1,1,10]
      do i = 1,nNodes
        map_unsorted(:,i) = [IO_intValue(fileContent(l+1+i),chunkPos,1),i]
      enddo
      exit
    endif
  enddo

  call math_sort(map_unsorted)
  allocate(FEM2DAMASK(minval(map_unsorted(1,:)):maxval(map_unsorted(1,:))),source=-1)
  do i = 1, nNodes
    FEM2DAMASK(map_unsorted(1,i)) = map_unsorted(2,i)
  enddo

end subroutine inputRead_mapNodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_elemNodes(nodes, &
                               nNode,fileContent)

  real(pReal), allocatable,  dimension(:,:), intent(out) :: nodes
  integer,                                   intent(in)  :: nNode
  character(len=*),            dimension(:), intent(in)  :: fileContent                             !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i,j,m,l

  allocate(nodes(3,nNode))

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'coordinates') then
      chunkPos = [4,1,10,11,30,31,50,51,70]
      do i=1,nNode
        m = discretization_Marc_FEM2DAMASK_node(IO_intValue(fileContent(l+1+i),chunkPos,1))
        do j = 1,3
          nodes(j,m) = mesh_unitlength * IO_floatValue(fileContent(l+1+i),chunkPos,j+1)
        enddo
      enddo
      exit
    endif
  enddo

end subroutine inputRead_elemNodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets element type (and checks if the whole mesh comprises of only one type)
!--------------------------------------------------------------------------------------------------
subroutine inputRead_elemType(elem, &
                              nElem,fileContent)

  type(tElement),                 intent(out) :: elem
  integer,                        intent(in)  :: nElem
  character(len=*), dimension(:), intent(in)  :: fileContent                                        !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i,j,t,l,remainingChunks

  t = -1
  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'connectivity') then
      j = 0
      do i=1,nElem                                                                                  ! read all elements
        chunkPos = IO_stringPos(fileContent(l+1+i+j))
        if (t == -1) then
          t = mapElemtype(IO_stringValue(fileContent(l+1+i+j),chunkPos,2))
          call elem%init(t)
        else
          if (t /= mapElemtype(IO_stringValue(fileContent(l+1+i+j),chunkPos,2))) call IO_error(191,el=t,ip=i)
        endif
        remainingChunks = elem%nNodes - (chunkPos(1) - 2)
        do while(remainingChunks > 0)
          j = j + 1
          chunkPos = IO_stringPos(fileContent(l+1+i+j))
          remainingChunks = remainingChunks - chunkPos(1)
        enddo
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
      case ( '125')                ! 155, 128 (need test)
        mapElemtype = 2            ! Two-dimensional Plane Strain triangle (155: cubic shape function, 125/128: second order isoparametric)
      case ( '11')
        mapElemtype = 3            ! Arbitrary Quadrilateral Plane-strain
      case ( '27')
        mapElemtype = 4            ! Plane Strain, Eight-node Distorted Quadrilateral
      case ( '54')
        mapElemtype = 5            ! Plane Strain, Eight-node Distorted Quadrilateral with reduced integration
      case ( '134')
        mapElemtype = 6            ! Three-dimensional Four-node Tetrahedron
      !case ( '157')               ! need test
      !  mapElemtype = 7           ! Three-dimensional, Low-order, Tetrahedron, Herrmann Formulations
      case ( '127')
        mapElemtype = 8            ! Three-dimensional Ten-node Tetrahedron
      case ( '136')
        mapElemtype = 9            ! Three-dimensional Arbitrarily Distorted Pentahedral
      case ( '117')                ! 123 (need test)
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


end subroutine inputRead_elemType


!--------------------------------------------------------------------------------------------------
!> @brief Stores node IDs
!--------------------------------------------------------------------------------------------------
function inputRead_connectivityElem(nElem,nNodes,fileContent)

  integer, intent(in) :: &
    nElem, &
    nNodes                                                                                          !< number of nodes per element
  character(len=*), dimension(:), intent(in) :: fileContent                                         !< file content, separated per lines

  integer, dimension(nNodes,nElem) :: &
    inputRead_connectivityElem

  integer, allocatable, dimension(:) :: chunkPos

  integer, dimension(1+nElem) :: contInts
  integer :: i,k,j,t,e,l,nNodesAlreadyRead

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 1) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'connectivity') then
      j = 0
      do i = 1,nElem
        chunkPos = IO_stringPos(fileContent(l+1+i+j))
        e = discretization_Marc_FEM2DAMASK_elem(IO_intValue(fileContent(l+1+i+j),chunkPos,1))
        if (e /= 0) then                                                                            ! disregard non CP elems
          do k = 1,chunkPos(1)-2
            inputRead_connectivityElem(k,e) = &
              discretization_Marc_FEM2DAMASK_node(IO_IntValue(fileContent(l+1+i+j),chunkPos,k+2))
          enddo
          nNodesAlreadyRead = chunkPos(1) - 2
          do while(nNodesAlreadyRead < nNodes)                                                      ! read on if not all nodes in one line
            j = j + 1
            chunkPos = IO_stringPos(fileContent(l+1+i+j))
            do k = 1,chunkPos(1)
              inputRead_connectivityElem(nNodesAlreadyRead+k,e) = &
                discretization_Marc_FEM2DAMASK_node(IO_IntValue(fileContent(l+1+i+j),chunkPos,k))
            enddo
            nNodesAlreadyRead = nNodesAlreadyRead + chunkPos(1)
          enddo
        endif
      enddo
      exit
    endif
  enddo

end function inputRead_connectivityElem


!--------------------------------------------------------------------------------------------------
!> @brief Store material ID
!--------------------------------------------------------------------------------------------------
subroutine inputRead_material(materialAt,&
                              nElem,nNodes,nameElemSet,mapElemSet,initialcondTableStyle,fileContent)

  integer, dimension(:), allocatable, intent(out) :: &
    materialAt
  integer, intent(in) :: &
    nElem, &
    nNodes, &                                                                                       !< number of nodes per element
    initialcondTableStyle
  character(len=*), dimension(:), intent(in) :: nameElemSet
  integer, dimension(:,:),        intent(in) :: mapElemSet                                         !< list of elements in elementSet
  character(len=*), dimension(:), intent(in) :: fileContent                                        !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos

  integer, dimension(1+nElem) :: contInts
  integer :: i,j,t,sv,myVal,e,nNodesAlreadyRead,l,k,m


  allocate(materialAt(nElem),source=0)

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if(chunkPos(1) < 2) cycle
    if(IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'initial' .and. &
       IO_lc(IO_stringValue(fileContent(l),chunkPos,2)) == 'state') then
      k = merge(2,1,initialcondTableStyle == 2)
      chunkPos = IO_stringPos(fileContent(l+k))
      sv = IO_IntValue(fileContent(l+k),chunkPos,1)                                                 ! figure state variable index
      if( (sv == 2)) then                                                                           ! state var 2 is used to identify material from material.yaml
        m = 1
        chunkPos = IO_stringPos(fileContent(l+k+m))
        do while (scan(IO_stringValue(fileContent(l+k+m),chunkPos,1),'+-',back=.true.)>1)           ! is noEfloat value?
          myVal = nint(IO_floatValue(fileContent(l+k+m),chunkPos,1))
          if (initialcondTableStyle == 2) m = m + 2
          contInts = continuousIntValues(fileContent(l+k+m+1:),nElem,nameElemSet,mapElemSet,size(nameElemSet)) ! get affected elements
          do i = 1,contInts(1)
            e = discretization_Marc_FEM2DAMASK_elem(contInts(1+i))
            materialAt(e) = myVal
          enddo
          if (initialcondTableStyle == 0) m = m + 1
        enddo
      endif
    endif
  enddo

  if(any(materialAt < 1)) call IO_error(180)

end subroutine inputRead_material


!--------------------------------------------------------------------------------------------------
!> @brief Calculates cell node coordinates from element node coordinates
!--------------------------------------------------------------------------------------------------
pure subroutine buildCells(connectivity,definition, &
                           elem,connectivity_elem)

  type(tCellNodeDefinition), dimension(:),    intent(out) :: definition                             ! definition of cell nodes for increasing number of parents
  integer,                   dimension(:,:,:),intent(out) :: connectivity

  type(tElement),                             intent(in)  :: elem                                   ! element definition
  integer,                   dimension(:,:),  intent(in)  :: connectivity_elem                      ! connectivity of the elements

  integer,dimension(:),     allocatable :: candidates_local
  integer,dimension(:,:),   allocatable :: parentsAndWeights,candidates_global

  integer :: e, n, c, p, s,i,m,j,nParentNodes,nCellNode,Nelem,candidateID

  Nelem = size(connectivity_elem,2)

!---------------------------------------------------------------------------------------------------
! initialize global connectivity to negative local connectivity
  connectivity = -spread(elem%cell,3,Nelem)                                                         ! local cell node ID

!---------------------------------------------------------------------------------------------------
! set connectivity of cell nodes that coincide with FE nodes (defined by 1 parent node)
! and renumber local (negative) to global (positive) node ID
  do e = 1, Nelem
    do c = 1, elem%NcellNodes
      realNode: if (count(elem%cellNodeParentNodeWeights(:,c) /= 0) == 1) then
        where(connectivity(:,:,e) == -c) connectivity(:,:,e) = connectivity_elem(c,e)
      endif realNode
    enddo
  enddo

  nCellNode = maxval(connectivity_elem)

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
    allocate(definition(nParentNodes-1)%parents(i,nParentNodes))
    allocate(definition(nParentNodes-1)%weights(i,nParentNodes))

    i = 1
    n = 1
    do while(n <= size(candidates_local)*Nelem)
      j=0
      parentsAndWeights(:,1) = candidates_global(1:nParentNodes,n+j)
      parentsAndWeights(:,2) = candidates_global(nParentNodes+1:nParentNodes*2,n+j)

      e = candidates_global(nParentNodes*2+1,n+j)
      c = candidates_global(nParentNodes*2+2,n+j)

      do while (n+j<= size(candidates_local)*Nelem)
        if (any(candidates_global(1:2*nParentNodes,n+j)/=candidates_global(1:2*nParentNodes,n))) exit
        where (connectivity(:,:,candidates_global(nParentNodes*2+1,n+j)) == -candidates_global(nParentNodes*2+2,n+j)) ! still locally defined
          connectivity(:,:,candidates_global(nParentNodes*2+1,n+j)) = nCellNode + 1                                   ! gets current new cell node id
        end where

        j = j+1
      enddo
      nCellNode = nCellNode + 1
      definition(nParentNodes-1)%parents(i,:) = parentsAndWeights(:,1)
      definition(nParentNodes-1)%weights(i,:) = parentsAndWeights(:,2)
      i = i + 1
      n = n+j
    enddo

  enddo

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
!> @brief Calculates cell node coordinates from element node coordinates
!--------------------------------------------------------------------------------------------------
pure function buildCellNodes(node_elem)

  real(pReal),               dimension(:,:), intent(in)  :: node_elem                               !< element nodes
  real(pReal),               dimension(:,:), allocatable :: buildCellNodes                          !< cell node coordinates

  integer :: i, j, k, n


  allocate(buildCellNodes(3,maxval(connectivity_cell)))
  n = size(node_elem,2)
  buildCellNodes(:,1:n) = node_elem                                                                 !< initial nodes coincide with element nodes

  do i = 1, size(cellNodeDefinition)
    do j = 1, size(cellNodeDefinition(i)%parents,1)
      n = n+1
      buildCellNodes(:,n) = 0.0_pReal
      do k = 1, size(cellNodeDefinition(i)%parents,2)
        buildCellNodes(:,n) = buildCellNodes(:,n) &
                            + buildCellNodes(:,cellNodeDefinition(i)%parents(j,k)) &
                            * real(cellNodeDefinition(i)%weights(j,k),pReal)
      enddo
      buildCellNodes(:,n) = buildCellNodes(:,n)/real(sum(cellNodeDefinition(i)%weights(j,:)),pReal)
    enddo
  enddo

end function buildCellNodes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP coordinates as center of cell
!--------------------------------------------------------------------------------------------------
pure function buildIPcoordinates(node_cell)

  real(pReal), dimension(:,:), intent(in)  :: node_cell                                             !< cell node coordinates
  real(pReal), dimension(:,:), allocatable :: buildIPcoordinates                                    !< cell-center/IP coordinates

  integer, dimension(:,:), allocatable :: connectivity_cell_reshaped
  integer :: i, n, NcellNodesPerCell,Ncells


  NcellNodesPerCell = size(connectivity_cell,1)
  Ncells = size(connectivity_cell,2)*size(connectivity_cell,3)
  connectivity_cell_reshaped = reshape(connectivity_cell,[NcellNodesPerCell,Ncells])

  allocate(buildIPcoordinates(3,Ncells))

  do i = 1, size(connectivity_cell_reshaped,2)
    buildIPcoordinates(:,i) = 0.0_pReal
    do n = 1, size(connectivity_cell_reshaped,1)
      buildIPcoordinates(:,i) = buildIPcoordinates(:,i) &
                              + node_cell(:,connectivity_cell_reshaped(n,i))
    enddo
    buildIPcoordinates(:,i) = buildIPcoordinates(:,i)/real(size(connectivity_cell_reshaped,1),pReal)
  enddo

end function buildIPcoordinates


!---------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume.
!> @details The IP volume is calculated differently depending on the cell type.
!> 2D cells assume an element depth of 1.0
!---------------------------------------------------------------------------------------------------
pure function IPvolume(elem,node)

  type(tElement),                intent(in) :: elem
  real(pReal), dimension(:,:),   intent(in) :: node

  real(pReal), dimension(elem%nIPs,size(connectivity_cell,3)) :: IPvolume
  real(pReal), dimension(3) :: x0,x1,x2,x3,x4,x5,x6,x7

  integer :: e,i


  do e = 1,size(connectivity_cell,3)
    do i = 1,elem%nIPs

      select case (elem%cellType)
        case (1)                                                                                    ! 2D 3node
          IPvolume(i,e) = math_areaTriangle(node(1:3,connectivity_cell(1,i,e)), &
                                            node(1:3,connectivity_cell(2,i,e)), &
                                            node(1:3,connectivity_cell(3,i,e)))

        case (2)                                                                                    ! 2D 4node
          IPvolume(i,e) = math_areaTriangle(node(1:3,connectivity_cell(1,i,e)), &                   ! assume planar shape, division in two triangles suffices
                                            node(1:3,connectivity_cell(2,i,e)), &
                                            node(1:3,connectivity_cell(3,i,e))) &
                        + math_areaTriangle(node(1:3,connectivity_cell(3,i,e)), &
                                            node(1:3,connectivity_cell(4,i,e)), &
                                            node(1:3,connectivity_cell(1,i,e)))
        case (3)                                                                                    ! 3D 4node
          IPvolume(i,e) = math_volTetrahedron(node(1:3,connectivity_cell(1,i,e)), &
                                              node(1:3,connectivity_cell(2,i,e)), &
                                              node(1:3,connectivity_cell(3,i,e)), &
                                              node(1:3,connectivity_cell(4,i,e)))
        case (4)                                                                                    ! 3D 8node
          ! J. Grandy, Efficient Calculation of Volume of Hexahedral Cells
          ! Lawrence Livermore National Laboratory
          ! https://www.osti.gov/servlets/purl/632793
          x0 = node(1:3,connectivity_cell(1,i,e))
          x1 = node(1:3,connectivity_cell(2,i,e))
          x2 = node(1:3,connectivity_cell(4,i,e))
          x3 = node(1:3,connectivity_cell(3,i,e))
          x4 = node(1:3,connectivity_cell(5,i,e))
          x5 = node(1:3,connectivity_cell(6,i,e))
          x6 = node(1:3,connectivity_cell(8,i,e))
          x7 = node(1:3,connectivity_cell(7,i,e))
          IPvolume(i,e) = dot_product((x7-x1)+(x6-x0),math_cross((x7-x2),        (x3-x0))) &
                        + dot_product((x6-x0),        math_cross((x7-x2)+(x5-x0),(x7-x4))) &
                        + dot_product((x7-x1),        math_cross((x5-x0),        (x7-x4)+(x3-x0)))
          IPvolume(i,e) = IPvolume(i,e)/12.0_pReal
      end select
    enddo
  enddo

end function IPvolume


!--------------------------------------------------------------------------------------------------
!> @brief calculation of IP interface areas
!--------------------------------------------------------------------------------------------------
pure function IPareaNormal(elem,nElem,node)

  type(tElement),                intent(in) :: elem
  integer,                       intent(in) :: nElem
  real(pReal), dimension(:,:),   intent(in) :: node

  real(pReal), dimension(3,elem%nIPneighbors,elem%nIPs,nElem) :: ipAreaNormal

  real(pReal), dimension (3,size(elem%cellFace,1)) :: nodePos
  integer :: e,i,f,n,m

  m = size(elem%cellFace,1)

  do e = 1,nElem
    do i = 1,elem%nIPs
      do f = 1,elem%nIPneighbors
        nodePos = node(1:3,connectivity_cell(elem%cellface(1:m,f),i,e))

        select case (elem%cellType)
          case (1,2)                                                                                ! 2D 3 or 4 node
            IPareaNormal(1,f,i,e) =   nodePos(2,2) - nodePos(2,1)                                   ! x_normal =  y_connectingVector
            IPareaNormal(2,f,i,e) = -(nodePos(1,2) - nodePos(1,1))                                  ! y_normal = -x_connectingVector
            IPareaNormal(3,f,i,e) = 0.0_pReal
          case (3)                                                                                  ! 3D 4node
            IPareaNormal(1:3,f,i,e) = math_cross(nodePos(1:3,2) - nodePos(1:3,1), &
                                                 nodePos(1:3,3) - nodePos(1:3,1))
          case (4)                                                                                  ! 3D 8node
            ! Get the normal of the quadrilateral face as the average of four normals of triangular
            ! subfaces. Since the face consists only of two triangles, the sum has to be divided
            ! by two. This procedure tries to compensate for probable non-planar cell surfaces
            IPareaNormal(1:3,f,i,e) = 0.0_pReal
            do n = 1, m
              IPareaNormal(1:3,f,i,e) = IPareaNormal(1:3,f,i,e) &
                                      + math_cross(nodePos(1:3,mod(n+0,m)+1) - nodePos(1:3,n), &
                                                   nodePos(1:3,mod(n+1,m)+1) - nodePos(1:3,n)) * 0.5_pReal
            enddo
        end select
      enddo
    enddo
  enddo

end function IPareaNormal


!--------------------------------------------------------------------------------------------------
!> @brief IP neighborhood
!--------------------------------------------------------------------------------------------------
function IPneighborhood(elem)

  type(tElement),            intent(in) :: elem                                                     ! definition of the element in use
  integer, dimension(3,size(elem%cellFace,2), &
                     size(connectivity_cell,2),size(connectivity_cell,3)) :: IPneighborhood                   ! neighboring IPs as [element ID, IP ID, face ID]

  integer, dimension(size(elem%cellFace,1)+3,&
                     size(elem%cellFace,2)*size(connectivity_cell,2)*size(connectivity_cell,3)) :: face
  integer, dimension(size(connectivity_cell,1))  :: myConnectivity
  integer, dimension(size(elem%cellFace,1))      :: face_unordered
  integer :: e,i,f,n,c,s

  c = 0
  do e = 1, size(connectivity_cell,3)
    do i = 1, size(connectivity_cell,2)
      myConnectivity = connectivity_cell(:,i,e)
      do f = 1, size(elem%cellFace,2)
        c = c + 1
        face_unordered = myConnectivity(elem%cellFace(:,f))
        do n = 1, size(face_unordered)
          face(n,c) = minval(face_unordered)
          face_unordered(minloc(face_unordered)) = huge(face_unordered)
        enddo
        face(n:n+3,c) = [e,i,f]
      enddo
  enddo; enddo

!--------------------------------------------------------------------------------------------------
! sort face definitions
  call math_sort(face,sortDim=1)
  do c=2, size(face,1)-4
    s = 1
    e = 1
    do while (e < size(face,2))
      e = e + 1
      if(any(face(:c,s) /= face(:c,e))) then
        if(e-1/=s) call math_sort(face(:,s:e-1),sortDim=c)
        s = e
      endif
    enddo
  enddo

  IPneighborhood = 0
  do c=1, size(face,2) - 1
    if(all(face(:n-1,c) == face(:n-1,c+1))) then
      IPneighborhood(:,face(n+2,c+1),face(n+1,c+1),face(n+0,c+1)) = face(n:n+3,c+0)
      IPneighborhood(:,face(n+2,c+0),face(n+1,c+0),face(n+0,c+0)) = face(n:n+3,c+1)
    endif
  enddo

end function IPneighborhood


!--------------------------------------------------------------------------------------------------
!> @brief return integer list corresponding to items in consecutive lines.
!! First integer in array is counter
!> @details ints concatenated by "c" as last char, range of a "to" b, or named set
!--------------------------------------------------------------------------------------------------
function continuousIntValues(fileContent,maxN,lookupName,lookupMap,lookupMaxN)

  character(len=*), dimension(:),   intent(in) :: fileContent                                       !< file content, separated per lines
  integer,                          intent(in) :: maxN
  integer,                          intent(in) :: lookupMaxN
  integer,          dimension(:,:), intent(in) :: lookupMap
  character(len=*), dimension(:),   intent(in) :: lookupName

  integer,           dimension(1+maxN)         :: continuousIntValues

  integer :: l,i,first,last
  integer, allocatable, dimension(:) :: chunkPos
  logical :: rangeGeneration

  continuousIntValues = 0
  rangeGeneration = .false.

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if (chunkPos(1) < 1) then                                                                       ! empty line
      exit
    elseif (verify(IO_stringValue(fileContent(l),chunkPos,1),'0123456789') > 0) then                ! a non-int, i.e. set name
      do i = 1, lookupMaxN                                                                          ! loop over known set names
        if (IO_stringValue(fileContent(l),chunkPos,1) == lookupName(i)) then                        ! found matching name
          continuousIntValues = lookupMap(:,i)                                                      ! return resp. entity list
          exit
        endif
      enddo
      exit
    elseif(containsRange(fileContent(l),chunkPos)) then
      first = IO_intValue(fileContent(l),chunkPos,1)
      last  = IO_intValue(fileContent(l),chunkPos,3)
      do i = first, last, sign(1,last-first)
        continuousIntValues(1) = continuousIntValues(1) + 1
        continuousIntValues(1+continuousIntValues(1)) = i
      enddo
      exit
    else
      do i = 1,chunkPos(1)-1                                                                        ! interpret up to second to last value
        continuousIntValues(1) = continuousIntValues(1) + 1
        continuousIntValues(1+continuousIntValues(1)) = IO_intValue(fileContent(l),chunkPos,i)
      enddo
      if ( IO_lc(IO_stringValue(fileContent(l),chunkPos,chunkPos(1))) /= 'c' ) then                 ! line finished, read last value
        continuousIntValues(1) = continuousIntValues(1) + 1
        continuousIntValues(1+continuousIntValues(1)) = IO_intValue(fileContent(l),chunkPos,chunkPos(1))
        exit
      endif
    endif
  enddo

end function continuousIntValues


!--------------------------------------------------------------------------------------------------
!> @brief return whether a line contains a range ('X to Y')
!--------------------------------------------------------------------------------------------------
logical function containsRange(str,chunkPos)

  character(len=*),      intent(in) :: str
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string

  containsRange = .False.
  if(chunkPos(1) == 3) then
     if(IO_lc(IO_stringValue(str,chunkPos,2)) == 'to') containsRange = .True.
  endif

end function containsRange


end module discretization_marc
