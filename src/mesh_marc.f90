!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
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
  
  type tCellNodeDefinition
    integer, dimension(:,:), allocatable :: parents
    integer, dimension(:,:), allocatable :: weights
  end type tCellNodeDefinition
  
  type(tCellNodeDefinition), dimension(:), allocatable :: cellNodeDefinition
  
  real(pReal), public, protected :: &
    mesh_unitlength                                                                                 !< physical length of one unit in mesh
 
  integer, dimension(:,:), allocatable, target :: &
    mesh_mapFEtoCPelem, &                                                                           !< [sorted FEid, corresponding CPid]
    mesh_mapFEtoCPnode                                                                              !< [sorted FEid, corresponding CPid]

!-------------------------------------------------------------------------------------------------- 
! DEPRECATED
  real(pReal), dimension(:,:,:), allocatable, public :: &
    mesh_ipCoordinates                                                                              !< IP x,y,z coordinates (after deformation!)
!--------------------------------------------------------------------------------------------------

 public :: &
   mesh_init, &
   mesh_FEasCP
 
 
contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the mesh by calling all necessary private routines the mesh module
!! Order and routines strongly depend on type of solver
!--------------------------------------------------------------------------------------------------
subroutine mesh_init(ip,el)

  integer, intent(in) :: el, ip
   
  real(pReal), dimension(:,:), allocatable :: &
   node0_elem, &                                                                                       !< node x,y,z coordinates (initially!)
   node0_cell
  type(tElement) :: elem

  integer :: nElems
  integer,     dimension(:),   allocatable :: &
    microstructureAt, &
    homogenizationAt
  integer:: &
    Nnodes                                                                                      !< total number of nodes in mesh
   
  real(pReal), dimension(:,:), allocatable :: &
    ip_reshaped
  integer,dimension(:,:,:), allocatable :: &
    connectivity_cell                                                                                !< cell connectivity for each element,ip/cell
  integer, dimension(:,:), allocatable :: &
    connectivity_elem
  real(pReal), dimension(:,:,:,:),allocatable :: &
    unscaledNormals

  write(6,'(/,a)')   ' <<<+-  mesh init  -+>>>'
 
  mesh_unitlength = numerics_unitlength                                                             ! set physical extent of a length unit in mesh

  call inputRead(elem,node0_elem,connectivity_elem,microstructureAt,homogenizationAt)
  nElems = size(connectivity_elem,2)

  if (debug_e < 1 .or. debug_e > nElems)    call IO_error(602,ext_msg='element')
  if (debug_i < 1 .or. debug_i > elem%nIPs) call IO_error(602,ext_msg='IP')
 
  FEsolving_execElem = [ 1,nElems ]                                                                 ! parallel loop bounds set to comprise all DAMASK elements
  allocate(FEsolving_execIP(2,nElems), source=1)                                                    ! parallel loop bounds set to comprise from first IP...
  FEsolving_execIP(2,:) = elem%nIPs
 
  allocate(calcMode(elem%nIPs,nElems),source=.false.)                                               ! pretend to have collected what first call is asking (F = I)
  calcMode(ip,mesh_FEasCP('elem',el)) = .true.                                                      ! first ip,el needs to be already pingponged to "calc"
 

  allocate(mesh_ipCoordinates(3,elem%nIPs,nElems),source=0.0_pReal)                                 ! deprecated

  allocate(cellNodeDefinition(elem%nNodes-1))
  allocate(connectivity_cell(elem%NcellNodesPerCell,elem%nIPs,nElems))
  call buildCells(connectivity_cell,cellNodeDefinition,&
                  elem,connectivity_elem)
  allocate(node0_cell(3,maxval(connectivity_cell)))
  call buildCellNodes(node0_cell,&
                      cellNodeDefinition,node0_elem)
  allocate(ip_reshaped(3,elem%nIPs*nElems),source=0.0_pReal)
  call buildIPcoordinates(ip_reshaped,reshape(connectivity_cell,[elem%NcellNodesPerCell,&
                          elem%nIPs*nElems]),node0_cell)

  call discretization_init(microstructureAt,homogenizationAt,&
                           ip_reshaped,&
                           node0_elem)
                           
  call writeGeometry(0,connectivity_elem,&
                     reshape(connectivity_cell,[elem%NcellNodesPerCell,elem%nIPs*nElems]),&
                     node0_cell,ip_reshaped)

!--------------------------------------------------------------------------------------------------
! geometry information required by the nonlocal CP model
  call geometry_plastic_nonlocal_setIPvolume(IPvolume(elem,node0_cell,connectivity_cell))
  unscaledNormals = IPareaNormal(elem,nElems,connectivity_cell,node0_cell)
  call geometry_plastic_nonlocal_setIParea(norm2(unscaledNormals,1))
  call geometry_plastic_nonlocal_setIPareaNormal(unscaledNormals/spread(norm2(unscaledNormals,1),1,3))
  call geometry_plastic_nonlocal_results
  
end subroutine mesh_init


!--------------------------------------------------------------------------------------------------
!> @brief Writes all information needed for the DADF5 geometry
!--------------------------------------------------------------------------------------------------
subroutine writeGeometry(elemType, &
                         connectivity_elem,connectivity_cell, &
                         coordinates_nodes,coordinates_points)

  integer,                     intent(in) :: elemType
  integer, dimension(:,:),     intent(in) :: &
    connectivity_elem, &
    connectivity_cell
  real(pReal), dimension(:,:), intent(in) :: &
    coordinates_nodes, &
    coordinates_points
  
  integer,     dimension(:,:), allocatable :: &
    connectivity_temp
  real(pReal), dimension(:,:), allocatable :: &
    coordinates_temp

#if defined(DAMASK_HDF5)
  call results_openJobFile
  call HDF5_closeGroup(results_addGroup('geometry'))
  
  connectivity_temp = connectivity_elem
  call results_writeDataset('geometry',connectivity_temp,'T_e',&
                            'connectivity of the elements','-')
                            
  connectivity_temp = connectivity_cell
  call results_writeDataset('geometry',connectivity_temp,'T_c', &
                            'connectivity of the cells','-')
                            
  coordinates_temp = coordinates_nodes
  call results_writeDataset('geometry',coordinates_temp,'x_n', &
                            'coordinates of the nodes','m')
                            
  coordinates_temp = coordinates_points
  call results_writeDataset('geometry',coordinates_temp,'x_p', &
                            'coordinates of the material points','m')
                      
  call results_closeJobFile
#endif  

end subroutine writeGeometry


subroutine inputRead(elem,node0_elem,connectivity_elem,microstructureAt,homogenizationAt)

  type(tElement), intent(out) :: elem
  real(pReal), dimension(:,:), allocatable, intent(out) :: &
   node0_elem                                                                                       !< node x,y,z coordinates (initially!)
  integer, dimension(:,:), allocatable, intent(out) :: &
    connectivity_elem

  integer,     dimension(:),   allocatable, intent(out) :: &
    microstructureAt, &
    homogenizationAt
   
  integer :: &
    fileFormatVersion, &
    hypoelasticTableStyle, &
    initialcondTableStyle, &
    nNodes, &
    nElems
  integer, parameter :: &
    FILEUNIT = 222
  integer, dimension(:), allocatable :: &
    matNumber                                                                                       !< material numbers for hypoelastic material
  character(len=pStringLen), dimension(:), allocatable :: inputFile                                 !< file content, separated per lines
  
  character(len=64), dimension(:), allocatable :: &
    nameElemSet
  integer, dimension(:,:), allocatable :: &
    mapElemSet                                                                                      !< list of elements in elementSet

  inputFile = IO_read_ASCII(trim(getSolverJobName())//trim(InputFileExtension))
  call inputRead_fileFormat(fileFormatVersion, &
                            inputFile)
  call inputRead_tableStyles(initialcondTableStyle,hypoelasticTableStyle, &
                             inputFile)
  if (fileFormatVersion > 12) &
    call inputRead_matNumber(matNumber, &
                             hypoelasticTableStyle,inputFile)
  call inputRead_NnodesAndElements(nNodes,nElems,&
                                   inputFile)
  
  call IO_open_inputFile(FILEUNIT)                                                                  ! ToDo: It would be better to use fileContent
  
  call inputRead_mapElemSets(nameElemSet,mapElemSet,&
                             FILEUNIT)

  allocate (mesh_mapFEtoCPelem(2,nElems), source=0)
  call inputRead_mapElems(hypoelasticTableStyle,nameElemSet,mapElemSet,fileFormatVersion,matNumber,FILEUNIT)

  allocate (mesh_mapFEtoCPnode(2,Nnodes), source=0)
  call inputRead_mapNodes(inputFile)

  call inputRead_elemType(elem, &
                          nElems,inputFile)
  call inputRead_elemNodes(node0_elem, &
                           Nnodes,inputFile)

  connectivity_elem = inputRead_connectivityElem(nElems,elem%nNodes,inputFile)

  call inputRead_microstructureAndHomogenization(microstructureAt,homogenizationAt, &
                                                 nElems,elem%nNodes,nameElemSet,mapElemSet,&
                                                 initialcondTableStyle,FILEUNIT)
  close(FILEUNIT)
   
end subroutine inputRead



!--------------------------------------------------------------------------------------------------
!> @brief Figures out version of Marc input file format
!--------------------------------------------------------------------------------------------------
subroutine inputRead_fileFormat(fileFormat,fileContent)
 
  integer,                                 intent(out) :: fileFormat
  character(len=pStringLen), dimension(:), intent(in)  :: fileContent                                !< file content, separated per lines
  
  integer, allocatable, dimension(:) :: chunkPos
  integer :: l
  
  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if ( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'version') then
      fileFormat = IO_intValue(fileContent(l),chunkPos,2)
      exit
    endif
  enddo

end subroutine inputRead_fileFormat


!--------------------------------------------------------------------------------------------------
!> @brief Figures out table styles for initial cond and hypoelastic
!--------------------------------------------------------------------------------------------------
subroutine inputRead_tableStyles(initialcond,hypoelastic,fileContent)
 
  integer,                                 intent(out) :: initialcond, hypoelastic
  character(len=pStringLen), dimension(:), intent(in)  :: fileContent                               !< file content, separated per lines

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

end subroutine inputRead_tableStyles


!--------------------------------------------------------------------------------------------------
!> @brief Figures out material number of hypoelastic material
!--------------------------------------------------------------------------------------------------
subroutine inputRead_matNumber(matNumber, &
                               tableStyle,fileContent)
 
  integer, allocatable,      dimension(:), intent(out) :: matNumber
  integer,                                 intent(in)  :: tableStyle
  character(len=pStringLen), dimension(:), intent(in)  :: fileContent                               !< file content, separated per lines


  integer, allocatable, dimension(:) :: chunkPos
  integer :: i, j, data_blocks, l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if ( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'hypoelastic') then
      if (len(trim(fileContent(l+1)))/=0) then
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
  
  integer,                                 intent(out) :: nNodes, nElems
  character(len=pStringLen), dimension(:), intent(in)  :: fileContent                               !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: l

  nNodes = 0
  nElems = 0

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if     (IO_lc(IO_StringValue(fileContent(l),chunkPos,1)) == 'sizing') then
      nElems = IO_IntValue (fileContent(l),chunkPos,3)
    elseif (IO_lc(IO_StringValue(fileContent(l),chunkPos,1)) == 'coordinates') then
      chunkPos = IO_stringPos(fileContent(l+1))
      nNodes = IO_IntValue (fileContent(l+1),chunkPos,2)
    endif
  enddo

end subroutine inputRead_NnodesAndElements


!--------------------------------------------------------------------------------------------------
!> @brief Count overall number of element sets in mesh.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_NelemSets(nElemSets,maxNelemInSet,&
                               fileUnit)
 
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

620 end subroutine inputRead_NelemSets


!--------------------------------------------------------------------------------------------------
!> @brief map element sets
!--------------------------------------------------------------------------------------------------
subroutine inputRead_mapElemSets(nameElemSet,mapElemSet,fileUnit)
 
  character(len=64), dimension(:),   allocatable, intent(out) :: nameElemSet
  integer,           dimension(:,:), allocatable, intent(out) :: mapElemSet
  integer,                                        intent(in)  :: fileUnit

  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) :: line
  integer :: elemSet, NelemSets, maxNelemInSet


  call inputRead_NelemSets(NelemSets,maxNelemInSet,fileUnit) 
  allocate(nameElemSet(NelemSets)); nameElemSet = 'n/a'
  allocate(mapElemSet(1+maxNelemInSet,NelemSets),source=0)
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

620 end subroutine inputRead_mapElemSets



!--------------------------------------------------------------------------------------------------
!> @brief Maps elements from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_mapElems(tableStyle,nameElemSet,mapElemSet,fileFormatVersion,matNumber,fileUnit)
 
  integer, intent(in) :: fileUnit,tableStyle,fileFormatVersion
  integer, dimension(:), intent(in) :: matNumber
  character(len=64), intent(in), dimension(:) :: nameElemSet
  integer, dimension(:,:), intent(in) :: &
    mapElemSet
 
  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) :: line, &
                        tmp
 
  integer, dimension(:), allocatable :: contInts
  integer :: i,cpElem 
 
  allocate(contInts(size(mesh_mapFEtoCPelem,2)+1))
 
  cpElem = 0
  contInts = 0
  rewind(fileUnit)
  do
    read (fileUnit,'(A300)',END=620) line
    chunkPos = IO_stringPos(line)
    Marc2016andEarlier: if (fileFormatVersion < 13) then
      if( IO_lc(IO_stringValue(line,chunkPos,1)) == 'hypoelastic' ) then
        skipLines: do i=1,3+TableStyle
          read (fileUnit,'(A300)') line
        enddo skipLines
        contInts = IO_continuousIntValues(fileUnit,size(mesh_mapFEtoCPelem,2),nameElemSet,&
                                          mapElemSet,size(nameElemSet))
        exit
      endif  
    else Marc2016andEarlier
      if ( IO_lc(IO_stringValue(line,chunkPos,1)) == 'connectivity') then
        read (fileUnit,'(A300)',END=620) line
        chunkPos = IO_stringPos(line)
        if(any(matNumber==IO_intValue(line,chunkPos,6))) then
          do 
            read (fileUnit,'(A300)',END=620) line
            chunkPos = IO_stringPos(line)
            tmp = IO_lc(IO_stringValue(line,chunkPos,1))
            if (verify(trim(tmp),"0123456789")/=0) then                                             ! found keyword
              exit
            else
              contInts(1) = contInts(1) + 1  
              read (tmp,*) contInts(contInts(1)+1)     
            endif
          enddo
        endif  
      endif
    endif Marc2016andEarlier
  enddo    
620 do i = 1,contInts(1)
    cpElem = cpElem+1
    mesh_mapFEtoCPelem(1,cpElem) = contInts(1+i)
    mesh_mapFEtoCPelem(2,cpElem) = cpElem
  enddo
 
call math_sort(mesh_mapFEtoCPelem)

end subroutine inputRead_mapElems


!--------------------------------------------------------------------------------------------------
!> @brief Maps node from FE ID to internal (consecutive) representation.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_mapNodes(fileContent)

  character(len=pStringLen), dimension(:), intent(in) :: fileContent                                !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i, l

  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'coordinates' ) then
      do i = 1,size(mesh_mapFEtoCPnode,2)
        mesh_mapFEtoCPnode(1,i) = IO_fixedIntValue (fileContent(l+1+i),[0,10],1)
        mesh_mapFEtoCPnode(2,i) = i
      enddo
      exit
    endif
  enddo

  call math_sort(mesh_mapFEtoCPnode)

end subroutine inputRead_mapNodes


!--------------------------------------------------------------------------------------------------
!> @brief store x,y,z coordinates of all nodes in mesh.
!--------------------------------------------------------------------------------------------------
subroutine inputRead_elemNodes(nodes, &
                               nNode,fileContent)

  real(pReal), allocatable,  dimension(:,:), intent(out) :: nodes 
  integer,                                   intent(in)  :: nNode
  character(len=pStringLen), dimension(:),   intent(in)  :: fileContent                             !< file content, separated per lines
  
  integer, dimension(5), parameter   :: node_ends = [0,10,30,50,70]
  integer, allocatable, dimension(:) :: chunkPos
  integer :: i,j,m,l

  allocate(nodes(3,nNode))

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

end subroutine inputRead_elemNodes


!--------------------------------------------------------------------------------------------------
!> @brief Gets element type (and checks if the whole mesh comprises of only one type)
!--------------------------------------------------------------------------------------------------
subroutine inputRead_elemType(elem, &
                              nElem,fileContent)

  type(tElement),                          intent(out) :: elem
  integer,                                 intent(in)  :: nElem
  character(len=pStringLen), dimension(:), intent(in)  :: fileContent                               !< file content, separated per lines

  integer, allocatable, dimension(:) :: chunkPos
  integer :: i,j,t,l,remainingChunks

  t = -1
  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'connectivity' ) then
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


end subroutine inputRead_elemType


!--------------------------------------------------------------------------------------------------
!> @brief Stores node IDs
!--------------------------------------------------------------------------------------------------
function inputRead_connectivityElem(nElem,nNodes,fileContent)
 
  integer, intent(in) :: &
    nElem, &
    nNodes                                                                                          !< number of nodes per element
  character(len=pStringLen), dimension(:), intent(in)  :: fileContent                               !< file content, separated per lines
    
  integer, dimension(nNodes,nElem) :: &
    inputRead_connectivityElem
 
  integer, allocatable, dimension(:) :: chunkPos
 
  integer, dimension(1+nElem) :: contInts
  integer :: i,k,j,t,e,l,nNodesAlreadyRead
 
  do l = 1, size(fileContent)
    chunkPos = IO_stringPos(fileContent(l))
    if( IO_lc(IO_stringValue(fileContent(l),chunkPos,1)) == 'connectivity' ) then
      j = 0
      do i = 1,nElem
        chunkPos = IO_stringPos(fileContent(l+1+i+j))
        e = mesh_FEasCP('elem',IO_intValue(fileContent(l+1+i+j),chunkPos,1))
        if (e /= 0) then                                                                            ! disregard non CP elems
          do k = 1,chunkPos(1)-2
            inputRead_connectivityElem(k,e) = &
              mesh_FEasCP('node',IO_IntValue(fileContent(l+1+i+j),chunkPos,k+2))
          enddo
          nNodesAlreadyRead = chunkPos(1) - 2
          do while(nNodesAlreadyRead < nNodes)                                                      ! read on if not all nodes in one line
            j = j + 1
            chunkPos = IO_stringPos(fileContent(l+1+i+j))
            do k = 1,chunkPos(1)
              inputRead_connectivityElem(nNodesAlreadyRead+k,e) = &
                mesh_FEasCP('node',IO_IntValue(fileContent(l+1+i+j),chunkPos,k))
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
!> @brief Stores homogenization and microstructure ID
!--------------------------------------------------------------------------------------------------
subroutine inputRead_microstructureAndHomogenization(microstructureAt,homogenizationAt, &
                                    nElem,nNodes,nameElemSet,mapElemSet,initialcondTableStyle,fileUnit)

  integer, dimension(:), allocatable, intent(out) :: &
    microstructureAt, &
    homogenizationAt
  integer, intent(in) :: &
    nElem, &
    nNodes, &                                                                                       !< number of nodes per element
    initialcondTableStyle, &
    fileUnit
  character(len=64), dimension(:), intent(in) :: &
    nameElemSet
  integer, dimension(:,:), intent(in) :: &
    mapElemSet                                                                                      !< list of elements in elementSet

  integer, allocatable, dimension(:) :: chunkPos
  character(len=300) line
 
  integer, dimension(1+nElem) :: contInts
  integer :: i,j,t,sv,myVal,e,nNodesAlreadyRead


  allocate(microstructureAt(nElem),source=0)
  allocate(homogenizationAt(nElem),source=0)

  rewind(fileUnit)
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
                    (fileUnit,nElem,nameElemSet,mapElemSet,size(nameElemSet))
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
 
630 end subroutine inputRead_microstructureAndHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief Calculates cell node coordinates from element node coordinates
!--------------------------------------------------------------------------------------------------
subroutine buildCells(connectivity_cell,cellNodeDefinition, &
                      elem,connectivity_elem)

  type(tCellNodeDefinition), dimension(:),    intent(out) :: cellNodeDefinition                     ! definition of cell nodes for increasing number of parents
  integer,                   dimension(:,:,:),intent(out) :: connectivity_cell
  
  type(tElement),                             intent(in)  :: elem                                   ! element definition
  integer,                   dimension(:,:),  intent(in)  :: connectivity_elem                      ! connectivity of the elements
  
  integer,dimension(:),     allocatable :: candidates_local
  integer,dimension(:,:),   allocatable :: parentsAndWeights,candidates_global
  
  integer :: e, n, c, p, s,i,m,j,nParentNodes,nCellNode,Nelem,candidateID
  
  Nelem = size(connectivity_elem,2)

!---------------------------------------------------------------------------------------------------
! initialize global connectivity to negative local connectivity
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
    allocate(cellNodeDefinition(nParentNodes-1)%parents(i,nParentNodes))
    allocate(cellNodeDefinition(nParentNodes-1)%weights(i,nParentNodes))
    
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
        where (connectivity_cell(:,:,candidates_global(nParentNodes*2+1,n+j)) == -candidates_global(nParentNodes*2+2,n+j)) ! still locally defined
          connectivity_cell(:,:,candidates_global(nParentNodes*2+1,n+j)) = nCellNode + 1                                   ! gets current new cell node id
        end where
       
        j = j+1
      enddo
      nCellNode = nCellNode + 1
      cellNodeDefinition(nParentNodes-1)%parents(i,:) = parentsAndWeights(:,1)
      cellNodeDefinition(nParentNodes-1)%weights(i,:) = parentsAndWeights(:,2)
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
subroutine buildCellNodes(node_cell, &
                          definition,node_elem)

  real(pReal),               dimension(:,:), intent(out) :: node_cell                               !< cell node coordinates
  type(tCellNodeDefinition), dimension(:),   intent(in)  :: definition                              !< cell node definition (weights and parents)
  real(pReal),               dimension(:,:), intent(in)  :: node_elem                               !< element nodes
  
  integer :: i, j, k, n
  
  n = size(node_elem,2)
  node_cell(:,1:n) = node_elem                                                                      !< initial nodes coincide with element nodes

  do i = 1, size(cellNodeDefinition,1)
    do j = 1, size(cellNodeDefinition(i)%parents,1)
      n = n+1
      node_cell(:,n) = 0.0_pReal
      do k = 1, size(cellNodeDefinition(i)%parents,2)
        node_cell(:,n) = node_cell(:,n) &
                       + node_cell(:,definition(i)%parents(j,k)) * real(definition(i)%weights(j,k),pReal)
      enddo
      node_cell(:,n) = node_cell(:,n)/real(sum(definition(i)%weights(j,:)),pReal)
    enddo
  enddo
  
end subroutine buildCellNodes


!--------------------------------------------------------------------------------------------------
!> @brief Calculates IP coordinates as center of cell
!--------------------------------------------------------------------------------------------------
subroutine buildIPcoordinates(IPcoordinates, &
                              connectivity_cell,node_cell)

  real(pReal), dimension(:,:), intent(out):: IPcoordinates                                          !< cell-center/IP coordinates
  integer,     dimension(:,:), intent(in) :: connectivity_cell                                      !< connectivity for each cell
  real(pReal), dimension(:,:), intent(in) :: node_cell                                              !< cell node coordinates
  
  integer :: i, n

  do i = 1, size(connectivity_cell,2)
    IPcoordinates(:,i) = 0.0_pReal
    do n = 1, size(connectivity_cell,1)
      IPcoordinates(:,i) = IPcoordinates(:,i) &
                         + node_cell(:,connectivity_cell(n,i))
    enddo
    IPcoordinates(:,i) = IPcoordinates(:,i)/real(size(connectivity_cell,1),pReal)
  enddo
  
end subroutine buildIPcoordinates


!---------------------------------------------------------------------------------------------------
!> @brief Calculates IP volume.
!> @details The IP volume is calculated differently depending on the cell type.
!> 2D cells assume an element depth of 1.0
!---------------------------------------------------------------------------------------------------
function IPvolume(elem,node,connectivity)
 
  type(tElement),                intent(in) :: elem
  real(pReal), dimension(:,:),   intent(in) :: node
  integer,     dimension(:,:,:), intent(in) :: connectivity
  
  real(pReal), dimension(elem%nIPs,size(connectivity,3)) :: IPvolume
  real(pReal), dimension(3) :: x0,x1,x2,x3,x4,x5,x6,x7

  integer :: e,i

  do e = 1,size(connectivity,3)
    do i = 1,elem%nIPs

      select case (elem%cellType)
        case (1)                                                                                    ! 2D 3node
          IPvolume(i,e) = math_areaTriangle(node(1:3,connectivity(1,i,e)), &
                                            node(1:3,connectivity(2,i,e)), &
                                            node(1:3,connectivity(3,i,e)))

        case (2)                                                                                    ! 2D 4node
          IPvolume(i,e) = math_areaTriangle(node(1:3,connectivity(1,i,e)), &                        ! assume planar shape, division in two triangles suffices
                                            node(1:3,connectivity(2,i,e)), &
                                            node(1:3,connectivity(3,i,e))) &
                        + math_areaTriangle(node(1:3,connectivity(3,i,e)), &
                                            node(1:3,connectivity(4,i,e)), &
                                            node(1:3,connectivity(1,i,e)))
        case (3)                                                                                    ! 3D 4node
          IPvolume(i,e) = math_volTetrahedron(node(1:3,connectivity(1,i,e)), &
                                              node(1:3,connectivity(2,i,e)), &
                                              node(1:3,connectivity(3,i,e)), &
                                              node(1:3,connectivity(4,i,e)))
        case (4)                                                                                    ! 3D 8node
          ! J. Grandy, Efficient Calculation of Volume of Hexahedral Cells
          ! Lawrence Livermore National Laboratory
          ! https://www.osti.gov/servlets/purl/632793
          x0 = node(1:3,connectivity(1,i,e))
          x1 = node(1:3,connectivity(2,i,e))
          x2 = node(1:3,connectivity(4,i,e))
          x3 = node(1:3,connectivity(3,i,e))
          x4 = node(1:3,connectivity(5,i,e))
          x5 = node(1:3,connectivity(6,i,e))
          x6 = node(1:3,connectivity(8,i,e))
          x7 = node(1:3,connectivity(7,i,e))
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
function IPareaNormal(elem,nElem,connectivity,node)

  type(tElement),                intent(in) :: elem
  integer,                       intent(in) :: nElem
  integer,     dimension(:,:,:), intent(in) :: connectivity
  real(pReal), dimension(:,:),   intent(in) :: node
  
  real(pReal), dimension(3,elem%nIPneighbors,elem%nIPs,nElem) :: ipAreaNormal

  real(pReal), dimension (3,size(elem%cellFace,1)) :: nodePos
  integer :: e,i,f,n,m
 
  m = size(elem%cellFace,1)

  do e = 1,nElem
    do i = 1,elem%nIPs
      do f = 1,elem%nIPneighbors
        nodePos = node(1:3,connectivity(elem%cellface(1:m,f),i,e))
        
        select case (elem%cellType)
          case (1,2)                                                                                ! 2D 3 or 4 node
            IPareaNormal(1,f,i,e) =   nodePos(2,2) - nodePos(2,1)                                   ! x_normal =  y_connectingVector
            IPareaNormal(2,f,i,e) = -(nodePos(1,2) - nodePos(1,1))                                  ! y_normal = -x_connectingVector
            IPareaNormal(3,f,i,e) = 0.0_pReal
          case (3)                                                                                  ! 3D 4node
            IPareaNormal(1:3,f,i,e) = math_cross(nodePos(1:3,2) - nodePos(1:3,1), &
                                                 nodePos(1:3,3) - nodePos(1:3,1))
          case (4)                                                                                  ! 3D 8node
            ! for this cell type we get the normal of the quadrilateral face as an average of
            ! four normals of triangular subfaces; since the face consists only of two triangles,
            ! the sum has to be divided by two; this whole prcedure tries to compensate for
            ! probable non-planar cell surfaces
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
