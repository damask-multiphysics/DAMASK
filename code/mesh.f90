!* $Id$
!##############################################################
 MODULE mesh     
!##############################################################

 use prec, only: pReal,pInt
 implicit none

! ---------------------------
! _Nelems    : total number of elements in mesh
! _NcpElems  : total number of CP elements in mesh
! _Nnodes    : total number of nodes in mesh
! _maxNnodes : max number of nodes in any CP element
! _maxNips   : max number of IPs in any CP element
! _maxNipNeighbors   : max number of IP neighbors in any CP element
! _maxNsharedElems : max number of CP elements sharing a node
!
! _element    : FEid, type(internal representation), material, texture, node indices
! _node       : x,y,z coordinates (initially!)
! _sharedElem : entryCount and list of elements containing node
!
! _mapFEtoCPelem : [sorted FEid, corresponding CPid]
! _mapFEtoCPnode : [sorted FEid, corresponding CPid]
!
! MISSING: these definitions should actually reside in the
! FE-solver specific part (different for MARC/ABAQUS)..!
! Hence, I suggest to prefix with "FE_"
!
! _Nnodes            : # nodes in a specific type of element (how we use it)
! _NoriginalNodes    : # nodes in a specific type of element (how it is originally defined by marc)
! _Nips              : # IPs in a specific type of element
! _NipNeighbors      : # IP neighbors in a specific type of element
! _ipNeighbor        : +x,-x,+y,-y,+z,-z list of intra-element IPs and
!     (negative) neighbor faces per own IP in a specific type of element
! _NfaceNodes        : # nodes per face in a specific type of element

! _nodeOnFace        : list of node indices on each face of a specific type of element
! _nodesAtIP         : map IP index to two node indices in a specific type of element
! _ipNeighborhood    : 6 or less neighboring IPs as [element_num, IP_index]
! _NsubNodes        : # subnodes required to fully define all IP volumes

!     order is +x,-x,+y,-y,+z,-z but meaning strongly depends on Elemtype
! ---------------------------
 integer(pInt) mesh_Nelems,mesh_NcpElems,mesh_NelemSets,mesh_maxNelemInSet
 integer(pInt) mesh_Nmaterials
 integer(pInt) mesh_Nnodes,mesh_maxNnodes,mesh_maxNips,mesh_maxNipNeighbors,mesh_maxNsharedElems,mesh_maxNsubNodes
 integer(pInt), dimension(2) :: mesh_maxValStateVar = 0_pInt
 character(len=64), dimension(:),   allocatable :: mesh_nameElemSet, &    ! names of elementSet
                                                   mesh_nameMaterial, &   ! names of material in solid section
                                                   mesh_mapMaterial       ! name of elementSet for material
 integer(pInt), dimension(:,:),     allocatable :: mesh_mapElemSet        ! list of elements in elementSet
 integer(pInt), dimension(:,:),     allocatable, target :: mesh_mapFEtoCPelem, mesh_mapFEtoCPnode
 integer(pInt), dimension(:,:),     allocatable :: mesh_element, mesh_sharedElem
 integer(pInt), dimension(:,:,:,:), allocatable :: mesh_ipNeighborhood

 real(pReal),   dimension(:,:,:),   allocatable :: mesh_subNodeCoord      ! coordinates of subnodes per element
 real(pReal),   dimension(:,:),     allocatable :: mesh_ipVolume          ! volume associated with IP
 real(pReal),   dimension(:,:,:),   allocatable :: mesh_ipArea, &         ! area of interface to neighboring IP
                                                   mesh_ipCenterOfGravity ! center of gravity of IP
 real(pReal),   dimension(:,:,:,:), allocatable :: mesh_ipAreaNormal      ! area normal of interface to neighboring IP
 real(pReal),                       allocatable :: mesh_node (:,:)
 
 integer(pInt), dimension(:,:,:,:), allocatable :: FE_nodesAtIP
 integer(pInt), dimension(:,:,:),   allocatable :: FE_ipNeighbor
 integer(pInt), dimension(:,:,:),   allocatable :: FE_subNodeParent
 integer(pInt), dimension(:,:,:,:), allocatable :: FE_subNodeOnIPFace

 integer(pInt) :: hypoelasticTableStyle
 integer(pInt) :: initialcondTableStyle
 integer(pInt), parameter :: FE_Nelemtypes = 7
 integer(pInt), parameter :: FE_maxNnodes = 8
 integer(pInt), parameter :: FE_maxNsubNodes = 56
 integer(pInt), parameter :: FE_maxNips = 27
 integer(pInt), parameter :: FE_maxNipNeighbors = 6
 integer(pInt), parameter :: FE_NipFaceNodes = 4
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_Nnodes = &
 (/8, & ! element 7
   4, & ! element 134
   4, & ! element 11
   4, & ! element 27
   4, & ! element 157
   6, & ! element 136
   8  & ! element 21
  /)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_NoriginalNodes = &
 (/8, & ! element 7
   4, & ! element 134
   4, & ! element 11
   4, & ! element 27
   4, & ! element 157
   6, & ! element 136
   20 & ! element 21
  /)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_Nips = &
 (/8, & ! element 7
   1, & ! element 134
   4, & ! element 11
   9, & ! element 27
   4, & ! element 157
   6, & ! element 136
   27 & ! element 21
  /)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_NipNeighbors = &
 (/6, & ! element 7
   4, & ! element 134
   4, & ! element 11
   4, & ! element 27
   6, & ! element 157
   6, & ! element 136
   6  & ! element 21
  /)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_NsubNodes = &
 (/19,& ! element 7
   0, & ! element 134
   5, & ! element 11
   12,& ! element 27
   0, & ! element 157
   15,& ! element 136
   56 & ! element 21   
  /)
 integer(pInt), dimension(FE_maxNipNeighbors,FE_Nelemtypes), parameter :: FE_NfaceNodes = &
 reshape((/&
  4,4,4,4,4,4, & ! element 7
  3,3,3,3,0,0, & ! element 134
  2,2,2,2,0,0, & ! element 11
  2,2,2,2,0,0, & ! element 27
  3,3,3,3,0,0, & ! element 157
  3,4,4,4,3,0, & ! element 136
  4,4,4,4,4,4  & ! element 21
   /),(/FE_maxNipNeighbors,FE_Nelemtypes/))
 integer(pInt), dimension(FE_NipFaceNodes,FE_maxNipNeighbors,FE_Nelemtypes), parameter :: FE_nodeOnFace = &
 reshape((/&
  1,2,3,4 , & ! element 7
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,0 , & ! element 134
  1,4,2,0 , &
  2,3,4,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element 11
  2,3,0,0 , &
  3,4,0,0 , &
  4,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,0,0 , & ! element 27
  2,3,0,0 , &
  3,4,0,0 , &
  4,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 157
  1,4,2,0 , &
  2,3,4,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  1,2,3,0 , & ! element 136
  1,4,5,2 , &
  2,5,6,3 , &
  1,3,6,4 , &
  4,6,5,0 , &
  0,0,0,0 , &
  1,2,3,4 , & ! element 21
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5   &
   /),(/FE_NipFaceNodes,FE_maxNipNeighbors,FE_Nelemtypes/))

 CONTAINS
! ---------------------------
! subroutine mesh_init()
! function mesh_FEtoCPelement(FEid)
! function mesh_build_ipNeighorhood()
! ---------------------------


!***********************************************************
! initialization 
!***********************************************************
 subroutine mesh_init (ip,element)

 use cpfem_interface
 use prec, only: pInt
 use IO, only: IO_error,IO_open_InputFile
 use FEsolving, only: parallelExecution, FEsolving_execElem, FEsolving_execIP, calcMode, lastMode
 
 implicit none
 
 integer(pInt), parameter :: fileUnit = 222
 integer(pInt) e,element,ip
 
 write(6,*)
 write(6,*) '<<<+-  mesh init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)

 call mesh_build_FEdata()                                      ! --- get properties of the different types of elements

 if (IO_open_inputFile(fileUnit)) then                         ! --- parse info from input file...

 select case (FEsolver)
     case ('Spectral')
                         call mesh_spectral_count_nodesAndElements(fileUnit)
                         call mesh_spectral_count_cpElements()
                         call mesh_spectral_map_elements()
                         call mesh_spectral_map_nodes()
                         call mesh_spectral_count_cpSizes()
                         call mesh_spectral_build_nodes(fileUnit)
                         call mesh_spectral_build_elements(fileUnit)
     
     case ('Marc')
                         call mesh_marc_get_tableStyles(fileUnit)
                         call mesh_marc_count_nodesAndElements(fileUnit)
                         call mesh_marc_count_elementSets(fileUnit)
                         call mesh_marc_map_elementSets(fileUnit)
                         call mesh_marc_count_cpElements(fileUnit)
                         call mesh_marc_map_elements(fileUnit)
                         call mesh_marc_map_nodes(fileUnit)
                         call mesh_marc_build_nodes(fileUnit)
                         call mesh_marc_count_cpSizes(fileunit)
                         call mesh_marc_build_elements(fileUnit)
     case ('Abaqus')
                         call mesh_abaqus_count_nodesAndElements(fileUnit)
                         call mesh_abaqus_count_elementSets(fileUnit)
                         call mesh_abaqus_count_materials(fileUnit)
                         call mesh_abaqus_map_elementSets(fileUnit)
                         call mesh_abaqus_map_materials(fileUnit)
                         call mesh_abaqus_count_cpElements(fileUnit)
                         call mesh_abaqus_map_elements(fileUnit)
                         call mesh_abaqus_map_nodes(fileUnit)
                         call mesh_abaqus_build_nodes(fileUnit)
                         call mesh_abaqus_count_cpSizes(fileunit)
                         call mesh_abaqus_build_elements(fileUnit)
   end select
   close (fileUnit)
   
   call mesh_build_sharedElems()
   call mesh_build_ipNeighborhood()
   call mesh_build_subNodeCoords()
   call mesh_build_ipVolumes()
   call mesh_build_ipAreas()
   call mesh_tell_statistics()

   parallelExecution = (parallelExecution .and. (mesh_Nelems == mesh_NcpElems))      ! plus potential killer from non-local constitutive
 else
   call IO_error(101) ! cannot open input file
 endif
 
 FEsolving_execElem = (/1,mesh_NcpElems/)
 allocate(FEsolving_execIP(2,mesh_NcpElems)); FEsolving_execIP = 1_pInt
 forall (e = 1:mesh_NcpElems) FEsolving_execIP(2,e) = FE_Nips(mesh_element(2,e))
 
 allocate(calcMode(mesh_maxNips,mesh_NcpElems))
 calcMode = .false.                                       ! pretend to have collected what first call is asking (F = I)
 calcMode(ip,mesh_FEasCP('elem',element)) = .true.        ! first ip,el needs to be already pingponged to "calc"
 lastMode = .true.                                        ! and its mode is already known...
 endsubroutine
 


!***********************************************************
! mapping of FE element types to internal representation
!***********************************************************
 function FE_mapElemtype(what)

 implicit none
 
 character(len=*), intent(in) :: what
 integer(pInt) FE_mapElemtype
  
 select case (what)
    case (  '7', 'C3D8')
      FE_mapElemtype = 1            ! Three-dimensional Arbitrarily Distorted Brick
    case ('134', 'C3D4')
      FE_mapElemtype = 2            ! Three-dimensional Four-node Tetrahedron
    case ( '11', 'CPE4')
      FE_mapElemtype = 3            ! Arbitrary Quadrilateral Plane-strain
    case ( '27', 'CPE8')
      FE_mapElemtype = 4            ! Plane Strain, Eight-node Distorted Quadrilateral
    case ('157')
      FE_mapElemtype = 5            ! Three-dimensional, Low-order, Tetrahedron, Herrmann Formulations
    case ('136', 'C3D6')
      FE_mapElemtype = 6            ! Three-dimensional Arbitrarily Distorted Pentahedral
    case ( '21', 'C3D20')
      FE_mapElemtype = 7            ! Three-dimensional Arbitrarily Distorted Quadratic Hexahedral
    case default 
      FE_mapElemtype = 0            ! unknown element --> should raise an error upstream..!
 endselect

 endfunction



!***********************************************************
! FE to CP id mapping by binary search thru lookup array
!
! valid questions are 'elem', 'node'
!***********************************************************
 function mesh_FEasCP(what,id)
 
 use prec, only: pInt
 use IO, only: IO_lc
 implicit none
 
 character(len=*), intent(in) :: what
 integer(pInt), intent(in) :: id
 integer(pInt), dimension(:,:), pointer :: lookupMap
 integer(pInt) mesh_FEasCP, lower,upper,center
 
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
 upper = size(lookupMap,2)
 
 ! check at bounds QUESTION is it valid to extend bounds by 1 and just do binary search w/o init check at bounds?
 if (lookupMap(1,lower) == id) then
   mesh_FEasCP = lookupMap(2,lower)
   return
 elseif (lookupMap(1,upper) == id) then
   mesh_FEasCP = lookupMap(2,upper)
   return
 endif
 
 ! binary search in between bounds
 do while (upper-lower > 1)
   center = (lower+upper)/2
   if (lookupMap(1,center) < id) then
     lower = center
   elseif (lookupMap(1,center) > id) then
     upper = center
   else
     mesh_FEasCP = lookupMap(2,center)
     exit
   endif
 enddo
 return
 
 endfunction


!***********************************************************
! find face-matching element of same type
!***********************************************************
 function mesh_faceMatch(face,elem)

 use prec, only: pInt
 implicit none

 integer(pInt) face,elem
 integer(pInt) mesh_faceMatch
 integer(pInt), dimension(FE_NfaceNodes(face,mesh_element(2,elem))) :: nodeMap
 integer(pInt) minN,NsharedElems,lonelyNode,faceNode,i,n,t
 
 minN = mesh_maxNsharedElems+1         ! init to worst case
 mesh_faceMatch = 0_pInt               ! intialize to "no match found"
 t = mesh_element(2,elem)              ! figure elemType

 do faceNode=1,FE_NfaceNodes(face,t)   ! loop over nodes on face
   nodeMap(faceNode) = mesh_FEasCP('node',mesh_element(4+FE_nodeOnFace(faceNode,face,t),elem)) ! CP id of face node
   NsharedElems = mesh_sharedElem(1,nodeMap(faceNode)) ! figure # shared elements for this node
   if (NsharedElems < minN) then
     minN = NsharedElems               ! remember min # shared elems
     lonelyNode = faceNode             ! remember most lonely node
   endif
 enddo
candidate: do i=1,minN  ! iterate over lonelyNode's shared elements
   mesh_faceMatch = mesh_sharedElem(1+i,nodeMap(lonelyNode)) ! present candidate elem
   if (mesh_faceMatch == elem) then    ! my own element ?
     mesh_faceMatch = 0_pInt           ! disregard
     cycle candidate
   endif
   do faceNode=1,FE_NfaceNodes(face,t) ! check remaining face nodes to match
     if (faceNode == lonelyNode) cycle ! disregard lonely node (matches anyway)
     n = nodeMap(faceNode)
     if (all(mesh_sharedElem(2:1+mesh_sharedElem(1,n),n) /= mesh_faceMatch)) then  ! no ref to candidate elem?
       mesh_faceMatch = 0_pInt         ! set to "no match" (so far)
       cycle candidate                 ! next candidate elem
     endif
   enddo
   exit                                ! surviving candidate
 enddo candidate
  
 return
 
 endfunction

 
!********************************************************************
! get properties of different types of finite elements
!
! assign globals:
! FE_nodesAtIP, FE_ipNeighbor, FE_subNodeParent, FE_subNodeOnIPFace
!********************************************************************
 subroutine mesh_build_FEdata ()
 
 use prec, only: pInt
 implicit none
 
 allocate(FE_nodesAtIP(2,2,FE_maxNips,FE_Nelemtypes)) ; FE_nodesAtIP = 0_pInt
 allocate(FE_ipNeighbor(FE_maxNipNeighbors,FE_maxNips,FE_Nelemtypes)) ; FE_ipNeighbor = 0_pInt
 allocate(FE_subNodeParent(FE_maxNips,FE_maxNsubNodes,FE_Nelemtypes)) ; FE_subNodeParent = 0_pInt
 allocate(FE_subNodeOnIPFace(FE_NipFaceNodes,FE_maxNipNeighbors,FE_maxNips,FE_Nelemtypes)) ; FE_subNodeOnIPFace = 0_pInt
 
 ! fill FE_nodesAtIP with data
 FE_nodesAtIP(:,:,:FE_Nips(1),1) = &  ! element 7
    reshape((/&
    1,0, 0,0,  &
    2,0, 0,0,  &
    4,0, 0,0,  &
    3,0, 0,0,  &
    5,0, 0,0,  &
    6,0, 0,0,  &
    8,0, 0,0,  &
    7,0, 0,0   &
    /),(/2,2,FE_Nips(1)/))
 FE_nodesAtIP(:,:,:FE_Nips(2),2) = &  ! element 134
    reshape((/&
    1,0, 0,0   &
    /),(/2,2,FE_Nips(2)/))
 FE_nodesAtIP(:,:,:FE_Nips(3),3) = &  ! element 11
    reshape((/&
    1,0, 0,0,  &
    2,0, 0,0,  &
    4,0, 0,0,  &
    3,0, 0,0   &
    /),(/2,2,FE_Nips(3)/))
 FE_nodesAtIP(:,:,:FE_Nips(4),4) = &  ! element 27
    reshape((/&
    1,0, 0,0,  &
    1,2, 0,0,  &
    2,0, 0,0,  &
    1,4, 0,0,  &
    1,3, 2,4,  &
    2,3, 0,0,  &
    4,0, 0,0,  &
    3,4, 0,0,  &
    3,0, 0,0   &
    /),(/2,2,FE_Nips(4)/))
 FE_nodesAtIP(:,:,:FE_Nips(5),5) = &  ! element 157
    reshape((/&
    1,0, 0,0,  &
    2,0, 0,0,  &
    3,0, 0,0,  &
    4,0, 0,0   &
    /),(/2,2,FE_Nips(5)/))
 FE_nodesAtIP(:,:,:FE_Nips(6),6) = &  ! element 136
    reshape((/&
    1,0, 0,0,  &
    2,0, 0,0,  &
    3,0, 0,0,  &
    4,0, 0,0,  &
    5,0, 0,0,  &
    6,0, 0,0   &
    /),(/2,2,FE_Nips(6)/))
  FE_nodesAtIP(:,:,:FE_Nips(7),7) = &  ! element 21
    reshape((/&
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
    /),(/2,2,FE_Nips(7)/))
 
 ! fill FE_ipNeighbor with data
 FE_ipNeighbor(:FE_NipNeighbors(1),:FE_Nips(1),1) = &  ! element 7
    reshape((/&
     2,-5, 3,-2, 5,-1,  &
    -3, 1, 4,-2, 6,-1,  &
     4,-5,-4, 1, 7,-1,  &
    -3, 3,-4, 2, 8,-1,  &
     6,-5, 7,-2,-6, 1,  &
    -3, 5, 8,-2,-6, 2,  &
     8,-5,-4, 5,-6, 3,  &
    -3, 7,-4, 6,-6, 4   &
    /),(/FE_NipNeighbors(1),FE_Nips(1)/))
 FE_ipNeighbor(:FE_NipNeighbors(2),:FE_Nips(2),2) = &  ! element 134
    reshape((/&
    -1,-2,-3,-4   &
    /),(/FE_NipNeighbors(2),FE_Nips(2)/))
 FE_ipNeighbor(:FE_NipNeighbors(3),:FE_Nips(3),3) = &  ! element 11
    reshape((/&
     2,-4, 3,-1,  &
    -2, 1, 4,-1,  &
     4,-4,-3, 1,  &
    -2, 3,-3, 2   &
    /),(/FE_NipNeighbors(3),FE_Nips(3)/))
 FE_ipNeighbor(:FE_NipNeighbors(4),:FE_Nips(4),4) = &  ! element 27
    reshape((/&
     2,-4, 4,-1,  &
     3, 1, 5,-1,  &
    -2, 2, 6,-1,  &
     5,-4, 7, 1,  &
     6, 4, 8, 2,  &
    -2, 5, 9, 3,  &
     8,-4,-3, 4,  &
     9, 7,-3, 5,  &
    -2, 8,-3, 6   &
    /),(/FE_NipNeighbors(4),FE_Nips(4)/))
 FE_ipNeighbor(:FE_NipNeighbors(5),:FE_Nips(5),5) = &  ! element 157
    reshape((/&
     2,-4, 3,-2, 4,-1,  &
     3,-2, 1,-3, 4,-1,  &
     1,-3, 2,-4, 4,-1,  &
     1,-3, 2,-4, 3,-2   &
    /),(/FE_NipNeighbors(5),FE_Nips(5)/))
 FE_ipNeighbor(:FE_NipNeighbors(6),:FE_Nips(6),6) = &  ! element 136
    reshape((/&
     2,-4, 3,-2, 4,-1,  &
    -3, 1, 3,-2, 5,-1,  &
     2,-4,-3, 1, 6,-1,  &
     5,-4, 6,-2,-5, 1,  &
    -3, 4, 6,-2,-5, 2,  &
     5,-4,-3, 4,-5, 3   &
    /),(/FE_NipNeighbors(6),FE_Nips(6)/))
 FE_ipNeighbor(:FE_NipNeighbors(7),:FE_Nips(7),7) = &  ! element 21
    reshape((/&
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
    /),(/FE_NipNeighbors(7),FE_Nips(7)/))
 
 ! fill FE_subNodeParent with data
 FE_subNodeParent(:FE_Nips(1),:FE_NsubNodes(1),1) = &  ! element 7
    reshape((/&
    1, 2, 0, 0, 0, 0, 0, 0,  &
    2, 3, 0, 0, 0, 0, 0, 0,  &
    3, 4, 0, 0, 0, 0, 0, 0,  &
    4, 1, 0, 0, 0, 0, 0, 0,  &
    1, 5, 0, 0, 0, 0, 0, 0,  &
    2, 6, 0, 0, 0, 0, 0, 0,  &
    3, 7, 0, 0, 0, 0, 0, 0,  &
    4, 8, 0, 0, 0, 0, 0, 0,  &
    5, 6, 0, 0, 0, 0, 0, 0,  &
    6, 7, 0, 0, 0, 0, 0, 0,  &
    7, 8, 0, 0, 0, 0, 0, 0,  &
    8, 5, 0, 0, 0, 0, 0, 0,  &
    1, 2, 3, 4, 0, 0, 0, 0,  &
    1, 2, 6, 5, 0, 0, 0, 0,  &
    2, 3, 7, 6, 0, 0, 0, 0,  &
    3, 4, 8, 7, 0, 0, 0, 0,  &
    1, 4, 8, 5, 0, 0, 0, 0,  &
    5, 6, 7, 8, 0, 0, 0, 0,  &
    1, 2, 3, 4, 5, 6, 7, 8   & 
    /),(/FE_Nips(1),FE_NsubNodes(1)/))
 !FE_subNodeParent(:FE_Nips(2),:FE_NsubNodes(2),2) = &  ! element 134
 !   reshape((/&
 !   *still to be defined*
 !   /),(/FE_Nips(2),FE_NsubNodes(2)/))
 FE_subNodeParent(:FE_Nips(3),:FE_NsubNodes(3),3) = &  ! element 11
    reshape((/&
    1, 2, 0, 0,  & 
    2, 3, 0, 0,  & 
    3, 4, 0, 0,  &
    4, 1, 0, 0,  & 
    1, 2, 3, 4   &
    /),(/FE_Nips(3),FE_NsubNodes(3)/))
 FE_subNodeParent(:FE_Nips(4),:FE_NsubNodes(4),4) = &  ! element 27
    reshape((/&
    1, 1, 2, 0, 0, 0, 0, 0, 0,  &
    1, 2, 2, 0, 0, 0, 0, 0, 0,  & 
    2, 2, 3, 0, 0, 0, 0, 0, 0,  & 
    2, 3, 3, 0, 0, 0, 0, 0, 0,  & 
    3, 3, 4, 0, 0, 0, 0, 0, 0,  & 
    3, 4, 4, 0, 0, 0, 0, 0, 0,  & 
    4, 4, 1, 0, 0, 0, 0, 0, 0,  & 
    4, 1, 1, 0, 0, 0, 0, 0, 0,  & 
    1, 1, 1, 1, 2, 2, 4, 4, 3,  & 
    2, 2, 2, 2, 1, 1, 3, 3, 4,  & 
    3, 3, 3, 3, 2, 2, 4, 4, 1,  & 
    4, 4, 4, 4, 1, 1, 3, 3, 2   &
    /),(/FE_Nips(4),FE_NsubNodes(4)/))
 !FE_subNodeParent(:FE_Nips(5),:FE_NsubNodes(5),5) = &  ! element 157
 !   reshape((/&
 !   *still to be defined*
 !   /),(/FE_Nips(5),FE_NsubNodes(5)/))
 FE_subNodeParent(:FE_Nips(6),:FE_NsubNodes(6),6) = &  ! element 136
    reshape((/&
    1, 2, 0, 0, 0, 0,  &
    2, 3, 0, 0, 0, 0,  &
    3, 1, 0, 0, 0, 0,  &
    1, 4, 0, 0, 0, 0,  &
    2, 5, 0, 0, 0, 0,  &
    3, 6, 0, 0, 0, 0,  &
    4, 5, 0, 0, 0, 0,  &
    5, 6, 0, 0, 0, 0,  &
    6, 4, 0, 0, 0, 0,  &
    1, 2, 3, 0, 0, 0,  &
    1, 2, 4, 5, 0, 0,  &
    2, 3, 5, 6, 0, 0,  &
    1, 3, 4, 6, 0, 0,  &
    4, 5, 6, 0, 0, 0,  &
    1, 2, 3, 4, 5, 6   &
    /),(/FE_Nips(6),FE_NsubNodes(6)/))
 FE_subNodeParent(:FE_Nips(7),:FE_NsubNodes(7),7) = &  ! element 21
    reshape((/&
    1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    1, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 4, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    1, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    2, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    3, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    4, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    5, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    5, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & 
    6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 2, 2, 4, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 2, 2, 2, 1, 1, 3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    3, 3, 3, 3, 2, 2, 4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    4, 4, 4, 4, 1, 1, 3, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 2, 2, 5, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 2, 2, 2, 1, 1, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 2, 2, 2, 3, 3, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    3, 3, 3, 3, 2, 2, 7, 7, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    3, 3, 3, 3, 4, 4, 7, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    4, 4, 4, 4, 3, 3, 8, 8, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    4, 4, 4, 4, 1, 1, 8, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 4, 4, 5, 5, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    5, 5, 5, 5, 1, 1, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 6, 6, 6, 2, 2, 5, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 6, 6, 6, 2, 2, 7, 7, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 7, 7, 3, 3, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 7, 7, 3, 3, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 8, 8, 4, 4, 7, 7, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 8, 8, 4, 4, 5, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    5, 5, 5, 5, 1, 1, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    5, 5, 5, 5, 6, 6, 8, 8, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    6, 6, 6, 6, 5, 5, 7, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    7, 7, 7, 7, 6, 6, 8, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    8, 8, 8, 8, 5, 5, 7, 7, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 5, 5, 5, 5, 3, 3, 6, 6, 8, 8, 7, &
    2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 3, 6, 6, 6, 6, 4, 4, 5, 5, 7, 7, 8, &
    3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4, 7, 7, 7, 7, 1, 1, 6, 6, 8, 8, 5, &
    4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 3, 3, 3, 3, 8, 8, 8, 8, 2, 2, 5, 5, 7, 7, 6, &
    5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 6, 6, 6, 6, 8, 8, 8, 8, 2, 2, 4, 4, 7, 7, 3, &
    6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 5, 5, 5, 5, 7, 7, 7, 7, 1, 1, 3, 3, 8, 8, 4, &
    7, 7, 7, 7, 7, 7, 7, 7, 3, 3, 3, 3, 6, 6, 6, 6, 8, 8, 8, 8, 2, 2, 4, 4, 5, 5, 1, &
    8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4, 4, 5, 5, 5, 5, 7, 7, 7, 7, 1, 1, 3, 3, 6, 6, 2  & 
    /),(/FE_Nips(7),FE_NsubNodes(7)/))
 
 ! fill FE_subNodeOnIPFace with data
 FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(1),:FE_Nips(1),1) = &  ! element 7
    reshape((/&
     9,21,27,22, & ! 1
     1,13,25,12, &
    12,25,27,21, &
     1, 9,22,13, &
    13,22,27,25, &
     1,12,21, 9, &
     2,10,23,14, & ! 2
     9,22,27,21, &
    10,21,27,23, &
     2,14,22, 9, &
    14,23,27,22, &
     2, 9,21,10, &
    11,24,27,21, & ! 3
     4,12,25,16, &
     4,16,24,11, &
    12,21,27,25, &
    16,25,27,24, &
     4,11,21,12, &
     3,15,23,10, & ! 4
    11,21,27,24, &
     3,11,24,15, &
    10,23,27,21, &
    15,24,27,23, &
     3,10,21,11, &
    17,22,27,26, & ! 5
     5,20,25,13, &
    20,26,27,25, &
     5,13,22,17, &
     5,17,26,20, &
    13,25,27,22, &
     6,14,23,18, & ! 6
    17,26,27,22, &
    18,23,27,26, &
     6,17,22,14, &
     6,18,26,17, &
    14,22,27,23, &
    19,26,27,24, & ! 7
     8,16,25,20, &
     8,19,24,16, &
    20,25,27,26, &
     8,20,26,19, &
    16,24,27,25, &
     7,18,23,15, & ! 8
    19,24,27,26, &
     7,15,24,19, &
    18,26,27,23, &
     7,19,26,18, &
    15,23,27,24  &
    /),(/FE_NipFaceNodes,FE_NipNeighbors(1),FE_Nips(1)/))
 FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(2),:FE_Nips(2),2) = &  ! element 134
    reshape((/&
     1, 1, 3, 2, & ! 1
     1, 1, 2, 4, &
     2, 2, 3, 4, &
     1, 1, 4, 3  &
    /),(/FE_NipFaceNodes,FE_NipNeighbors(2),FE_Nips(2)/))
 FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(3),:FE_Nips(3),3) = &  ! element 11
    reshape((/&
     5, 9, 0, 0, & ! 1
     1, 8, 0, 0, &
     8, 9, 0, 0, &
     1, 5, 0, 0, &
     2, 6, 0, 0, & ! 2
     5, 9, 0, 0, &
     6, 9, 0, 0, &
     2, 5, 0, 0, &
     3, 6, 0, 0, & ! 3
     7, 9, 0, 0, &
     3, 7, 0, 0, &
     6, 9, 0, 0, &
     7, 9, 0, 0, & ! 4
     4, 8, 0, 0, &
     4, 7, 0, 0, &
     8, 9, 0, 0  &
    /),(/FE_NipFaceNodes,FE_NipNeighbors(3),FE_Nips(3)/))
 FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(4),:FE_Nips(4),4) = &  ! element 27
    reshape((/&
     9,17, 0, 0, & ! 1
     1,16, 0, 0, &
    16,17, 0, 0, &
     1, 9, 0, 0, &
    10,18, 0, 0, & ! 2
     9,17, 0, 0, &
    17,18, 0, 0, &
     9,10, 0, 0, &
     2,11, 0, 0, & ! 3
    10,18, 0, 0, &
    11,18, 0, 0, &
     2,10, 0, 0, &
    17,20, 0, 0, & ! 4
    15,16, 0, 0, &
    15,20, 0, 0, &
    16,17, 0, 0, &
    18,19, 0, 0, & ! 5
    17,20, 0, 0, &
    19,20, 0, 0, &
    17,18, 0, 0, &
    11,12, 0, 0, & ! 6
    18,19, 0, 0, &
    12,19, 0, 0, &
    11,18, 0, 0, &
    14,20, 0, 0, & ! 7
     4,15, 0, 0, &
     4,14, 0, 0, &
    15,20, 0, 0, &
    13,19, 0, 0, & ! 8
    14,20, 0, 0, &
    13,14, 0, 0, &
    19,20, 0, 0, &
     3,12, 0, 0, & ! 9
    13,19, 0, 0, &
     3,13, 0, 0, &
    12,19, 0, 0  &
    /),(/FE_NipFaceNodes,FE_NipNeighbors(4),FE_Nips(4)/))
 !FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(5),:FE_Nips(5),5) = &  ! element 157
 !   reshape((/&
 !   *still to be defined*
 !   /),(/FE_NipFaceNodes,FE_NipNeighbors(5),FE_Nips(5)/))
 FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(6),:FE_Nips(6),6) = &  ! element 136
    reshape((/&
     7,16,21,17, & ! 1
     1,10,19, 9, &
     9,19,21,16, &
     1, 7,17,10, &
    10,17,21,19, &
     1, 9,16, 7, &
     2, 8,18,11, & ! 2
     7,17,21,16, &
     8,16,21,18, &
     2,11,17, 7, &
    11,18,21,17, &
     2, 7,16, 8, &
     8,18,21,16, & ! 3
     3, 9,19,12, &
     3,12,18, 8, &
     9,16,21,19, &
    12,19,21,18, &
     3, 8,16, 9, &
    13,17,21,20, & ! 4
     4,15,19,10, &
    15,20,21,19, &
     4,10,17,13, &
     4,13,20,15, &
    10,19,21,17, &
     5,11,18,14, & ! 5
    13,20,21,17, &
    14,18,21,20, &
     5,13,17,11, &
     5,14,20,13, &
    11,17,21,18, &
    14,20,21,18, & ! 6
     6,12,19,15, &
     6,14,18,12, &
    15,19,21,20, &
     6,15,20,14, &
    12,18,21,19  &
    /),(/FE_NipFaceNodes,FE_NipNeighbors(6),FE_Nips(6)/))
 FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(7),:FE_Nips(7),7) = &  ! element 21
    reshape((/&
     9,33,57,37, &  ! 1
     1,17,44,16, &
    33,16,44,57, &
     1, 9,37,17, &
    17,37,57,44, &
     1,16,33, 9, &
    10,34,58,38, &  !  2
     9,37,57,33, &
    34,33,57,58, &
     9,10,38,37, &
    37,38,58,57, &
     9,33,34,10, &
     2,11,39,18, &  !  3
    10,38,58,34, &
    11,34,58,39, &
    10, 2,18,38, &
    38,18,39,58, &
    10,34,11, 2, &
    33,36,60,57, &  !  4
    16,44,43,15, &
    36,15,43,60, &
    16,33,57,44, &
    44,57,60,43, &
    16,15,36,33, &
    34,35,59,58, &  !  5
    33,57,60,36, &
    35,36,60,59, &
    33,34,58,57, &
    57,58,59,60, &
    33,36,35,34, &
    11,12,40,39, &  !  6
    34,58,59,35, &
    12,35,59,40, &
    34,11,39,58, &
    58,39,40,59, &
    34,35,12,11, &
    36,14,42,60, &  !  7
    15,43,20, 4, &
    14, 4,20,42, &
    15,36,60,43, &
    43,60,42,20, &
    15, 4,14,36, &
    35,13,41,59, &  !  8
    36,60,42,14, &
    13,14,42,41, &
    36,35,59,60, &
    60,59,41,42, &
    36,14,13,35, &
    12, 3,19,40, &  !  9
    35,59,41,13, &
     3,13,41,19, &
    35,12,40,59, &
    59,40,19,41, &
    35,13, 3,12, &
    37,57,61,45, &  ! 10
    17,21,52,44, &
    57,44,52,61, &
    17,37,45,21, &
    21,45,61,52, &
    17,44,57,37, &
    38,58,62,46, &  ! 11
    37,45,61,57, &
    58,57,61,62, &
    37,38,46,45, &
    45,46,62,61, &
    37,57,58,38, &
    18,39,47,22, &  ! 12
    38,46,62,58, &
    39,58,62,47, &
    38,18,22,46, &
    46,22,47,62, &
    38,58,39,18, &
    57,60,64,61, &  ! 13
    44,52,51,43, &
    60,43,51,64, &
    44,57,61,52, &
    52,61,64,51, &
    44,43,60,57, &
    58,59,63,62, &  ! 14
    57,61,64,60, &
    59,60,64,63, &
    57,58,62,61, &
    61,62,63,64, &
    57,60,59,58, &
    39,40,48,47, &  ! 15
    58,62,63,59, &
    40,59,63,48, &
    58,39,47,62, &
    62,47,48,63, &
    58,59,40,39, &
    60,42,50,64, &  ! 16
    43,51,24,20, &
    42,20,24,50, &
    43,60,64,51, &
    51,64,50,24, &
    43,20,42,60, &
    59,41,49,63, &  ! 17
    60,64,50,42, &
    41,42,50,49, &
    60,59,63,64, &
    64,63,49,50, &
    60,42,41,59, &
    40,19,23,48, &  ! 18
    59,63,49,41, &
    19,41,49,23, &
    59,40,48,63, &
    63,48,23,49, &
    59,41,19,40, &
    45,61,53,25, &  ! 19
    21, 5,32,52, &
    61,52,32,53, &
    21,45,25, 5, &
     5,25,53,32, &
    21,52,61,45, &
    46,62,54,26, &  ! 20
    45,25,53,61, &
    62,61,53,54, &
    45,46,26,25, &
    25,26,54,53, &
    45,61,62,46, &
    22,47,27, 6, &  ! 21
    46,26,54,62, &
    47,62,54,27, &
    46,22, 6,26, &
    26, 6,27,54, &
    46,62,47,22, &
    61,64,56,53, &  ! 22
    52,32,31,51, &
    64,51,31,56, &
    52,61,53,32, &
    32,53,56,31, &
    52,51,64,61, &
    62,63,55,54, &  ! 23
    61,53,56,64, &
    63,64,56,55, &
    61,62,54,53, &
    53,54,55,56, &
    61,64,63,62, &
    47,48,28,27, &  ! 24
    62,54,55,63, &
    48,63,55,28, &
    62,47,27,54, &
    54,27,28,55, &
    62,63,48,47, &
    64,50,30,56, &  ! 25
    51,31, 8,24, &
    50,24, 8,30, &
    51,64,56,31, &
    31,56,30, 8, &
    51,24,50,64, &
    63,49,29,55, &  ! 26
    64,56,30,50, &
    49,50,30,29, &
    64,63,55,56, &
    56,55,29,30, &
    64,50,49,63, &
    48,23, 7,28, &  ! 27
    63,55,29,49, &
    23,49,29, 7, &
    63,48,28,55, &
    55,28, 7,29, &
    63,49,23,48  &
    /),(/FE_NipFaceNodes,FE_NipNeighbors(7),FE_Nips(7)/))
 
 return
 
 endsubroutine


!********************************************************************
! figure out table styles (Marc only)
!
! initialcondTableStyle, hypoelasticTableStyle
!********************************************************************
 subroutine mesh_marc_get_tableStyles (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 6
 integer(pInt), dimension (1+2*maxNchunks) :: pos

 integer(pInt) unit
 character(len=300) line

 initialcondTableStyle = 0_pInt
 hypoelasticTableStyle = 0_pInt
 
610 FORMAT(A300)

 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_stringValue(line,pos,1)) == 'table' .and. pos(1) .GT. 5) then
     initialcondTableStyle = IO_intValue(line,pos,4)
     hypoelasticTableStyle = IO_intValue(line,pos,5)
     exit
   endif
 enddo

620 return
 
 endsubroutine

 
!********************************************************************
! count overall number of nodes and elements in mesh
!
! mesh_Nelems, mesh_Nnodes
!********************************************************************
 subroutine mesh_spectral_count_nodesAndElements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 integer(pInt) a,b,c,i

 integer(pInt) unit
 character(len=1024) line

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt

 rewind(unit)
 do 
   read(unit,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   pos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,pos,1)) == 'resolution') then
		 do i = 2,6,2
			 select case (IO_lc(IO_stringValue(line,pos,i)))
				 case('a')
						 a = IO_intValue(line,pos,i+1)
				 case('b')
						 b = IO_intValue(line,pos,i+1)
				 case('c')
						 c = IO_intValue(line,pos,i+1)
       end select
     enddo
		 mesh_Nelems = 2**(a+b+c)
		 mesh_Nnodes = (1+2**a)*(1+2**b)*(1+2**c)
     exit
   endif
 enddo

100 return
 
 endsubroutine

!********************************************************************
! count overall number of nodes and elements in mesh
!
! mesh_Nelems, mesh_Nnodes
!********************************************************************
 subroutine mesh_marc_count_nodesAndElements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 4
 integer(pInt), dimension (1+2*maxNchunks) :: pos

 integer(pInt) unit
 character(len=300) line

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt

610 FORMAT(A300)

 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,pos,1)) == 'sizing') then
       mesh_Nelems = IO_IntValue (line,pos,3)
       mesh_Nnodes = IO_IntValue (line,pos,4)
     exit
   endif
 enddo

620 return
 
 endsubroutine

!********************************************************************
! count overall number of nodes and elements in mesh
!
! mesh_Nelems, mesh_Nnodes
!********************************************************************
 subroutine mesh_abaqus_count_nodesAndElements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit
 logical inPart

 mesh_Nnodes = 0
 mesh_Nelems = 0
 
610 FORMAT(A300)

 inPart = .false.
 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.
   
   if (inPart) then
     select case ( IO_lc(IO_stringValue(line,pos,1)))
       case('*node')
          if( &
              IO_lc(IO_stringValue(line,pos,2)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,pos,2)) /= 'print'    .and. &
              IO_lc(IO_stringValue(line,pos,2)) /= 'file'     .and. &
              IO_lc(IO_stringValue(line,pos,2)) /= 'response' &
             ) &
            mesh_Nnodes = mesh_Nnodes + IO_countDataLines(unit)
       case('*element')
          if( &
              IO_lc(IO_stringValue(line,pos,2)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,pos,2)) /= 'matrix'   .and. &
              IO_lc(IO_stringValue(line,pos,2)) /= 'response' &
             ) then
            mesh_Nelems = mesh_Nelems + IO_countDataLines(unit)
          endif
     endselect
   endif
 enddo

620 return
 
 endsubroutine

 
!********************************************************************
! count overall number of element sets in mesh
!
! mesh_NelemSets, mesh_maxNelemInSet
!********************************************************************
 subroutine mesh_marc_count_elementSets (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos

 integer(pInt) unit
 character(len=300) line

 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = 0_pInt

610 FORMAT(A300)

 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,pos,1)) == 'define' .and. &
        IO_lc(IO_StringValue(line,pos,2)) == 'element' ) then
     mesh_NelemSets = mesh_NelemSets + 1_pInt
     mesh_maxNelemInSet = max(mesh_maxNelemInSet, &
                              IO_countContinousIntValues(unit))
   endif
 enddo

620 return
 
 endsubroutine


!********************************************************************
! count overall number of element sets in mesh
!
! mesh_NelemSets, mesh_maxNelemInSet
!********************************************************************
 subroutine mesh_abaqus_count_elementSets (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit
 logical inPart
 
 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = mesh_Nelems               ! have to be conservative, since Abaqus allows for recursive definitons
 
610 FORMAT(A300)

 inPart = .false.
 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.
   
   if ( inPart .and. IO_lc(IO_stringValue(line,pos,1)) == '*elset' ) &
     mesh_NelemSets = mesh_NelemSets + 1_pInt
 enddo

620 return
 
 endsubroutine


!********************************************************************
! count overall number of solid sections sets in mesh (Abaqus only)
!
! mesh_Nmaterials
!********************************************************************
 subroutine mesh_abaqus_count_materials (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit
 logical inPart
 
 mesh_Nmaterials = 0_pInt
 
610 FORMAT(A300)

 inPart = .false.
 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.

   if ( inPart .and. &
        IO_lc(IO_StringValue(line,pos,1)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,pos,2)) == 'section' ) &
     mesh_Nmaterials = mesh_Nmaterials + 1_pInt
 enddo

620 return

 endsubroutine


!********************************************************************
! count overall number of cpElements in mesh
!
! mesh_NcpElems
!********************************************************************
 subroutine mesh_spectral_count_cpElements ()

 use prec, only: pInt
 implicit none

 mesh_NcpElems = mesh_Nelems
 return
 
 endsubroutine
 

!********************************************************************
! count overall number of cpElements in mesh
!
! mesh_NcpElems
!********************************************************************
 subroutine mesh_marc_count_cpElements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension (1+2*maxNchunks) :: pos

 integer(pInt) unit,i
 character(len=300) line

 mesh_NcpElems = 0_pInt

610 FORMAT(A300)

 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_stringValue(line,pos,1)) == 'hypoelastic') then
       do i=1,3+hypoelasticTableStyle  ! Skip 3 or 4 lines
         read (unit,610,END=620) line
       enddo
       mesh_NcpElems = mesh_NcpElems + IO_countContinousIntValues(unit)
     exit
   endif
 enddo

620 return
 
 endsubroutine
 

!********************************************************************
! count overall number of cpElements in mesh
!
! mesh_NcpElems
!********************************************************************
 subroutine mesh_abaqus_count_cpElements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i
 logical materialFound
 character (len=64) materialName,elemSetName
 
 mesh_NcpElems = 0
 
610 FORMAT(A300)

 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_stringValue(line,pos,1)) )
     case('*material')
       materialName = IO_extractValue(IO_lc(IO_stringValue(line,pos,2)),'name')    ! extract name=value
       materialFound = materialName /= ''                                   ! valid name?
     case('*user')
       if (IO_lc(IO_StringValue(line,pos,2)) == 'material' .and. materialFound) then
         do i = 1,mesh_Nmaterials                                           ! look thru material names
           if (materialName == mesh_nameMaterial(i)) exit                   ! found one
         enddo
         if (i <= mesh_Nmaterials) then                                     ! found one?
           elemSetName = mesh_mapMaterial(i)                                ! take corresponding elemSet
           do i = 1,mesh_NelemSets                                          ! look thru all elemSet definitions
             if (elemSetName == mesh_nameElemSet(i)) &                      ! matched?
               mesh_NcpElems = mesh_NcpElems + mesh_mapElemSet(1,i)         ! add those elem count
           enddo
         endif
         materialFound = .false.
       endif
   endselect
 enddo

620 return
 
 endsubroutine


!********************************************************************
! map element sets
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!********************************************************************
 subroutine mesh_marc_map_elementSets (unit)

 use prec, only: pInt
 use IO

 implicit none

 integer(pInt), parameter :: maxNchunks = 4
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,elemSet

 allocate (mesh_nameElemSet(mesh_NelemSets))                     ; mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1+mesh_maxNelemInSet,mesh_NelemSets)) ; mesh_mapElemSet = 0_pInt

610 FORMAT(A300)

 elemSet = 0_pInt
 rewind(unit)
 do
   read (unit,610,END=640) line
   pos = IO_stringPos(line,maxNchunks)
   if( (IO_lc(IO_stringValue(line,pos,1)) == 'define' ) .and. &
       (IO_lc(IO_stringValue(line,pos,2)) == 'element' ) ) then
      elemSet = elemSet+1
      mesh_nameElemSet(elemSet) = IO_stringValue(line,pos,4)
      mesh_mapElemSet(:,elemSet) = IO_continousIntValues(unit,mesh_maxNelemInSet,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
   endif
 enddo
 
640 return

 endsubroutine


!********************************************************************
! Build element set mapping 
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!********************************************************************
 subroutine mesh_abaqus_map_elementSets (unit)

 use prec, only: pInt
 use IO

 implicit none

 integer(pInt), parameter :: maxNchunks = 4
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,elemSet
 logical inPart

 allocate (mesh_nameElemSet(mesh_NelemSets))                     ; mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1+mesh_maxNelemInSet,mesh_NelemSets)) ; mesh_mapElemSet = 0_pInt

610 FORMAT(A300)

 elemSet = 0_pInt
 inPart = .false.
 rewind(unit)
 do
   read (unit,610,END=640) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.
   
   if ( inPart .and. IO_lc(IO_stringValue(line,pos,1)) == '*elset' ) then
     elemSet = elemSet + 1_pInt
     mesh_nameElemSet(elemSet)  = IO_extractValue(IO_lc(IO_stringValue(line,pos,2)),'elset')
     mesh_mapElemSet(:,elemSet) = IO_continousIntValues(unit,mesh_Nelems,mesh_nameElemSet,mesh_mapElemSet,elemSet-1)
   endif
 enddo

640 return

 endsubroutine


!********************************************************************
! map solid section (Abaqus only)
!
! allocate globals: mesh_nameMaterial, mesh_mapMaterial
!********************************************************************
 subroutine mesh_abaqus_map_materials (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 20
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,count
 logical inPart,materialFound
 character(len=64) elemSetName,materialName
 
 allocate (mesh_nameMaterial(mesh_Nmaterials)) ; mesh_nameMaterial = ''
 allocate (mesh_mapMaterial(mesh_Nmaterials)) ;  mesh_mapMaterial = ''

610 FORMAT(A300)

 count = 0
 inPart = .false.
 rewind(unit)
 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.

   if ( inPart .and. &
        IO_lc(IO_StringValue(line,pos,1)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,pos,2)) == 'section' ) then

     elemSetName = ''
     materialName = ''

     do i = 3,pos(1)
       if (IO_extractValue(IO_lc(IO_stringValue(line,pos,i)),'elset') /= '') &
         elemSetName = IO_extractValue(IO_lc(IO_stringValue(line,pos,i)),'elset')
       if (IO_extractValue(IO_lc(IO_stringValue(line,pos,i)),'material') /= '') &
         materialName = IO_extractValue(IO_lc(IO_stringValue(line,pos,i)),'material')
     enddo

     if (elemSetName /= '' .and. materialName /= '') then
       count = count + 1_pInt
       mesh_nameMaterial(count) = materialName         ! name of material used for this section
       mesh_mapMaterial(count)  = elemSetName          ! mapped to respective element set
     endif       
   endif
 enddo

620 return

 endsubroutine



!********************************************************************
! map nodes from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPnode
!********************************************************************
 subroutine mesh_spectral_map_nodes ()

 use prec, only: pInt

 implicit none
 integer(pInt) i

 allocate (mesh_mapFEtoCPnode(2,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

 forall (i = 1:mesh_Nnodes) &
   mesh_mapFEtoCPnode(:,i) = i

 return
 
 endsubroutine



!********************************************************************
! map nodes from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPnode
!********************************************************************
 subroutine mesh_marc_map_nodes (unit)

 use prec, only: pInt
 use math, only: qsort
 use IO

 implicit none

 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt) unit,i

 allocate (mesh_mapFEtoCPnode(2,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

610 FORMAT(A300)

 node_count = 0_pInt

 rewind(unit)
 do
   read (unit,610,END=650) line
   pos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'coordinates' ) then
     read (unit,610,END=650) line                                         ! skip crap line
     do i = 1,mesh_Nnodes
       read (unit,610,END=650) line
       mesh_mapFEtoCPnode(1,i) = IO_fixedIntValue (line,(/0,10/),1)
       mesh_mapFEtoCPnode(2,i) = i
     enddo
     exit
   endif
 enddo

650 call qsort(mesh_mapFEtoCPnode,1,size(mesh_mapFEtoCPnode,2))

 return
 
 endsubroutine



!********************************************************************
! map nodes from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPnode
!********************************************************************
 subroutine mesh_abaqus_map_nodes (unit)

 use prec, only: pInt
 use math, only: qsort
 use IO

 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,count,cpNode
 logical inPart

 allocate (mesh_mapFEtoCPnode(2,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

610 FORMAT(A300)

 cpNode = 0_pInt
 inPart = .false.
 rewind(unit)
 do
   read (unit,610,END=650) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.

   if( inPart .and. &
       IO_lc(IO_stringValue(line,pos,1)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,pos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'response' ) &
   ) then
     count = IO_countDataLines(unit)
     do i = 1,count
       backspace(unit)
     enddo
     do i = 1,count
       read (unit,610,END=650) line
       pos = IO_stringPos(line,maxNchunks)
       cpNode = cpNode + 1_pInt
       mesh_mapFEtoCPnode(1,cpNode) = IO_intValue(line,pos,1)
       mesh_mapFEtoCPnode(2,cpNode) = cpNode
     enddo
   endif
 enddo

650 call qsort(mesh_mapFEtoCPnode,1,size(mesh_mapFEtoCPnode,2))

 return

 endsubroutine


!********************************************************************
! map elements from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPelem
!********************************************************************
 subroutine mesh_spectral_map_elements ()

 use prec, only: pInt

 implicit none
 integer(pInt) i

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

 forall (i = 1:mesh_NcpElems) &
   mesh_mapFEtoCPelem(:,i) = i

 return
 
 endsubroutine



!********************************************************************
! map elements from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPelem
!********************************************************************
 subroutine mesh_marc_map_elements (unit)

 use prec, only: pInt
 use math, only: qsort
 use IO

 implicit none

 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt), dimension (1+mesh_NcpElems) :: contInts
 integer(pInt) unit,i,cpElem

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

610 FORMAT(A300)

 cpElem = 0_pInt
 rewind(unit)
 do
   read (unit,610,END=660) line
   pos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'hypoelastic' ) then
     do i=1,3+hypoelasticTableStyle                                       ! skip three (or four if new table style!) lines
       read (unit,610,END=660) line 
     enddo
     contInts = IO_continousIntValues(unit,mesh_NcpElems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
     do i = 1,contInts(1)
       cpElem = cpElem+1
       mesh_mapFEtoCPelem(1,cpElem) = contInts(1+i)
       mesh_mapFEtoCPelem(2,cpElem) = cpElem
     enddo
   endif
 enddo

660 call qsort(mesh_mapFEtoCPelem,1,size(mesh_mapFEtoCPelem,2))           ! should be mesh_NcpElems

 return
 endsubroutine


!********************************************************************
! map elements from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPelem
!********************************************************************
 subroutine mesh_abaqus_map_elements (unit)

 use prec, only: pInt
 use math, only: qsort
 use IO

 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,j,cpElem
 logical materialFound
 character (len=64) materialName,elemSetName

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

610 FORMAT(A300)

 cpElem = 0_pInt
 rewind(unit)
 do 
   read (unit,610,END=660) line
   pos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_stringValue(line,pos,1)) )
     case('*material')
       materialName = IO_extractValue(IO_lc(IO_stringValue(line,pos,2)),'name')    ! extract name=value
       materialFound = materialName /= ''                                   ! valid name?
     case('*user')
       if (IO_lc(IO_stringValue(line,pos,2)) == 'material' .and. materialFound) then
         do i = 1,mesh_Nmaterials                                           ! look thru material names
           if (materialName == mesh_nameMaterial(i)) exit                   ! found one
         enddo
         if (i <= mesh_Nmaterials) then                                     ! found one?
           elemSetName = mesh_mapMaterial(i)                                ! take corresponding elemSet
           do i = 1,mesh_NelemSets                                          ! look thru all elemSet definitions
             if (elemSetName == mesh_nameElemSet(i)) then                   ! matched?
               do j = 1,mesh_mapElemSet(1,i)
                 cpElem = cpElem + 1_pInt
                 mesh_mapFEtoCPelem(1,cpElem) = mesh_mapElemSet(1+j,i)      ! store FE id
                 mesh_mapFEtoCPelem(2,cpElem) = cpElem                      ! store our id
               enddo
             endif
           enddo
         endif
         materialFound = .false.
       endif
   endselect
 enddo

660 call qsort(mesh_mapFEtoCPelem,1,size(mesh_mapFEtoCPelem,2))             ! should be mesh_NcpElems

 return
 endsubroutine


!********************************************************************
! get maximum count of nodes, IPs, IP neighbors, and subNodes
! among cpElements
!
! _maxNnodes, _maxNips, _maxNipNeighbors, _maxNsubNodes
!********************************************************************
subroutine mesh_spectral_count_cpSizes ()
 
 use prec, only: pInt
 implicit none

 integer(pInt) t
 
 t = FE_mapElemtype('C3D8R')                   ! fake 3D hexahedral 8 node 1 IP element

 mesh_maxNnodes =       FE_Nnodes(t)
 mesh_maxNips =         FE_Nips(t)
 mesh_maxNipNeighbors = FE_NipNeighbors(t)
 mesh_maxNsubNodes =    FE_NsubNodes(t)
 
 endsubroutine


!********************************************************************
! get maximum count of nodes, IPs, IP neighbors, and subNodes
! among cpElements
!
! _maxNnodes, _maxNips, _maxNipNeighbors, _maxNsubNodes
!********************************************************************
subroutine mesh_marc_count_cpSizes (unit)
 
 use prec, only: pInt
 use IO
 implicit none
 
 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,t,e

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNsubNodes    = 0_pInt
 
610 FORMAT(A300)
 rewind(unit)
 do
   read (unit,610,END=630) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'connectivity' ) then
     read (unit,610,END=630) line  ! Garbage line
     do i=1,mesh_Nelems            ! read all elements
       read (unit,610,END=630) line
       pos = IO_stringPos(line,maxNchunks)  ! limit to id and type
       e = mesh_FEasCP('elem',IO_intValue(line,pos,1))
       if (e /= 0) then
         t = FE_mapElemtype(IO_stringValue(line,pos,2))
         mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(t))
         mesh_maxNips =         max(mesh_maxNips,FE_Nips(t))
         mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(t))
         mesh_maxNsubNodes =    max(mesh_maxNsubNodes,FE_NsubNodes(t))
         call IO_skipChunks(unit,FE_NoriginalNodes(t)-(pos(1)-2))        ! read on if FE_Nnodes exceeds node count present on current line
       endif
     enddo
     exit
   endif
 enddo
 
630 return
 endsubroutine


!********************************************************************
! get maximum count of nodes, IPs, IP neighbors, and subNodes
! among cpElements
!
! _maxNnodes, _maxNips, _maxNipNeighbors, _maxNsubNodes
!********************************************************************
 subroutine mesh_abaqus_count_cpSizes (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,count,t
 logical inPart

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNsubNodes    = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(unit)
 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,2)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.

   if( inPart .and. &
       IO_lc(IO_stringValue(line,pos,1)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,pos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_stringValue(line,pos,2),'type'))  ! remember elem type
     count = IO_countDataLines(unit)
     do i = 1,count
       backspace(unit)
     enddo
     do i = 1,count
       read (unit,610,END=620) line
       pos = IO_stringPos(line,maxNchunks)                                  ! limit to 64 nodes max
       if (mesh_FEasCP('elem',IO_intValue(line,pos,1)) /= 0) then                                                     ! disregard non CP elems
         mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(t))
         mesh_maxNips =         max(mesh_maxNips,FE_Nips(t))
         mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(t))
         mesh_maxNsubNodes =    max(mesh_maxNsubNodes,FE_NsubNodes(t))
       endif
     enddo
   endif
 enddo
 
620 return

 endsubroutine


!********************************************************************
! store x,y,z coordinates of all nodes in mesh
!
! allocate globals:
! _node
!********************************************************************
 subroutine mesh_spectral_build_nodes (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 integer(pInt) a,b,c,n,i
 real(pReal)   x,y,z
 logical gotResolution,gotDimension

 integer(pInt) unit
 character(len=64) tag
 character(len=1024) line

 a = 1_pInt
 b = 1_pInt
 c = 1_pInt
 x = 1.0_pReal 
 y = 1.0_pReal 
 z = 1.0_pReal 

 gotResolution = .false.
 gotDimension =  .false.

 rewind(unit)
 do 
   read(unit,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   pos = IO_stringPos(line,maxNchunks)

   select case ( IO_lc(IO_StringValue(line,pos,1)) )
     case ('resolution')
         gotResolution = .true.
				 do i = 2,6,2
					 tag = IO_lc(IO_stringValue(line,pos,i))
					 select case (tag)
						 case('a')
								 a = 1+2**IO_intValue(line,pos,i+1)
						 case('b')
								 b = 1+2**IO_intValue(line,pos,i+1)
						 case('c')
								 c = 1+2**IO_intValue(line,pos,i+1)
					 end select
				 enddo
     case ('dimension')
         gotDimension = .true.
				 do i = 2,6,2
					 tag = IO_lc(IO_stringValue(line,pos,i))
					 select case (tag)
						 case('x')
								 x = IO_floatValue(line,pos,i+1)
						 case('y')
								 y = IO_floatValue(line,pos,i+1)
						 case('z')
								 z = IO_floatValue(line,pos,i+1)
					 end select
				 enddo
	 end select
   if (gotDimension .and. gotResolution) exit
 enddo

! --- sanity checks ---

 if (.not. gotDimension .or. .not. gotResolution) call IO_error(102)
 if (a < 2 .or. b < 2 .or. c < 2) call IO_error(103)
 if (x <= 0.0_pReal .or. y <= 0.0_pReal .or. z <= 0.0_pReal) call IO_error(104)
 
 forall (n = 0:mesh_Nnodes-1)
   mesh_node(1,n+1) = x * dble(mod(n,a)     / (a-1.0_pReal))
   mesh_node(2,n+1) = y * dble(mod(n/a,b)   / (b-1.0_pReal))
   mesh_node(3,n+1) = z * dble(mod(n/a/b,c) / (c-1.0_pReal))
 end forall

100 return

 endsubroutine


!********************************************************************
! store x,y,z coordinates of all nodes in mesh
!
! allocate globals:
! _node
!********************************************************************
 subroutine mesh_marc_build_nodes (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), dimension(5), parameter :: node_ends = (/0,10,30,50,70/)
 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,j,m

 allocate ( mesh_node (3,mesh_Nnodes) ); mesh_node = 0_pInt

610 FORMAT(A300)

 rewind(unit)
 do
   read (unit,610,END=670) line
   pos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'coordinates' ) then
     read (unit,610,END=670) line                                         ! skip crap line
     do i=1,mesh_Nnodes
       read (unit,610,END=670) line
       m = mesh_FEasCP('node',IO_fixedIntValue(line,node_ends,1))
       forall (j = 1:3) mesh_node(j,m) = IO_fixedNoEFloatValue (line,node_ends,j+1)
     enddo
     exit
   endif
 enddo

670 return

 endsubroutine


!********************************************************************
! store x,y,z coordinates of all nodes in mesh
!
! allocate globals:
! _node
!********************************************************************
 subroutine mesh_abaqus_build_nodes (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 4
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt) unit,i,j,m,count
 logical inPart
 
 allocate ( mesh_node (3,mesh_Nnodes) ); mesh_node = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(unit)
 do
   read (unit,610,END=670) line
   pos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.

   if( inPart .and. &
       IO_lc(IO_stringValue(line,pos,1)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,pos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'response' ) &
   ) then
     count = IO_countDataLines(unit)                          ! how many nodes are defined here?
     do i = 1,count
       backspace(unit)                                        ! rewind to first entry
     enddo
     do i = 1,count
       read (unit,610,END=670) line
       pos = IO_stringPos(line,maxNchunks)
       m = mesh_FEasCP('node',IO_intValue(line,pos,1))
       forall (j=1:3) mesh_node(j,m) = IO_floatValue(line,pos,j+1)
     enddo
   endif
 enddo

670 return

 endsubroutine


!********************************************************************
! store FEid, type, mat, tex, and node list per element
!
! allocate globals:
! _element
!********************************************************************
 subroutine mesh_spectral_build_elements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 integer(pInt) a,b,c,e,i,homog
 logical gotResolution,gotDimension,gotHomogenization

 integer(pInt) unit
 character(len=1024) line

 a = 1_pInt
 b = 1_pInt
 c = 1_pInt
 
 gotResolution =      .false.
 gotDimension =       .false.
 gotHomogenization =  .false.

 rewind(unit)
 do 
   read(unit,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   pos = IO_stringPos(line,maxNchunks)

   select case ( IO_lc(IO_StringValue(line,pos,1)) )
     case ('dimension')
         gotDimension = .true.

     case ('homogenization')
         gotHomogenization = .true.
				 homog = IO_intValue(line,pos,2)

     case ('resolution')
         gotResolution = .true.
				 do i = 2,6,2
					 select case (IO_lc(IO_stringValue(line,pos,i)))
						 case('a')
								 a = 2**IO_intValue(line,pos,i+1)
						 case('b')
								 b = 2**IO_intValue(line,pos,i+1)
						 case('c')
								 c = 2**IO_intValue(line,pos,i+1)
					 end select
				 enddo
	 end select
   if (gotDimension .and. gotHomogenization .and. gotResolution) exit
 enddo

100 allocate (mesh_element (4+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

 e = 0_pInt
 do while (e < mesh_NcpElems)
   read(unit,'(a1024)',END=110) line
   if (IO_isBlank(line)) cycle                             ! skip empty lines
   pos = IO_stringPos(line,1)
  
   e = e+1                                                 ! valid element entry
	 mesh_element ( 1,e) = e                                 ! FE id
	 mesh_element ( 2,e) = FE_mapElemtype('C3D8R')           ! elem type
	 mesh_element ( 3,e) = homog                             ! homogenization
	 mesh_element ( 4,e) = IO_IntValue(line,pos,1)           ! microstructure
	 mesh_element ( 5,e) = (e-1) + (e-1)/a + (e-1)/a/b*(a+1) ! base node
	 mesh_element ( 6,e) = mesh_element ( 5,e) + 1
	 mesh_element ( 7,e) = mesh_element ( 5,e) + (a+1) + 1
	 mesh_element ( 8,e) = mesh_element ( 5,e) + (a+1)
	 mesh_element ( 9,e) = mesh_element ( 5,e) + (a+1)*(b+1) ! second floor base node
	 mesh_element (10,e) = mesh_element ( 9,e) + 1
	 mesh_element (11,e) = mesh_element ( 9,e) + (a+1) + 1
	 mesh_element (12,e) = mesh_element ( 9,e) + (a+1)
 enddo
 
110 return

 endsubroutine



!********************************************************************
! store FEid, type, mat, tex, and node list per element
!
! allocate globals:
! _element
!********************************************************************
 subroutine mesh_marc_build_elements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 66
 integer(pInt), dimension (1+2*maxNchunks) :: pos
 character(len=300) line

 integer(pInt), dimension(1+mesh_NcpElems) :: contInts
 integer(pInt) unit,i,j,sv,val,e

 allocate (mesh_element (4+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

610 FORMAT(A300)

 rewind(unit)
 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'connectivity' ) then
     read (unit,610,END=620) line  ! Garbage line
     do i = 1,mesh_Nelems
       read (unit,610,END=620) line
       pos = IO_stringPos(line,maxNchunks)  ! limit to 64 nodes max (plus ID, type)
       e = mesh_FEasCP('elem',IO_intValue(line,pos,1))
       if (e /= 0) then       ! disregard non CP elems
         mesh_element (1,e) = IO_IntValue (line,pos,1)                     ! FE id
         mesh_element (2,e) = FE_mapElemtype(IO_StringValue(line,pos,2))   ! elem type
           forall (j = 1:FE_Nnodes(mesh_element(2,e))) &
             mesh_element(j+4,e) = IO_IntValue(line,pos,j+2)                 ! copy FE ids of nodes
           call IO_skipChunks(unit,FE_NoriginalNodes(mesh_element(2,e))-(pos(1)-2))        ! read on if FE_Nnodes exceeds node count present on current line
       endif
     enddo
     exit
   endif
 enddo
 
620 rewind(unit)                                     ! just in case "initial state" apears before "connectivity"
 read (unit,610,END=620) line
 do
   pos = IO_stringPos(line,2)
   if( (IO_lc(IO_stringValue(line,pos,1)) == 'initial') .and. &
       (IO_lc(IO_stringValue(line,pos,2)) == 'state') ) then
     if (initialcondTableStyle == 2) read (unit,610,END=620) line          ! read extra line for new style     
     read (unit,610,END=630) line                                          ! read line with index of state var
     pos = IO_stringPos(line,1)
     sv = IO_IntValue(line,pos,1)                                          ! figure state variable index
     if( (sv == 2).or.(sv == 3) ) then                                     ! only state vars 2 and 3 of interest
       read (unit,610,END=620) line                                        ! read line with value of state var
       pos = IO_stringPos(line,1)
       do while (scan(IO_stringValue(line,pos,1),'+-',back=.true.)>1)      ! is noEfloat value?
         val = NINT(IO_fixedNoEFloatValue(line,(/0,20/),1))                ! state var's value
         mesh_maxValStateVar(sv-1) = max(val,mesh_maxValStateVar(sv-1))    ! remember max val of homogenization and microstructure index
         if (initialcondTableStyle == 2) then
           read (unit,610,END=630) line                                    ! read extra line     
           read (unit,610,END=630) line                                    ! read extra line     
         endif
         contInts = IO_continousIntValues(unit,mesh_Nelems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)  ! get affected elements
         do i = 1,contInts(1)
           e = mesh_FEasCP('elem',contInts(1+i))
           mesh_element(1+sv,e) = val
         enddo
         if (initialcondTableStyle == 0) read (unit,610,END=620) line      ! ignore IP range for old table style
         read (unit,610,END=630) line
         pos = IO_stringPos(line,1)
       enddo
     endif
   else   
     read (unit,610,END=630) line
   endif
 enddo

630 return

 endsubroutine



!********************************************************************
! store FEid, type, mat, tex, and node list per element
!
! allocate globals:
! _element
!********************************************************************
 subroutine mesh_abaqus_build_elements (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), parameter :: maxNchunks = 65
 integer(pInt), dimension (1+2*maxNchunks) :: pos

 integer(pInt) unit,i,j,count,e,t,homog,micro
 logical inPart,materialFound
 character (len=64) materialName,elemSetName
 character(len=300) line

 allocate (mesh_element (4+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(unit)
 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,2)
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,pos,1)) == '*end' .and. &
        IO_lc(IO_stringValue(line,pos,2)) == 'part' ) inPart = .false.

   if( inPart .and. &
       IO_lc(IO_stringValue(line,pos,1)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,pos,2)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,pos,2)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_stringValue(line,pos,2),'type'))  ! remember elem type
     count = IO_countDataLines(unit)
     do i = 1,count
       backspace(unit)
     enddo
     do i = 1,count
       read (unit,610,END=620) line
       pos = IO_stringPos(line,maxNchunks)                                  ! limit to 64 nodes max
       e = mesh_FEasCP('elem',IO_intValue(line,pos,1))
       if (e /= 0) then                                                     ! disregard non CP elems
         mesh_element(1,e) = IO_intValue(line,pos,1)                        ! FE id
         mesh_element(2,e) = t                                              ! elem type
         forall (j=1:FE_Nnodes(t)) &
           mesh_element(4+j,e) = IO_intValue(line,pos,1+j)                  ! copy FE ids of nodes to position 5:
         call IO_skipChunks(unit,FE_NoriginalNodes(t)-(pos(1)-1))           ! read on (even multiple lines) if FE_NoriginalNodes exceeds required node count
       endif
     enddo
   endif
 enddo

 
620 rewind(unit) ! just in case "*material" definitions apear before "*element"

 materialFound = .false.
 do 
   read (unit,610,END=630) line
   pos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_StringValue(line,pos,1)))
     case('*material')
       materialName = IO_extractValue(IO_lc(IO_StringValue(line,pos,2)),'name')    ! extract name=value
       materialFound = materialName /= ''                                   ! valid name?
     case('*user')
       if ( IO_lc(IO_StringValue(line,pos,2)) == 'material' .and. &
            materialFound ) then
         do i = 1,mesh_Nmaterials                                           ! look thru material names
           if (materialName == mesh_nameMaterial(i)) exit                   ! found one
         enddo
         if (i <= mesh_Nmaterials) then                                     ! found one?
           elemSetName = mesh_mapMaterial(i)                                ! take corresponding elemSet
           read (unit,610,END=630) line                                     ! read homogenization and microstructure
           pos = IO_stringPos(line,2)
           homog = NINT(IO_floatValue(line,pos,1))
           micro = NINT(IO_floatValue(line,pos,2))
           do i = 1,mesh_NelemSets                                          ! look thru all elemSet definitions
             if (elemSetName == mesh_nameElemSet(i)) then                   ! matched?
               do j = 1,mesh_mapElemSet(1,i)
                 e = mesh_FEasCP('elem',mesh_mapElemSet(1+j,i))
                 mesh_element(3,e) = homog                                  ! store homogenization
                 mesh_element(4,e) = micro                                  ! store microstructure
                 mesh_maxValStateVar(1) = max(mesh_maxValStateVar(1),homog)
                 mesh_maxValStateVar(2) = max(mesh_maxValStateVar(2),micro)
               enddo
             endif
           enddo
         endif
         materialFound = .false.
       endif
   endselect
 enddo

630 return

 endsubroutine


!********************************************************************
! get maximum count of shared elements among cpElements and
! build list of elements shared by each node in mesh
!
! _maxNsharedElems
! _sharedElem
!********************************************************************
 subroutine mesh_build_sharedElems ()

 use prec, only: pInt
 use IO
 implicit none

 integer(pint) e,t,n,j
 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt), dimension (:), allocatable :: node_seen

 allocate(node_seen(maxval(FE_Nnodes)))
 node_count = 0_pInt

 do e = 1,mesh_NcpElems

   t = mesh_element(2,e)                                                  ! my type
   node_seen = 0_pInt                                                     ! reset node duplicates
   do j = 1,FE_Nnodes(t)                                                  ! check each node of element
     n = mesh_FEasCP('node',mesh_element(4+j,e))
     if (all(node_seen /= n)) node_count(n) = node_count(n) + 1_pInt      ! if FE node not yet encountered -> count it
     node_seen(j) = n                                                     ! remember this node to be counted already
   enddo

 enddo

 mesh_maxNsharedElems = maxval(node_count)                                ! most shared node
 
 allocate ( mesh_sharedElem( 1+mesh_maxNsharedElems,mesh_Nnodes) )
 mesh_sharedElem = 0_pInt

 do e = 1,mesh_NcpElems
   node_seen = 0_pInt
   do j = 1,FE_Nnodes(mesh_element(2,e))
     n = mesh_FEasCP('node',mesh_element(4+j,e))
     if (all(node_seen /= n)) then
       mesh_sharedElem(1,n) = mesh_sharedElem(1,n) + 1                    ! count for each node the connected elements
       mesh_sharedElem(1+mesh_sharedElem(1,n),n) = e                      ! store the respective element id
     endif
     node_seen(j) = n
   enddo
 enddo

 deallocate (node_seen)

 return

 endsubroutine
 
 

!***********************************************************
! build up of IP neighborhood
!
! allocate globals
! _ipNeighborhood
!***********************************************************
 subroutine mesh_build_ipNeighborhood()
 
 use prec, only: pInt
 implicit none
 
 integer(pInt) e,t,i,j,k,n
 integer(pInt) neighbor,neighboringElem,neighboringIP,matchingElem
 integer(pInt), dimension(2) :: linkedNode = 0_pInt

 allocate(mesh_ipNeighborhood(2,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ; mesh_ipNeighborhood = 0_pInt
 
 do e = 1,mesh_NcpElems                                                   ! loop over cpElems
   t = mesh_element(2,e)                                                  ! get elemType
   do i = 1,FE_Nips(t)                                                    ! loop over IPs of elem
     do n = 1,FE_NipNeighbors(t)                                          ! loop over neighbors of IP
       neighbor = FE_ipNeighbor(n,i,t)
       if (neighbor > 0) then                                             ! intra-element IP
         neighboringElem = e
         neighboringIP   = neighbor
       else                                                               ! neighboring element's IP
         neighboringElem = 0_pInt
         neighboringIP = 0_pInt
         matchingElem = mesh_faceMatch(-neighbor,e)                       ! get CP elem id of face match
         if (matchingElem > 0 .and. mesh_element(2,matchingElem) == t) then  ! found match of same type?
           if (FE_nodesAtIP(2,1,i,t) == 0) then                           ! single linked node
 matchNode1: do j = 1,FE_Nnodes(t)                                        ! check against all neighbor's nodes
               if (mesh_element(4+FE_nodesAtIP(1,1,i,t),e)==mesh_element(4+j,matchingElem)) then
                 linkedNode(1) = j                                        ! which neighboring node matches my first nodeAtIP (indexed globally)
                 linkedNode(2) = 0_pInt
                 exit matchNode1
               endif
             enddo matchNode1
 matchFace1: do j = 1,FE_Nips(t)
               if ( (linkedNode(1) == FE_nodesAtIP(1,1,j,t)) .and. &
                    (FE_nodesAtIP(2,1,j,t) == 0) ) then
                 neighboringElem = matchingElem
                 neighboringIP = j
                 exit matchFace1
               endif
             enddo matchFace1
           else                                                           ! double linked node
 matchNode2: do j = 1,FE_Nnodes(t)                                        ! check against all neighbor's nodes
               if (mesh_element(4+FE_nodesAtIP(1,1,i,t),e)==mesh_element(4+j,matchingElem)) linkedNode(1) = j   ! which neighboring node matches my first nodeAtIP (indexed globally)
               if (mesh_element(4+FE_nodesAtIP(2,1,i,t),e)==mesh_element(4+j,matchingElem)) linkedNode(2) = j   ! which neighboring node matches my second nodeAtIP (indexed globally)
             enddo matchNode2
 matchFace2: do j = 1,FE_Nips(t)
               if ((linkedNode(1) == FE_nodesAtIP(1,1,j,t) .and. linkedNode(2) == FE_nodesAtIP(2,1,j,t)) .or. &
                   (linkedNode(1) == FE_nodesAtIP(2,1,j,t) .and. linkedNode(2) == FE_nodesAtIP(1,1,j,t)) .or. &
                   (linkedNode(1) == FE_nodesAtIP(1,2,j,t) .and. linkedNode(2) == FE_nodesAtIP(2,2,j,t)) .or. &
                   (linkedNode(1) == FE_nodesAtIP(2,2,j,t) .and. linkedNode(2) == FE_nodesAtIP(1,2,j,t))) then
                 neighboringElem = matchingElem
                 neighboringIP = j
                 exit matchFace2
               endif
             enddo matchFace2
           endif
         endif
       endif
       mesh_ipNeighborhood(1,n,i,e) = neighboringElem
       mesh_ipNeighborhood(2,n,i,e) = neighboringIP
     enddo
   enddo
 enddo
 
 return

 endsubroutine



!***********************************************************
! assignment of coordinates for subnodes in each cp element
!
! allocate globals
! _subNodeCoord
!***********************************************************
 subroutine mesh_build_subNodeCoords()
 
 use prec, only: pInt,pReal
 implicit none

 integer(pInt) e,t,n,p
 
 allocate(mesh_subNodeCoord(3,mesh_maxNnodes+mesh_maxNsubNodes,mesh_NcpElems)) ; mesh_subNodeCoord = 0.0_pReal
 
 do e = 1,mesh_NcpElems                   ! loop over cpElems
   t = mesh_element(2,e)                  ! get elemType
   do n = 1,FE_Nnodes(t)
     mesh_subNodeCoord(:,n,e) = mesh_node(:,mesh_FEasCP('node',mesh_element(4+n,e))) ! loop over nodes of this element type
   enddo
   do n = 1,FE_NsubNodes(t)               ! now for the true subnodes
     do p = 1,FE_Nips(t)                  ! loop through possible parent nodes
       if (FE_subNodeParent(p,n,t) > 0) & ! valid parent node
         mesh_subNodeCoord(:,n+FE_Nnodes(t),e) = &
         mesh_subNodeCoord(:,n+FE_Nnodes(t),e) + &
         mesh_node(:,mesh_FEasCP('node',mesh_element(4+FE_subNodeParent(p,n,t),e))) ! add up parents
     enddo
     mesh_subNodeCoord(:,n+FE_Nnodes(t),e) = mesh_subNodeCoord(:,n+FE_Nnodes(t),e) / count(FE_subNodeParent(:,n,t) > 0)
   enddo
 enddo 

 return
 
 endsubroutine


!***********************************************************
! calculation of IP volume
!
! allocate globals
! _ipVolume
!***********************************************************
 subroutine mesh_build_ipVolumes()
 
 use prec, only: pInt
 use math, only: math_volTetrahedron
 implicit none
 
 integer(pInt) e,f,t,i,j,k,n
 integer(pInt), parameter :: Ntriangles = FE_NipFaceNodes-2                     ! each interface is made up of this many triangles
 logical(pInt), dimension(mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNode      ! flagList to find subnodes determining center of grav
 real(pReal), dimension(3,mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNodePos   ! coordinates of subnodes determining center of grav
 real(pReal), dimension(3,FE_NipFaceNodes) :: nPos                              ! coordinates of nodes on IP face
 real(pReal), dimension(Ntriangles,FE_NipFaceNodes) :: volume                   ! volumes of possible tetrahedra
 real(pReal), dimension(3) :: centerOfGravity

 allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems)) ;                      mesh_ipVolume = 0.0_pReal
 allocate(mesh_ipCenterOfGravity(3,mesh_maxNips,mesh_NcpElems)) ;  mesh_ipCenterOfGravity = 0.0_pReal
 
 do e = 1,mesh_NcpElems                                    ! loop over cpElems
   t = mesh_element(2,e)                                   ! get elemType
   do i = 1,FE_Nips(t)                                     ! loop over IPs of elem
     gravityNode = .false.                                 ! reset flagList
     gravityNodePos = 0.0_pReal                            ! reset coordinates
     do f = 1,FE_NipNeighbors(t)                           ! loop over interfaces of IP
       do n = 1,FE_NipFaceNodes                            ! loop over nodes on interface
         gravityNode(FE_subNodeOnIPFace(n,f,i,t)) = .true.
         gravityNodePos(:,FE_subNodeOnIPFace(n,f,i,t)) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       enddo
     enddo
     
     do j = 1,mesh_maxNnodes+mesh_maxNsubNodes-1           ! walk through entire flagList except last
       if (gravityNode(j)) then                            ! valid node index
         do k = j+1,mesh_maxNnodes+mesh_maxNsubNodes       ! walk through remainder of list
           if (gravityNode(k) .and. all(abs(gravityNodePos(:,j) - gravityNodePos(:,k)) < 1.0e-100_pReal)) then   ! found duplicate
             gravityNode(j) = .false.                      ! delete first instance
             gravityNodePos(:,j) = 0.0_pReal
             exit                                          ! continue with next suspect
           endif
         enddo
       endif
     enddo
     centerOfGravity = sum(gravityNodePos,2)/count(gravityNode)
     
     do f = 1,FE_NipNeighbors(t)         ! loop over interfaces of IP and add tetrahedra which connect to CoG
       forall (n = 1:FE_NipFaceNodes) nPos(:,n) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       forall (n = 1:FE_NipFaceNodes, j = 1:Ntriangles) &  ! start at each interface node and build valid triangles to cover interface
         volume(j,n) = math_volTetrahedron(nPos(:,n), &    ! calc volume of respective tetrahedron to CoG
                                           nPos(:,1+mod(n+j-1,FE_NipFaceNodes)), &
                                           nPos(:,1+mod(n+j-0,FE_NipFaceNodes)), &
                                           centerOfGravity)
       mesh_ipVolume(i,e) = mesh_ipVolume(i,e) + sum(volume)    ! add contribution from this interface
     enddo
     mesh_ipVolume(i,e) = mesh_ipVolume(i,e) / FE_NipFaceNodes  ! renormalize with interfaceNodeNum due to loop over them
     mesh_ipCenterOfGravity(:,i,e) = centerOfGravity
   enddo
 enddo
 return
 
 endsubroutine


!***********************************************************
! calculation of IP interface areas
!
! allocate globals
! _ipArea, _ipAreaNormal
!***********************************************************
 subroutine mesh_build_ipAreas()
 
 use prec, only: pInt,pReal
 use math
 implicit none
 
 integer(pInt) e,f,t,i,j,n
 integer(pInt), parameter :: Ntriangles = FE_NipFaceNodes-2     ! each interface is made up of this many triangles
 real(pReal), dimension (3,FE_NipFaceNodes) :: nPos             ! coordinates of nodes on IP face
 real(pReal), dimension(3,Ntriangles,FE_NipFaceNodes) :: normal
 real(pReal), dimension(Ntriangles,FE_NipFaceNodes)   :: area

 allocate(mesh_ipArea(mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ;         mesh_ipArea       = 0.0_pReal
 allocate(mesh_ipAreaNormal(3,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ; mesh_ipAreaNormal = 0.0_pReal
 
 do e = 1,mesh_NcpElems                                         ! loop over cpElems
   t = mesh_element(2,e)                                        ! get elemType
   do i = 1,FE_Nips(t)                                          ! loop over IPs of elem
     do f = 1,FE_NipNeighbors(t)                                ! loop over interfaces of IP 
	     forall (n = 1:FE_NipFaceNodes) nPos(:,n) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       forall (n = 1:FE_NipFaceNodes, j = 1:Ntriangles)         ! start at each interface node and build valid triangles to cover interface
         normal(:,j,n) = math_vectorproduct(nPos(:,1+mod(n+j-1,FE_NipFaceNodes)) - nPos(:,n), &    ! calc their normal vectors
                                            nPos(:,1+mod(n+j-0,FE_NipFaceNodes)) - nPos(:,n))
         area(j,n) = dsqrt(sum(normal(:,j,n)*normal(:,j,n)))                                       ! and area
       end forall
	     forall (n = 1:FE_NipFaceNodes, j = 1:Ntriangles, area(j,n) > 0.0_pReal) &
	       normal(:,j,n) = normal(:,j,n) / area(j,n)            ! make unit normal
       
       mesh_ipArea(f,i,e) = sum(area) / (FE_NipFaceNodes*2.0_pReal)                   ! area of parallelograms instead of triangles
	     mesh_ipAreaNormal(:,f,i,e) = sum(sum(normal,3),2) / count(area > 0.0_pReal)  ! average of all valid normals
     enddo
   enddo
 enddo
 return
 
 endsubroutine


!***********************************************************
! write statistics regarding input file parsing
! to the output file
! 
!***********************************************************
 subroutine mesh_tell_statistics()

 use prec, only: pInt
 use math, only: math_range
 use IO, only: IO_error

 implicit none

 integer(pInt), dimension (:,:), allocatable :: mesh_HomogMicro
 character(len=64) fmt

 integer(pInt) i,e,n,f,t

 if (mesh_maxValStateVar(1) == 0) call IO_error(110) ! no materials specified
 if (mesh_maxValStateVar(2) == 0) call IO_error(120) ! no textures specified
   
 allocate (mesh_HomogMicro(mesh_maxValStateVar(1),mesh_maxValStateVar(2))); mesh_HomogMicro = 0_pInt
 do i=1,mesh_NcpElems
   mesh_HomogMicro(mesh_element(3,i),mesh_element(4,i)) = &
   mesh_HomogMicro(mesh_element(3,i),mesh_element(4,i)) + 1 ! count combinations of homogenization and microstructure
 enddo
 
!$OMP CRITICAL (write2out)

! write(6,*)
! write(6,*) 'Input Parser: IP COORDINATES'
! write(6,'(a5,x,a5,3(x,a12))') 'elem','IP','x','y','z'
! do e = 1,mesh_NcpElems
!   do i = 1,FE_Nips(mesh_element(2,e))
!     write (6,'(i5,x,i5,3(x,f12.8))') e, i, mesh_ipCenterOfGravity(:,i,e)
!   enddo
! enddo 
! write(6,*)
! write(6,*) "Input Parser: IP NEIGHBORHOOD"
! write(6,*)
! write(6,"(a10,x,a10,x,a10,x,a3,x,a13,x,a13)") "elem","IP","neighbor","","elemNeighbor","ipNeighbor"
! do e = 1,mesh_NcpElems                  ! loop over cpElems
!   t = mesh_element(2,e)                 ! get elemType
!   do i = 1,FE_Nips(t)                   ! loop over IPs of elem
!     do n = 1,FE_NipNeighbors(t)         ! loop over neighbors of IP
!       write (6,"(i10,x,i10,x,i10,x,a3,x,i13,x,i13)") e,i,n,'-->',mesh_ipNeighborhood(1,n,i,e),mesh_ipNeighborhood(2,n,i,e)
!     enddo
!   enddo
! enddo
! write (6,*)
!  write (6,*) "Input Parser: ELEMENT VOLUME"
!  write (6,*)
!  write (6,"(a13,x,e15.8)") "total volume", sum(mesh_ipVolume)
!  write (6,*)
!  write (6,"(a5,x,a5,x,a15,x,a5,x,a15,x,a16)") "elem","IP","volume","face","area","-- normal --"
!  do e = 1,mesh_NcpElems
!    do i = 1,FE_Nips(mesh_element(2,e))
!      write (6,"(i5,x,i5,x,e15.8)") e,i,mesh_IPvolume(i,e)
!      do f = 1,FE_NipNeighbors(mesh_element(2,e))
!        write (6,"(i33,x,e15.8,x,3(f6.3,x))") f,mesh_ipArea(f,i,e),mesh_ipAreaNormal(:,f,i,e)
!      enddo
!    enddo
!  enddo
!  write (6,*)
! write (6,*) "Input Parser: SUBNODE COORDINATES"
! write (6,*)
! write(6,'(a5,x,a5,x,a15,x,a15,x,a20,3(x,a12))') 'elem','IP','IP neighbor','IPFaceNodes','subNodeOnIPFace','x','y','z'
! do e = 1,mesh_NcpElems                  ! loop over cpElems
!   t = mesh_element(2,e)                 ! get elemType
!   do i = 1,FE_Nips(t)                   ! loop over IPs of elem
!     do f = 1,FE_NipNeighbors(t)         ! loop over interfaces of IP
!       do n = 1,FE_NipFaceNodes          ! loop over nodes on interface
!         write(6,'(i5,x,i5,x,i15,x,i15,x,i20,3(x,f12.8))') e,i,f,n,FE_subNodeOnIPFace(n,f,i,t),&
!                                                            mesh_subNodeCoord(1,FE_subNodeOnIPFace(n,f,i,t),e),&
!                                                            mesh_subNodeCoord(2,FE_subNodeOnIPFace(n,f,i,t),e),&
!                                                            mesh_subNodeCoord(3,FE_subNodeOnIPFace(n,f,i,t),e)
!       enddo
!     enddo
!   enddo
! enddo
 write (6,*)
 write (6,*) "Input Parser: STATISTICS"
 write (6,*)
 write (6,*) mesh_Nelems,           " : total number of elements in mesh"
 write (6,*) mesh_NcpElems,         " : total number of CP elements in mesh"
 write (6,*) mesh_Nnodes,           " : total number of nodes in mesh"
 write (6,*) mesh_maxNnodes,        " : max number of nodes in any CP element"
 write (6,*) mesh_maxNips,          " : max number of IPs in any CP element"
 write (6,*) mesh_maxNipNeighbors,  " : max number of IP neighbors in any CP element"
 write (6,*) mesh_maxNsubNodes,     " : max number of (additional) subnodes in any CP element"
 write (6,*) mesh_maxNsharedElems,  " : max number of CP elements sharing a node"
 write (6,*)
 write (6,*) "Input Parser: HOMOGENIZATION/MICROSTRUCTURE"
 write (6,*)
 write (6,*) mesh_maxValStateVar(1), " : maximum homogenization index"
 write (6,*) mesh_maxValStateVar(2), " : maximum microstructure index"
 write (6,*)
 write (fmt,"(a,i5,a)") "(9(x),a2,x,",mesh_maxValStateVar(2),"(i8))"
 write (6,fmt) "+-",math_range(mesh_maxValStateVar(2))
 write (fmt,"(a,i5,a)") "(i8,x,a2,x,",mesh_maxValStateVar(2),"(i8))"
 do i=1,mesh_maxValStateVar(1)      ! loop over all (possibly assigned) homogenizations
   write (6,fmt) i,"| ",mesh_HomogMicro(i,:) ! loop over all (possibly assigned) microstrcutures
 enddo
 write (6,*)
 call flush(6)
!$OMP END CRITICAL (write2out)

 return
 
 endsubroutine

 
 END MODULE mesh
 
