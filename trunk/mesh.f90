
!##############################################################
 MODULE mesh     
!##############################################################

 use prec, only: pReal,pInt
 implicit none

! ---------------------------
! _Nelems    : total number of elements in mesh
! _NcpElems  : total number of CP elements in mesh
! _Nnodes    : total number of nodes in mesh
! _maxNnodes : max number of nodes in any element
! _maxNips   : max number of IPs in any element
! _maxNsharedElems : max number of elements sharing a node
!
! _element    : FEid, type, material, texture, node indices
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
! _mapElementtype   : map MARC/ABAQUS elemtype to 1-maxN 
!
! _Nnodes            : # nodes in a specific type of element
! _Nips              : # IPs in a specific type of element
! _NipNeighbors      : # IP neighbors in a specific type of element
! _ipNeighbor        : +x,-x,+y,-y,+z,-z list of intra-element IPs and
!     (negative) neighbor faces per own IP in a specific type of element
! _NfaceNodes        : # nodes per face in a specific type of element
! _nodeOnFace        : list of node indices on each face of a specific type of element
! _ipAtNode          : map node index to IP index in a specific type of element
! _nodeAtIP          : map IP index to node index in a specific type of element
! _ipNeighborhood    : 6 or less neighboring IPs as [element_num, IP_index]
!     order is +x,-x,+y,-y,+z,-z but meaning strongly depends on Elemtype
! ---------------------------
 integer(pInt) mesh_Nelems,mesh_NcpElems,mesh_Nnodes,mesh_maxNnodes,mesh_maxNips,mesh_maxNsharedElems
 integer(pInt), dimension(:,:), allocatable, target :: mesh_mapFEtoCPelem,mesh_mapFEtoCPnode
 integer(pInt), dimension(:,:), allocatable :: mesh_element, mesh_sharedElem
 integer(pInt), dimension(:,:,:,:), allocatable :: mesh_ipNeighborhood
 real(pReal), allocatable :: mesh_node (:,:)

 integer(pInt), parameter :: FE_Nelemtypes = 1
 integer(pInt), parameter :: FE_maxNnodes = 8
 integer(pInt), parameter :: FE_maxNips = 8
 integer(pInt), parameter :: FE_maxNneighbors = 6
 integer(pInt), parameter :: FE_maxNfaceNodes = 4
 integer(pInt), parameter :: FE_maxNfaces = 6
 integer(pInt), dimension(200) :: FE_mapElemtype 
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_Nnodes = &
 (/8/)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_Nips = &
 (/8/)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_NipNeighbors = &
 (/ 6 /)
 integer(pInt), dimension(FE_maxNfaces,FE_Nelemtypes), parameter :: FE_NfaceNodes = &
 reshape((/&
  4,4,4,4,4,4 & ! element 7
   /),(/FE_maxNfaces,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNips,FE_Nelemtypes), parameter :: FE_nodeAtIP = &
 reshape((/&
  1,2,4,3,5,6,8,7 & ! element 7
   /),(/FE_maxNips,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNnodes,FE_Nelemtypes), parameter :: FE_ipAtNode = &
 reshape((/&
  1,2,4,3,5,6,8,7 & ! element 7
   /),(/FE_maxNnodes,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNfaceNodes,FE_maxNfaces,FE_Nelemtypes), parameter :: FE_nodeOnFace = &
 reshape((/&
  1,2,3,4 , & ! element 7
  2,1,5,6 , &
  3,2,6,7 , &
  3,4,8,7 , &
  4,1,5,8 , &
  8,7,6,5 &
   /),(/FE_maxNfaceNodes,FE_maxNfaces,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNneighbors,FE_maxNips,FE_Nelemtypes), parameter :: FE_ipNeighbor = &
 reshape((/&
   2,-5, 3,-2, 5,-1 , & ! element 7
  -3, 1, 4,-2, 6,-1 , &
   4,-5,-4, 1, 7,-1 , &
  -3, 3,-4, 2, 8,-1 , &
   6,-5, 7,-2,-6, 1 , &
  -3, 5, 8,-2,-6, 2 , &
   8,-5,-4, 5,-6, 3 , &
  -3, 7,-4, 6,-6, 4   &
   /),(/FE_maxNneighbors,FE_maxNips,FE_Nelemtypes/))
   
 CONTAINS
! ---------------------------
! subroutine mesh_init()
! function mesh_FEtoCPelement(FEid)
! function mesh_build_ipNeighorhood()
! subroutine mesh_parse_inputFile()
! ---------------------------


 
!***********************************************************
! find face-matching element of same type
!
! 
!***********************************************************
 FUNCTION mesh_faceMatch(face,elem)

 use prec, only: pInt
 implicit none

 integer(pInt) face,elem
 integer(pInt) mesh_faceMatch
 integer(pInt), dimension(FE_NfaceNodes(face,mesh_element(2,elem))) :: nodeMapFE
 integer(pInt) minN,NsharedElems,lonelyNode,faceNode,i,j,t
 
 mesh_faceMatch = 0_pInt   ! intialize to "no match found"
 t = mesh_element(2,elem)  ! figure elemType
 do faceNode=1,FE_NfaceNodes(face,t)  ! loop over nodes on face
   nodeMapFE(faceNode) = mesh_element(4+FE_nodeOnFace(faceNode,face,t),elem) ! FE id of face node
   NsharedElems = mesh_sharedElem(1,nodeMapFE(faceNode))  ! figure # shared elements for this node
   if (NsharedElems < minN) then
     minN = NsharedElems       ! remember min # shared elems
     lonelyNode = faceNode     ! remember most lonely node
   endif
 end do
candidate: do i=1,minN  ! iterate over lonelyNode's shared elements
   mesh_faceMatch = mesh_sharedElem(1+i,nodeMapFE(lonelyNode)) ! present candidate elem
   if (mesh_faceMatch == elem) then   ! my own element ?
     mesh_faceMatch = 0_pInt          ! disregard
     cycle candidate
   endif
   do faceNode=1,FE_NfaceNodes(face,t) ! check remaining face nodes to match
     if (faceNode == lonelyNode) cycle ! disregard lonely node (matches anyway)
     NsharedElems = mesh_sharedElem(1,nodeMapFE(faceNode))  ! how many shared elems for checked node?
     do j=1,NsharedElems               ! iterate over other node's elements
       if (all(mesh_sharedElem(2:1+NsharedElems,nodeMapFE(faceNode)) /= mesh_faceMatch)) then  ! no ref to candidate elem?
         mesh_faceMatch = 0_pInt       ! set to "no match" (so far)
         cycle candidate               ! next candidate elem
       endif
     end do
   end do
 end do candidate
  
 return
 
 END FUNCTION


!***********************************************************
! build up of IP neighborhood
!***********************************************************
 SUBROUTINE mesh_build_ipNeighborhood()
 
 use prec, only: pInt
 implicit none
 
 integer(pInt) e,t,i,j,k,n
 integer(pInt) neighbor,neighboringElem,neighboringIP,matchingElem,faceNode,linkingNode
 
 do e = 1,mesh_NcpElems            ! loop over cpElems
   t = mesh_element(2,e)           ! get elemType
   do i = 1,FE_Nips(t)             ! loop over IPs of elem
     do n = 1,FE_NipNeighbors(t)   ! loop over neighbors of IP
       neighbor = FE_ipNeighbor(n,i,t)
       if (neighbor > 0) then ! intra-element IP
           neighboringElem = e
           neighboringIP   = neighbor
       else                   ! neighboring element's IP
         neighboringElem = 0_pInt
         neighboringIP = 0_pInt
         matchingElem = mesh_faceMatch(-neighbor,e)
         if (matchingElem > 0) then
matchFace: do j = 1,FE_NfaceNodes(-neighbor,t)  ! count over nodes on matching face
             faceNode = FE_nodeOnFace(j,-neighbor,t)  ! get face node id
             if (i == FE_ipAtNode(faceNode,t)) then   ! ip linked to face node is me?
               linkingNode = mesh_element(4+faceNode,e) ! FE id of this facial node
               do k = 1,FE_Nnodes(t)   ! loop over nodes in matching element
                 if (linkingNode == mesh_element(4+k,matchingElem)) then
                   neighboringElem = matchingElem
                   neighboringIP = FE_ipAtNode(j,t)
                   exit matchFace
                 endif
               end do
             endif
           end do matchFace
         endif
       endif
       mesh_ipNeighborhood(1,n,i,e) = neighboringElem
       mesh_ipNeighborhood(2,n,i,e) = neighboringIP
     end do
   end do
 end do
 
 return
 
 END SUBROUTINE
 
 
!***********************************************************
! FE to CP id mapping by binary search thru lookup array
!
! valid questions are 'elem', 'node'
!***********************************************************
 FUNCTION mesh_FEasCP(what,id)
 
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
 end select
 
 lower = 1_pInt
 upper = size(lookupMap,2)
 
 ! check at bounds
 if (lookupMap(1,lower) == id) then
   mesh_FEasCP = lookupMap(2,lower)
   return
 elseif (lookupMap(1,upper) == id) then
   mesh_FEasCP = lookupMap(2,upper)
   return
 endif
 
 ! binary search in between bounds
 do while (upper-lower > 0)
   center = (lower+upper)/2
   if (lookupMap(1,center) < id) then
     lower = center
   elseif (lookupMap(1,center) > id) then
     upper = center
   else
     mesh_FEasCP = lookupMap(2,center)
     exit
   end if
 end do
 return
 
 END FUNCTION
 

!***********************************************************
! initialization 
!***********************************************************
 SUBROUTINE mesh_init ()

 mesh_Nelems = 0_pInt
 mesh_NcpElems = 0_pInt
 mesh_Nnodes = 0_pInt
 mesh_maxNips = 0_pInt
 mesh_maxNnodes = 0_pInt
 mesh_maxNsharedElems = 0_pInt
 call mesh_parse_inputFile ()

 END SUBROUTINE
 
 
!***********************************************************
! parsing of input file 
!***********************************************************
 SUBROUTINE mesh_parse_inputFile()

 
 END SUBROUTINE
 
 END MODULE mesh
 