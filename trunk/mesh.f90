
!##############################################################
 MODULE mesh     
!##############################################################

 use prec, only: pReal,pInt
 implicit none

! ---------------------------
! _Nelems    : total number of elements in mesh
! _NcpElems  : total number of CP elements in mesh
! _NelemTypes: total number of element types in mesh
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
 integer(pInt) mesh_Nelems,mesh_NcpElems,mesh_NelemTypes
 integer(pInt) mesh_Nnodes,mesh_maxNnodes,mesh_maxNips,mesh_maxNsharedElems
 integer(pInt), dimension(:,:), allocatable, target :: mesh_mapFEtoCPelem,mesh_mapFEtoCPnode
 integer(pInt), dimension(:,:), allocatable :: mesh_element, mesh_sharedElem
 integer(pInt), dimension(:,:,:,:), allocatable :: mesh_ipNeighborhood
 real(pReal), allocatable :: mesh_node (:,:)

 integer(pInt), parameter :: FE_Nelemtypes = 2
 integer(pInt), parameter :: FE_maxNnodes = 8
 integer(pInt), parameter :: FE_maxNips = 8
 integer(pInt), parameter :: FE_maxNneighbors = 6
 integer(pInt), parameter :: FE_maxNfaceNodes = 4
 integer(pInt), parameter :: FE_maxNfaces = 6
 integer(pInt), dimension(200):: FE_mapElemtype 
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_Nnodes = &
 (/8, & ! element 7
   4  & ! element 134
  /)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_Nips = &
 (/8, & ! element 7
   1  & ! element 134
  /)
 integer(pInt), dimension(FE_Nelemtypes), parameter :: FE_NipNeighbors = &
 (/6, & ! element 7
   4  & ! element 134
  /)
 integer(pInt), dimension(FE_maxNfaces,FE_Nelemtypes), parameter :: FE_NfaceNodes = &
 reshape((/&
  4,4,4,4,4,4, & ! element 7
  3,3,3,3,0,0  & ! element 134
   /),(/FE_maxNfaces,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNips,FE_Nelemtypes), parameter :: FE_nodeAtIP = &
 reshape((/&
  1,2,4,3,5,6,8,7, & ! element 7
  1,0,0,0,0,0,0,0  & ! element 134
   /),(/FE_maxNips,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNnodes,FE_Nelemtypes), parameter :: FE_ipAtNode = &
 reshape((/&
  1,2,4,3,5,6,8,7, & ! element 7
  1,1,1,1,0,0,0,0  & ! element 134
   /),(/FE_maxNnodes,FE_Nelemtypes/))
 integer(pInt), dimension(FE_maxNfaceNodes,FE_maxNfaces,FE_Nelemtypes), parameter :: FE_nodeOnFace = &
 reshape((/&
  1,2,3,4 , & ! element 7
  2,1,5,6 , &
  3,2,6,7 , &
  3,4,8,7 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,0 , & ! element 134
  1,4,2,0 , &
  2,3,4,0 , &
  1,3,4,0 , &
  0,0,0,0 , &
  0,0,0,0 &
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
  -3, 7,-4, 6,-6, 4 , &
  -1,-2,-3,-4, 0, 0 , & ! element 134
   0, 0, 0, 0, 0, 0 , &
   0, 0, 0, 0, 0, 0 , &
   0, 0, 0, 0, 0, 0 , &
   0, 0, 0, 0, 0, 0 , &
   0, 0, 0, 0, 0, 0 , &
   0, 0, 0, 0, 0, 0 , &
   0, 0, 0, 0, 0, 0   &
   /),(/FE_maxNneighbors,FE_maxNips,FE_Nelemtypes/))
   
 CONTAINS
! ---------------------------
! subroutine mesh_init()
! function mesh_FEtoCPelement(FEid)
! function mesh_build_ipNeighorhood()
! subroutine mesh_parse_inputFile()
! ---------------------------


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

 FE_mapElemtype = 1 ! MISSING this should be zero...
 FE_mapElemtype(  7) = 1
 FE_mapElemtype(134) = 2

! call to various subrountes to parse the stuff from the input file...
 END SUBROUTINE
 
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
 do while (upper-lower > 1)
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
 

 
!********************************************************************
! Build node mapping from FEM to CP
!********************************************************************
 SUBROUTINE mesh_build_nodeMapping (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt) unit,i,j,Nnodes,cur_node
 integer(pInt), dimension (133) :: pos
 character*264 line

610 FORMAT(A264)

 rewind(unit)
 allocate ( mesh_mapFEtoCPnode(2,mesh_Nnodes) )
 mesh_mapFEtoCPnode(:,:) = 0_pInt
 node_count(:) = 0_pInt

 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'coordinates' ) then
     read (unit,610,END=620) line  ! skip crap line
     do i=1,mesh_Nnodes
       read (unit,610,END=620) line
       mesh_mapFEtoCPnode(1,i) = IO_fixedIntValue (line,(/0,10/),1)
       mesh_mapFEtoCPnode(2,i) = i
     end do
   end if
 end do
620 continue

 do i=2,mesh_Nnodes
   if( mesh_mapFEtoCPnode(1,i).lt.mesh_mapFEtoCPnode(1,i-1) )then
     write(*,*) 'Need to sort node'
   end if
 end do

 rewind(unit)
 do
   read (unit,610,END=630) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'connectivity' ) then
     read (unit,610,END=620) line  ! Garbage line
	 do i=1,mesh_Nelems
	   read (unit,610,END=620) line
	   pos = IO_stringPos(line,66)  ! limit to 64 nodes max (plus ID, type)
	   Nnodes = FE_Nnodes(FE_mapElemtype(IO_intValue(line,pos,2)))
	   do j=1,Nnodes
	     cur_node = IO_IntValue (line,pos,j+2)
	     node_count( mesh_FEasCP('node',cur_node) )= node_count( mesh_FEasCP('node',cur_node) )+1
	   end do
     end do
   end if
 end do
630 continue

 mesh_maxNsharedElems = MAXVAL(node_count)

 return
 END SUBROUTINE

!********************************************************************
! Build ele mapping from FEM to CP
!********************************************************************
 SUBROUTINE mesh_build_CPeleMapping (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer unit, i,cur_CPele,start_ele,end_ele
 character*264 line
 integer(pInt), dimension (41) :: pos


610 FORMAT(A264)
 rewind(unit)

 allocate ( mesh_mapFEtoCPelem(2,mesh_NcpElems) )
 mesh_mapFEtoCPelem(:,:) = 0_pInt
 cur_CPele = 0_pInt

 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'hypoelastic' ) then
     do i=1,4
       read (unit,610,END=620) line
     end do
     pos = IO_stringPos(line,20)
	 if( IO_lc(IO_Stringvalue(line,pos,2)).eq.'to' )then
       start_ele = IO_IntValue(line,pos,1)
	   end_ele = IO_IntValue(line,pos,3)
	   do i=start_ele,end_ele
	     cur_CPele = cur_CPele+1
	     mesh_mapFEtoCPelem(1,cur_CPele) = i
		 mesh_mapFEtoCPelem(2,cur_CPele) = cur_CPele
       end do
     else
       do while( IO_lc(IO_Stringvalue(line,pos,pos(1))).eq.'c' )
	     do i=1,pos(1)-1
	   	   cur_CPele = cur_CPele+1
           mesh_mapFEtoCPelem(1,cur_CPele) = IO_IntValue(line,pos,i)
		   mesh_mapFEtoCPelem(2,cur_CPele) = cur_CPele
	     end do
		 read (unit,610,END=620) line
		 pos = IO_stringPos(line,20)
       end do
	   do i=1,pos(1)
	   	 cur_CPele = cur_CPele+1
         mesh_mapFEtoCPelem(1,cur_CPele) = IO_IntValue(line,pos,i)
		 mesh_mapFEtoCPelem(2,cur_CPele) = cur_CPele
	   end do
     end if
   end if
 end do
620 continue

 do i=2,mesh_NcpElems
   if( mesh_mapFEtoCPelem(1,i).lt.mesh_mapFEtoCPelem(1,i-1) )then
     write(*,*) 'Need to sort ele'
   end if
 end do

 return
 END SUBROUTINE


!********************************************************************
!********************************************************************
 SUBROUTINE mesh_build_Sharedelems (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer unit
 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt), dimension (41) :: pos
 integer i,j,FE_node,CP_node,Nnodes,CP_elem
 character*264 line

610 FORMAT(A264)
 rewind(unit)

 allocate ( mesh_sharedElem( 1+mesh_maxNsharedElems,mesh_Nnodes) )
 mesh_sharedElem(:,:) = 0_pInt
 node_count(:) = 0_pInt

 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'connectivity' ) then
     read (unit,610,END=620) line  ! Garbage line
	 do i=1,mesh_Nelems
	   read (unit,610,END=620) line
	   pos = IO_stringPos(line,66)  ! limit to 64 nodes max (plus ID, type)
	   Nnodes = FE_Nnodes(FE_mapElemtype(IO_intValue(line,pos,2)))
	   CP_elem = mesh_FEasCP('elem',IO_IntValue(line,pos,1))
	   if( CP_elem.ne.0 )then
	     do j=1,Nnodes
	       FE_node = IO_IntValue (line,pos,j+2)
           CP_node = mesh_FEasCP('node',FE_node)
	       node_count( CP_node )= node_count( CP_node )+1
		   mesh_sharedElem(node_count(CP_node)+1,CP_node) = CP_elem
	     end do
	   end if
     end do
   end if
 end do
620 continue

 do i=1,mesh_Nnodes
   mesh_sharedElem(1,i) = node_count(i)
 end do

 return
 END SUBROUTINE

!********************************************************************
!********************************************************************
 SUBROUTINE mesh_build_nodeCoord (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer unit,i,j,m
 integer(pInt), dimension(3) :: pos
 integer(pInt), dimension(5), parameter :: node_ends = (/0,10,30,50,70/)
 character*264 line

 rewind(unit)
 allocate ( mesh_node (3,mesh_Nnodes) )
 mesh_node(:,:) = 0_pInt

610 FORMAT(A264)

 do while(.true.)
   read (unit,610,END=620) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'coordinates' ) then
     read (unit,610,END=620) line ! skip crap line
     do i=1,mesh_Nnodes
       read (unit,610,END=620) line
	   m = mesh_FEasCP('node',IO_fixedIntValue (line,node_ends,1))
       do j=1,3
         mesh_node(j,m) = IO_fixedNoEFloatValue (line,node_ends,j+1)
       end do
     end do
   end if

 end do
620 continue

 return
 END SUBROUTINE

!********************************************************************
!********************************************************************
 SUBROUTINE mesh_build_element (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer unit
 integer FE_node,Nnodes,i,j,sv,ele,val,start_ele,end_ele
 integer(pInt), dimension (41) :: pos
 logical not_found
 character*264 line,line2

 rewind(unit)
 allocate ( mesh_element (mesh_Nelems,4+mesh_maxNnodes) )
 mesh_element(:,:) = 0_pInt

610 FORMAT(A264)

 do
   read (unit,610,END=620) line
   pos = IO_stringPos(line,2)

   if( IO_lc(IO_stringValue(line,pos,1)) == 'connectivity' ) then
     read (unit,610,END=620) line  ! Garbage line
	 do i=1,mesh_Nelems
	   read (unit,610,END=620) line
	   pos = IO_stringPos(line,66)  ! limit to 64 nodes max (plus ID, type)
	   Nnodes = FE_Nnodes(FE_mapElemtype(IO_intValue(line,pos,2)))
	   mesh_element (i,1) = IO_IntValue (line,pos,1)
	   mesh_element (i,2) = IO_IntValue (line,pos,2)
       do j=1,Nnodes
         FE_node = IO_IntValue (line,pos,j+2)
!         CP_node = mesh_FEasCP('node',FE_node)
         mesh_element(i,j+4) = FE_node
	   end do
     end do
   end if

   if( (IO_lc(IO_stringValue(line,pos,1)) == 'initial').and.    &
       (IO_lc(IO_stringValue(line,pos,2)) == 'state') ) then
     read (unit,610,END=620) line
	 pos = IO_stringPos(line,66)
	 sv = IO_IntValue (line,pos,1)
	 if( (sv.ne.2).and.(sv.ne.3) )then
	   write(*,*) 'Major PROBLEM!!  -> Invalid state variable found'
	   write(*,*) sv
	 end if

     read (unit,610,END=620) line
	 val = NINT(IO_fixedNoEFloatValue (line,(/0,20/),1))

	 line2 = line

     do while( line.eq.line2 )

       line2 = line
       read (unit,610,END=620) line
	   pos = IO_stringPos(line,20)
	   if ( IO_lc(IO_Stringvalue(line,pos,2)).eq.'to' )then
         start_ele = IO_IntValue(line,pos,1)
	     end_ele = IO_IntValue(line,pos,3)
	     do i=start_ele,end_ele
	       not_found = .true.
	       j=1
	       do while( not_found )
	         if( mesh_element(j,1).eq.i )then
		       not_found = .false.
		       ele = i
		     end if
		     j=j+1
	       end do
           mesh_element(ele,sv+1) = val 
         end do        
	   else

         do while( IO_lc(IO_Stringvalue(line,pos,pos(1))).eq.'c' )
	       do i=1,pos(1)-1
	         ele = IO_IntValue(line,pos,i)
	         not_found = .true.
	         j=1
	         do while( not_found )
	           if( mesh_element(j,1).eq.ele )then
		         not_found = .false.
		         ele = j
		       end if
		       j=j+1
	         end do
             mesh_element(ele,sv+1) = val 
	       end do
           read (unit,610,END=620) line
	       pos = IO_stringPos(line,20)
         end do
         do i=1,pos(1)
	       ele = IO_IntValue(line,pos,i)
	       not_found = .true.
	       j=1
	       do while( not_found )
	         if( mesh_element(j,1).eq.ele )then
		       not_found = .false.
		       ele = j
		     end if
		     j=j+1
	       end do
           mesh_element(ele,sv+1) = val 
	     end do

       end if
       read (unit,610,END=620) line  ! Garbage line
       read (unit,610,END=620) line
     end do
 
   end if
   pos = IO_stringPos(line,20)
   if( (IO_lc(IO_stringValue(line,pos,1)) == 'initial').and.    &
       (IO_lc(IO_stringValue(line,pos,2)) == 'state') ) then
	  backspace(unit)
   end if

 end do

620 continue

 return
 END SUBROUTINE

 
!********************************************************************
! Get global variables (#ele, #nodes, #ele types, #CP ele, max # nodes per element)
!********************************************************************
 SUBROUTINE mesh_get_globals (unit)

 use prec, only: pInt
 use IO
 implicit none

 integer(pInt) unit,i,pos(41),Nnodes
 character*264 line

610 FORMAT(A264)

 rewind(unit)
 mesh_NelemTypes = 0_pInt

 do 
   read (unit,610,END=620) line
   pos = IO_stringPos(line,20)

   select case ( IO_lc(IO_Stringvalue(line,pos,1)))
     case('sizing')
       mesh_Nelems = IO_IntValue (line,pos,3)
       mesh_Nnodes = IO_IntValue (line,pos,4)
     case('elements')
       mesh_NelemTypes = mesh_NelemTypes+1
	   Nnodes = FE_Nnodes(FE_mapElemtype(IO_intValue(line,pos,2)))
	   if( Nnodes.gt.mesh_maxNnodes )then
	     mesh_maxNnodes = Nnodes
	   end if
     case('hypoelastic')
       do i=1,4
         read (unit,610,END=620) line
       end do
       pos = IO_stringPos(line,20)
       if( IO_lc(IO_Stringvalue(line,pos,2)).eq.'to' )then
         mesh_NcpElems = IO_IntValue(line,pos,3)-IO_IntValue(line,pos,1)+1
       else
	     mesh_NcpElems = mesh_NcpElems + pos(1)
		 do while( IO_lc(IO_Stringvalue(line,pos,pos(1))).eq.'c' )
	       mesh_NcpElems = mesh_NcpElems - 1 ! Counted the c character from the line
		   read (unit,610,END=620) line
		   pos = IO_stringPos(line,20)
           mesh_NcpElems = mesh_NcpElems + pos(1)
		 end do
       end if

   end select

 end do
620 continue


 return
 END SUBROUTINE

 
 END MODULE mesh
 