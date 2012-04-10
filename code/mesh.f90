! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
!* $Id$
!##############################################################
 MODULE mesh     
!##############################################################

 use prec, only: pReal, pInt
 implicit none
 private
 
 integer(pInt), public :: &
   mesh_NcpElems, &                                                                                 ! total number of CP elements in mesh
   mesh_NelemSets, &
   mesh_maxNelemInSet, &
   mesh_Nmaterials, &
   mesh_Nnodes, &                                                                                   ! total number of nodes in mesh
   mesh_maxNnodes, &                                                                                ! max number of nodes in any CP element
   mesh_maxNips, &                                                                                  ! max number of IPs in any CP element
   mesh_maxNipNeighbors, &                                                                          ! max number of IP neighbors in any CP element
   mesh_maxNsharedElems, &                                                                          ! max number of CP elements sharing a node
   mesh_maxNsubNodes

 integer(pInt), dimension(:,:), allocatable, public :: &
   mesh_element, &                                                                                  ! FEid, type(internal representation), material, texture, node indices
   mesh_sharedElem, &                                                                               ! entryCount and list of elements containing node
   mesh_nodeTwins                                                                                   ! node twins are surface nodes that lie exactly on opposite sides of the mesh (surfaces nodes with equal coordinate values in two dimensions)
 
 integer(pInt), dimension(:,:,:,:), allocatable, public :: &
   mesh_ipNeighborhood                                                                              ! 6 or less neighboring IPs as [element_num, IP_index]

 real(pReal), dimension(:,:), allocatable, public :: &
   mesh_ipVolume, &                                                                                 ! volume associated with IP (initially!)
   mesh_node0, &                                                                                    ! node x,y,z coordinates (initially!)
   mesh_node                                                                                        ! node x,y,z coordinates (after deformation! ONLY FOR MARC!!!)
 
 real(pReal), dimension(:,:,:), allocatable, public :: &
   mesh_ipCenterOfGravity, &                                                                        ! center of gravity of IP (after deformation!)
   mesh_ipArea                                                                                      ! area of interface to neighboring IP (initially!)
 
 real(pReal),dimension(:,:,:,:), allocatable, public :: & 
   mesh_ipAreaNormal                                                                                ! area normal of interface to neighboring IP (initially!)
    
 logical, dimension(3), public :: mesh_periodicSurface                                              ! flag indicating periodic outer surfaces (used for fluxes)
                                                              
 integer(pInt), private :: &
   mesh_Nelems, &                                                                                   ! total number of elements in mesh
   hypoelasticTableStyle, &
   initialcondTableStyle
    
 integer(pInt), dimension(2), private :: &
   mesh_maxValStateVar = 0_pInt
             
 character(len=64), dimension(:), allocatable, private :: &
   mesh_nameElemSet, &                                                                              ! names of elementSet
   mesh_nameMaterial, &                                                                             ! names of material in solid section
   mesh_mapMaterial                                                                                 ! name of elementSet for material
     
 integer(pInt), dimension(:,:), allocatable, private :: &
   mesh_mapElemSet                                                                                  ! list of elements in elementSet
     
 integer(pInt), dimension(:,:), allocatable, target, private :: &
   mesh_mapFEtoCPelem, &                                                                            ! [sorted FEid, corresponding CPid]
   mesh_mapFEtoCPnode                                                                               ! [sorted FEid, corresponding CPid]
   
 real(pReal),dimension(:,:,:), allocatable, private :: &
   mesh_subNodeCoord                                                                                ! coordinates of subnodes per element
 
 logical, private :: noPart                                                                         ! for cases where the ABAQUS input file does not use part/assembly information
 

! Thee definitions should actually reside in the FE-solver specific part (different for MARC/ABAQUS)
! Hence, I suggest to prefix with "FE_"

 integer(pInt), parameter, public :: &
   FE_Nelemtypes       = 10_pInt, &
   FE_maxNnodes        = 8_pInt, &
   FE_maxNsubNodes     = 56_pInt, &
   FE_maxNips          = 27_pInt, &
   FE_maxNipNeighbors  = 6_pInt, &
   FE_maxmaxNnodesAtIP = 8_pInt, &                                                                  ! max number of (equivalent) nodes attached to an IP
   FE_NipFaceNodes     = 4_pInt
                      
 integer(pInt), dimension(FE_Nelemtypes), parameter, public :: FE_Nnodes = &                        ! nodes in a specific type of element (how we use it) 
 int([8, & ! element 7
      4, & ! element 134
      4, & ! element 11
      4, & ! element 27
      4, & ! element 157
      6, & ! element 136
      8, & ! element 21
      8, & ! element 117
      8, & ! element 57 (c3d20r == c3d8 --> copy of 7)
      3  & ! element 155, 125, 128
  ],pInt)
 integer(pInt), dimension(FE_Nelemtypes), parameter, public :: FE_Nips = &                          ! IPs in a specific type of element
 int([8, & ! element 7
      1, & ! element 134
      4, & ! element 11
      9, & ! element 27
      4, & ! element 157
      6, & ! element 136
      27,& ! element 21
      1, & ! element 117
      8, & ! element 57 (c3d20r == c3d8 --> copy of 7)
      3  & ! element 155, 125, 128
  ],pInt)
 integer(pInt), dimension(FE_Nelemtypes), parameter, public :: FE_NipNeighbors = &                  !IP neighbors in a specific type of element
 int([6, & ! element 7
      4, & ! element 134
      4, & ! element 11
      4, & ! element 27
      6, & ! element 157
      6, & ! element 136
      6, & ! element 21
      6, & ! element 117
      6, & ! element 57 (c3d20r == c3d8 --> copy of 7)
      4  & ! element 155, 125, 128
  ],pInt)
 integer(pInt), dimension(FE_Nelemtypes), parameter, private  :: FE_NsubNodes = &                   ! subnodes required to fully define all IP volumes
 int([19,& ! element 7                                                                              ! order is +x,-x,+y,-y,+z,-z but meaning strongly depends on Elemtype
      0, & ! element 134
      5, & ! element 11
      12,& ! element 27
      0, & ! element 157
      15,& ! element 136
      56,& ! element 21
      0, & ! element 117
      19,& ! element 57 (c3d20r == c3d8 --> copy of 7)
      4  & ! element 155, 125, 128
  ],pInt)
 integer(pInt), dimension(FE_Nelemtypes), parameter, private :: FE_NoriginalNodes = &               ! nodes in a specific type of element (how it is originally defined by marc)
 int([8, & ! element 7
      4, & ! element 134
      4, & ! element 11
      8, & ! element 27
      4, & ! element 157
      6, & ! element 136
      20,& ! element 21
      8, & ! element 117
      20,& ! element 57 (c3d20r == c3d8 --> copy of 7)
      6  & ! element 155, 125, 128
  ],pInt)
 integer(pInt), dimension(FE_maxNipNeighbors,FE_Nelemtypes), parameter, private :: FE_NfaceNodes = &! nodes per face in a specific type of element
 reshape(int([&
  4,4,4,4,4,4, & ! element 7
  3,3,3,3,0,0, & ! element 134
  2,2,2,2,0,0, & ! element 11
  2,2,2,2,0,0, & ! element 27
  3,3,3,3,0,0, & ! element 157
  3,4,4,4,3,0, & ! element 136
  4,4,4,4,4,4, & ! element 21
  4,4,4,4,4,4, & ! element 117
  4,4,4,4,4,4, & ! element 57 (c3d20r == c3d8 --> copy of 7)
  2,2,2,0,0,0  & ! element 155, 125, 128
  ],pInt),[FE_maxNipNeighbors,FE_Nelemtypes])   
 integer(pInt), dimension(FE_Nelemtypes), parameter, private :: FE_maxNnodesAtIP = &                ! map IP index to two node indices in a specific type of element
 int([1, & ! element 7
      4, & ! element 134
      1, & ! element 11
      2, & ! element 27
      1, & ! element 157
      1, & ! element 136
      4, & ! element 21
      8, & ! element 117
      1, & ! element 57 (c3d20r == c3d8 --> copy of 7)
      1  & ! element 155, 125, 128
  ],pInt)
 integer(pInt), dimension(FE_NipFaceNodes,FE_maxNipNeighbors,FE_Nelemtypes), parameter, private :: &
                                                                           FE_nodeOnFace = &        ! List of node indices on each face of a specific type of element
 reshape(int([&
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
  8,7,6,5 , &
  1,2,3,4 , & ! element 117
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,3,4 , & ! element 57 (c3d20r == c3d8 --> copy of 7)
  2,1,5,6 , &
  3,2,6,7 , &
  4,3,7,8 , &
  4,1,5,8 , &
  8,7,6,5 , &
  1,2,0,0 , & ! element 155,125,128
  2,3,0,0 , &
  3,1,0,0 , &
  0,0,0,0 , &
  0,0,0,0 , &
  0,0,0,0   &
  ],pInt),[FE_NipFaceNodes,FE_maxNipNeighbors,FE_Nelemtypes])
 
 integer(pInt), dimension(:,:,:), allocatable, private :: &
   FE_nodesAtIP, &                                                                                  ! map IP index to two node indices in a specific type of element
   FE_ipNeighbor, &                                                                                 ! +x,-x,+y,-y,+z,-z list of intra-element IPs and(negative) neighbor faces per own IP in a specific type of element
   FE_subNodeParent
  
 integer(pInt), dimension(:,:,:,:), allocatable, private :: &
   FE_subNodeOnIPFace
       
 public  :: mesh_init, &
            mesh_FEasCP, &
            mesh_build_subNodeCoords, &
            mesh_build_ipVolumes, &
            mesh_build_ipCoordinates, &
            mesh_regrid
 private :: FE_mapElemtype, &
            mesh_faceMatch, &
            mesh_build_FEdata, &
            mesh_marc_get_tableStyles, &
            mesh_get_damaskOptions, &
            mesh_spectral_count_nodesAndElements, &
            mesh_marc_count_nodesAndElements, &
            mesh_abaqus_count_nodesAndElements, &
            mesh_abaqus_count_elementSets, &
            mesh_abaqus_count_materials, &
            mesh_spectral_count_cpElements, &
            mesh_abaqus_count_cpElements, &
            mesh_marc_map_elementSets, &
            mesh_abaqus_map_elementSets, &
            mesh_abaqus_map_materials, &
            mesh_spectral_map_nodes, &
            mesh_marc_map_nodes, &
            mesh_abaqus_map_nodes, &
            mesh_marc_map_elements, &
            mesh_abaqus_map_elements, &
            mesh_spectral_count_cpSizes, &
            mesh_marc_count_cpSizes, &
            mesh_abaqus_count_cpSizes, &
            mesh_spectral_build_nodes, &
            mesh_marc_build_nodes, &
            mesh_abaqus_build_nodes, &
            mesh_spectral_build_elements, &
            mesh_marc_build_elements, &
            mesh_abaqus_build_elements, &
            mesh_build_ipNeighborhood, &
            mesh_build_ipAreas, &
            mesh_build_nodeTwins, &
            mesh_tell_statistics
contains


!***********************************************************
! initialization 
!***********************************************************
subroutine mesh_init(ip,element)

 use DAMASK_interface
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO,        only: IO_error, &
                      IO_open_InputFile, &
                      IO_abaqus_hasNoPart
 use FEsolving, only: parallelExecution, &
                      FEsolving_execElem, &
                      FEsolving_execIP, &
                      calcMode, &
                      lastMode, &
                      FEmodelGeometry
 
 implicit none
 integer(pInt), parameter :: fileUnit = 222_pInt
 integer(pInt) :: e, element, ip
 
 !$OMP CRITICAL (write2out)
   write(6,*)
   write(6,*) '<<<+-  mesh init  -+>>>'
   write(6,*) '$Id$'
#include "compilation_info.f90"
 !$OMP END CRITICAL (write2out)

 call mesh_build_FEdata                                                                             ! get properties of the different types of elements

 call IO_open_inputFile(fileUnit,FEmodelGeometry)                                                   ! parse info from input file...

 select case (FEsolver)
   case ('Spectral')
                       call mesh_spectral_count_nodesAndElements(fileUnit)
                       call mesh_spectral_count_cpElements
                       call mesh_spectral_map_elements
                       call mesh_spectral_map_nodes
                       call mesh_spectral_count_cpSizes
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
                       noPart = IO_abaqus_hasNoPart(fileUnit)
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
 call mesh_get_damaskOptions(fileUnit)
 close (fileUnit)
 
 call mesh_build_subNodeCoords
 call mesh_build_ipCoordinates
 call mesh_build_ipVolumes
 call mesh_build_ipAreas
 call mesh_build_nodeTwins
 call mesh_build_sharedElems
 call mesh_build_ipNeighborhood
 call mesh_tell_statistics

 parallelExecution = (parallelExecution .and. (mesh_Nelems == mesh_NcpElems))      ! plus potential killer from non-local constitutive
 
 FEsolving_execElem = [ 1_pInt,mesh_NcpElems]
 allocate(FEsolving_execIP(2_pInt,mesh_NcpElems)); FEsolving_execIP = 1_pInt
 forall (e = 1_pInt:mesh_NcpElems) FEsolving_execIP(2,e) = FE_Nips(mesh_element(2,e))
 
 allocate(calcMode(mesh_maxNips,mesh_NcpElems))
 calcMode = .false.                                       ! pretend to have collected what first call is asking (F = I)
 calcMode(ip,mesh_FEasCP('elem',element)) = .true.        ! first ip,el needs to be already pingponged to "calc"
 lastMode = .true.                                        ! and its mode is already known...

end subroutine mesh_init
 

!***********************************************************
! mapping of FE element types to internal representation
!***********************************************************
integer(pInt) function FE_mapElemtype(what)
 
 use IO, only: IO_lc

 implicit none
 character(len=*), intent(in) :: what
  
 select case (IO_lc(what))
    case (  '7', &
            'c3d8')
      FE_mapElemtype = 1_pInt            ! Three-dimensional Arbitrarily Distorted Brick
    case ('134', &
          'c3d4')
      FE_mapElemtype = 2_pInt            ! Three-dimensional Four-node Tetrahedron
    case ( '11', &
           'cpe4')
      FE_mapElemtype = 3_pInt            ! Arbitrary Quadrilateral Plane-strain
    case ( '27', &
           'cpe8')
      FE_mapElemtype = 4_pInt           ! Plane Strain, Eight-node Distorted Quadrilateral
    case ('157')
      FE_mapElemtype = 5_pInt            ! Three-dimensional, Low-order, Tetrahedron, Herrmann Formulations
    case ('136', &
          'c3d6')
      FE_mapElemtype = 6_pInt            ! Three-dimensional Arbitrarily Distorted Pentahedral
    case ( '21', &
           'c3d20')
      FE_mapElemtype = 7_pInt            ! Three-dimensional Arbitrarily Distorted quadratic hexahedral
    case ( '117', &
           '123', &
           'c3d8r')
      FE_mapElemtype = 8_pInt            ! Three-dimensional Arbitrarily Distorted linear hexahedral with reduced integration
    case ( '57', &
           'c3d20r')
      FE_mapElemtype = 9_pInt            ! Three-dimensional Arbitrarily Distorted quad hexahedral with reduced integration
    case ( '155', &
           '125', &
           '128')
      FE_mapElemtype = 10_pInt           ! Two-dimensional Plane Strain triangle (155: cubic shape function, 125/128: second order isoparametric)
    case default 
      FE_mapElemtype = 0_pInt            ! unknown element --> should raise an error upstream..!
 endselect

 end function FE_mapElemtype



!***********************************************************
! FE to CP id mapping by binary search thru lookup array
!
! valid questions are 'elem', 'node'
!***********************************************************
integer(pInt) function mesh_FEasCP(what,myID)
 
 use IO, only: IO_lc

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
 
 ! check at bounds QUESTION is it valid to extend bounds by 1 and just do binary search w/o init check at bounds?
 if (lookupMap(1_pInt,lower) == myID) then
   mesh_FEasCP = lookupMap(2_pInt,lower)
   return
 elseif (lookupMap(1_pInt,upper) == myID) then
   mesh_FEasCP = lookupMap(2_pInt,upper)
   return
 endif
 
 ! binary search in between bounds
 do while (upper-lower > 1_pInt)
   center = (lower+upper)/2_pInt
   if (lookupMap(1_pInt,center) < myID) then
     lower = center
   elseif (lookupMap(1_pInt,center) > myID) then
     upper = center
   else
     mesh_FEasCP = lookupMap(2_pInt,center)
     exit
   endif
 enddo
 
end function mesh_FEasCP


!***********************************************************
! find face-matching element of same type
!***********************************************************
subroutine mesh_faceMatch(elem, face ,matchingElem, matchingFace)

implicit none
!*** output variables
integer(pInt), intent(out) ::     matchingElem, &                                                   ! matching CP element ID
                                  matchingFace                                                      ! matching FE face ID 

!*** input variables
integer(pInt), intent(in) ::      face, &                                                           ! FE face ID
                                  elem                                                              ! FE elem ID

!*** local variables
integer(pInt), dimension(FE_NfaceNodes(face,mesh_element(2,elem))) :: &
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
myType = mesh_element(2_pInt,elem)                                                                  ! figure elemType

do n = 1_pInt,FE_NfaceNodes(face,myType)                                                            ! loop over nodes on face
  myFaceNodes(n) = mesh_FEasCP('node',mesh_element(4_pInt+FE_nodeOnFace(n,face,myType),elem))       ! CP id of face node
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
    candidateType = mesh_element(2_pInt,candidateElem)                                              ! figure elemType of candidate
checkCandidateFace: do candidateFace = 1_pInt,FE_maxNipNeighbors                                    ! check each face of candidate
      if (FE_NfaceNodes(candidateFace,candidateType) /= FE_NfaceNodes(face,myType) & ! incompatible face
          .or. (candidateElem == elem .and. candidateFace == face)) then  ! this is my face
        cycle checkCandidateFace
      endif
      checkTwins = .false.
      do n = 1_pInt,FE_NfaceNodes(candidateFace,candidateType)           ! loop through nodes on face
        candidateFaceNode = mesh_FEasCP('node', mesh_element(4_pInt+FE_nodeOnFace(n,candidateFace,candidateType),candidateElem))
        if (all(myFaceNodes /= candidateFaceNode)) then             ! candidate node does not match any of my face nodes
          checkTwins = .true.                                       ! perhaps the twin nodes do match
          exit
        endif
      enddo
      if(checkTwins) then
checkCandidateFaceTwins: do dir = 1_pInt,3_pInt
          do n = 1_pInt,FE_NfaceNodes(candidateFace,candidateType)       ! loop through nodes on face
            candidateFaceNode = mesh_FEasCP('node', mesh_element(4+FE_nodeOnFace(n,candidateFace,candidateType),candidateElem))
            if (all(myFaceNodes /= mesh_nodeTwins(dir,candidateFaceNode))) then ! node twin does not match either
              if (dir == 3_pInt) then
                cycle checkCandidateFace
              else
                cycle checkCandidateFaceTwins                       ! try twins in next dimension
              endif
            endif
          enddo
          exit checkCandidateFaceTwins
        enddo checkCandidateFaceTwins
      endif
      matchingFace = candidateFace
      matchingElem = candidateElem
      exit checkCandidate                                           ! found my matching candidate
    enddo checkCandidateFace
  endif
enddo checkCandidate

deallocate(element_seen)

end subroutine mesh_faceMatch

 
!********************************************************************
! get properties of different types of finite elements
!
! assign globals:
! FE_nodesAtIP, FE_ipNeighbor, FE_subNodeParent, FE_subNodeOnIPFace
!********************************************************************
subroutine mesh_build_FEdata

 implicit none
 allocate(FE_nodesAtIP(FE_maxmaxNnodesAtIP,FE_maxNips,FE_Nelemtypes)) ; FE_nodesAtIP = 0_pInt
 allocate(FE_ipNeighbor(FE_maxNipNeighbors,FE_maxNips,FE_Nelemtypes)) ; FE_ipNeighbor = 0_pInt
 allocate(FE_subNodeParent(FE_maxNips,FE_maxNsubNodes,FE_Nelemtypes)) ; FE_subNodeParent = 0_pInt
 allocate(FE_subNodeOnIPFace(FE_NipFaceNodes,FE_maxNipNeighbors,FE_maxNips,FE_Nelemtypes)) ; FE_subNodeOnIPFace = 0_pInt
 
 ! fill FE_nodesAtIP with data
 FE_nodesAtIP(1:FE_maxNnodesAtIP(1),1:FE_Nips(1),1) = &  ! element 7
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3,  &
    5,  &
    6,  &
    8,  &
    7   &
    ],pInt),[FE_maxNnodesAtIP(1),FE_Nips(1)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(2),1:FE_Nips(2),2) = &  ! element 134
    reshape(int([&
    1,2,3,4   &
    ],pInt),[FE_maxNnodesAtIP(2),FE_Nips(2)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(3),1:FE_Nips(3),3) = &  ! element 11
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(3),FE_Nips(3)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(4),1:FE_Nips(4),4) = &  ! element 27
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
    ],pInt),[FE_maxNnodesAtIP(4),FE_Nips(4)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(5),1:FE_Nips(5),5) = &  ! element 157
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4   &
    ],pInt),[FE_maxNnodesAtIP(5),FE_Nips(5)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(6),1:FE_Nips(6),6) = &  ! element 136
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4,  &
    5,  &
    6   &
    ],pInt),[FE_maxNnodesAtIP(6),FE_Nips(6)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(7),1:FE_Nips(7),7) = &  ! element 21
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
    ],pInt),[FE_maxNnodesAtIP(7),FE_Nips(7)])

! FE_nodesAtIP(:,:FE_Nips(8),8) = &  ! element 117 (c3d8r --> single IP per element, so no need for this mapping)
!    reshape((/&
!    1,2,3,4,5,6,7,8   &
!    /),(/FE_maxNnodesAtIP(8),FE_Nips(8)/))

 FE_nodesAtIP(1:FE_maxNnodesAtIP(9),1:FE_Nips(9),9) = &  ! element 57 (c3d20r == c3d8 --> copy of 7)
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3,  &
    5,  &
    6,  &
    8,  &
    7   &
    ],pInt),[FE_maxNnodesAtIP(9),FE_Nips(9)])

 FE_nodesAtIP(1:FE_maxNnodesAtIP(10),1:FE_Nips(10),10) = &  ! element 155, 125, 128
    reshape(int([&
    1,  &
    2,  &
    3   &
    ],pInt),[FE_maxNnodesAtIP(10),FE_Nips(10)])
 
 ! *** FE_ipNeighbor ***
 ! is a list of the neighborhood of each IP.
 ! It is sorted in (local) +x,-x, +y,-y, +z,-z direction.
 ! Positive integers denote an intra-FE IP identifier.
 ! Negative integers denote the interface behind which the neighboring (extra-FE) IP will be located.

 FE_ipNeighbor(1:FE_NipNeighbors(1),1:FE_Nips(1),1) = &  ! element 7
    reshape(int([&
     2,-5, 3,-2, 5,-1,  &
    -3, 1, 4,-2, 6,-1,  &
     4,-5,-4, 1, 7,-1,  &
    -3, 3,-4, 2, 8,-1,  &
     6,-5, 7,-2,-6, 1,  &
    -3, 5, 8,-2,-6, 2,  &
     8,-5,-4, 5,-6, 3,  &
    -3, 7,-4, 6,-6, 4   &
    ],pInt),[FE_NipNeighbors(1),FE_Nips(1)])

 FE_ipNeighbor(1:FE_NipNeighbors(2),1:FE_Nips(2),2) = &  ! element 134
    reshape(int([&
    -1,-2,-3,-4   &
    ],pInt),[FE_NipNeighbors(2),FE_Nips(2)])

 FE_ipNeighbor(1:FE_NipNeighbors(3),1:FE_Nips(3),3) = &  ! element 11
    reshape(int([&
     2,-4, 3,-1,  &
    -2, 1, 4,-1,  &
     4,-4,-3, 1,  &
    -2, 3,-3, 2   &
    ],pInt),[FE_NipNeighbors(3),FE_Nips(3)])

 FE_ipNeighbor(1:FE_NipNeighbors(4),1:FE_Nips(4),4) = &  ! element 27
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
    ],pInt),[FE_NipNeighbors(4),FE_Nips(4)])

 FE_ipNeighbor(1:FE_NipNeighbors(5),1:FE_Nips(5),5) = &  ! element 157
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
     3,-2, 1,-3, 4,-1,  &
     1,-3, 2,-4, 4,-1,  &
     1,-3, 2,-4, 3,-2   &
    ],pInt),[FE_NipNeighbors(5),FE_Nips(5)])

 FE_ipNeighbor(1:FE_NipNeighbors(6),1:FE_Nips(6),6) = &  ! element 136
    reshape(int([&
     2,-4, 3,-2, 4,-1,  &
    -3, 1, 3,-2, 5,-1,  &
     2,-4,-3, 1, 6,-1,  &
     5,-4, 6,-2,-5, 1,  &
    -3, 4, 6,-2,-5, 2,  &
     5,-4,-3, 4,-5, 3   &
    ],pInt),[FE_NipNeighbors(6),FE_Nips(6)])

 FE_ipNeighbor(1:FE_NipNeighbors(7),1:FE_Nips(7),7) = &  ! element 21
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
    ],pInt),[FE_NipNeighbors(7),FE_Nips(7)])

FE_ipNeighbor(1:FE_NipNeighbors(8),1:FE_Nips(8),8) = &  ! element 117
    reshape(int([&
    -3,-5,-4,-2,-6,-1   &
    ],pInt),[FE_NipNeighbors(8),FE_Nips(8)])

 FE_ipNeighbor(1:FE_NipNeighbors(9),1:FE_Nips(9),9) = &  ! element 57 (c3d20r == c3d8 --> copy of 7)
    reshape(int([&
     2,-5, 3,-2, 5,-1,  &
    -3, 1, 4,-2, 6,-1,  &
     4,-5,-4, 1, 7,-1,  &
    -3, 3,-4, 2, 8,-1,  &
     6,-5, 7,-2,-6, 1,  &
    -3, 5, 8,-2,-6, 2,  &
     8,-5,-4, 5,-6, 3,  &
    -3, 7,-4, 6,-6, 4   &
    ],pInt),[FE_NipNeighbors(9),FE_Nips(9)])

 FE_ipNeighbor(1:FE_NipNeighbors(10),1:FE_Nips(10),10) = &  ! element 155, 125, 128
    reshape(int([&
     2,-3, 3,-1,  &
    -2, 1, 3,-1,  &
     2,-3,-2, 1   &
    ],pInt),[FE_NipNeighbors(10),FE_Nips(10)])
 
 ! *** FE_subNodeParent ***
 ! lists the group of nodes for which the center of gravity
 ! corresponds to the location of a each subnode.
 ! example: face-centered subnode with faceNodes 1,2,3,4 to be used in,
 !          e.g., a 8 IP grid, would be encoded:
 !          1, 2, 3, 4, 0, 0, 0, 0

 FE_subNodeParent(1:FE_Nips(1),1:FE_NsubNodes(1),1) = &  ! element 7
    reshape(int([&
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
    ],pInt),(/FE_Nips(1),FE_NsubNodes(1)/))

!FE_subNodeParent(:FE_Nips(2),:FE_NsubNodes(2),2)      ! element 134 has no subnodes

 FE_subNodeParent(1:FE_Nips(3),1:FE_NsubNodes(3),3) = &  ! element 11
    reshape(int([&
    1, 2, 0, 0,  & 
    2, 3, 0, 0,  & 
    3, 4, 0, 0,  &
    4, 1, 0, 0,  & 
    1, 2, 3, 4   &
    ],pInt),[FE_Nips(3),FE_NsubNodes(3)])

 FE_subNodeParent(1:FE_Nips(4),1:FE_NsubNodes(4),4) = &  ! element 27
    reshape(int([&
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
    ],pInt),[FE_Nips(4),FE_NsubNodes(4)])

 !FE_subNodeParent(:FE_Nips(5),:FE_NsubNodes(5),5) = &  ! element 157
 !   reshape((/&
 !   *still to be defined*
 !   ],pInt),(/FE_Nips(5),FE_NsubNodes(5)/))

 FE_subNodeParent(1:FE_Nips(6),1:FE_NsubNodes(6),6) = &  ! element 136
    reshape(int([&
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
    ],pInt),[FE_Nips(6),FE_NsubNodes(6)])

 FE_subNodeParent(1:FE_Nips(7),1:FE_NsubNodes(7),7) = &  ! element 21
    reshape(int([&
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
    ],pInt),[FE_Nips(7),FE_NsubNodes(7)])

!FE_subNodeParent(:FE_Nips(8),:FE_NsubNodes(8),8)      ! element 117 has no subnodes

 FE_subNodeParent(1:FE_Nips(9),1:FE_NsubNodes(9),9) = &  ! element 57 (c3d20r == c3d8 --> copy of 7)
    reshape(int([&
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
    ],pInt),[FE_Nips(9),FE_NsubNodes(9)])

 FE_subNodeParent(1:FE_Nips(10),1:FE_NsubNodes(10),10) = &  ! element 155, 125, 128
    reshape(int([&
    1, 2, 0,  & 
    2, 3, 0,  & 
    3, 1, 0,  &
    1, 2, 3   &
    ],pInt),[FE_Nips(10),FE_NsubNodes(10)])
 
 ! *** FE_subNodeOnIPFace ***
 ! indicates which subnodes make up the interfaces enclosing the IP volume.
 ! The sorting convention is such that the outward pointing normal
 ! follows from a right-handed traversal of the face node list.
 ! For two-dimensional elements, which only have lines as "interface"
 ! one nevertheless has to specify each interface by a closed path,
 ! e.g., 1,2, 2,1, assuming the line connects nodes 1 and 2.
 ! This will result in zero ipVolume and interfaceArea, but is not
 ! detrimental at the moment since non-local constitutive laws are
 ! currently not foreseen in 2D cases.
 
 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(1),1:FE_Nips(1),1) = &  ! element 7
    reshape(int([&
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
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(1),FE_Nips(1)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(2),1:FE_Nips(2),2) = &  ! element 134
    reshape(int([&
     1, 1, 3, 2, & ! 1
     1, 1, 2, 4, &
     2, 2, 3, 4, &
     1, 1, 4, 3  &
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(2),FE_Nips(2)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(3),1:FE_Nips(3),3) = &  ! element 11
    reshape(int([&
     5, 9, 9, 5 , & ! 1
     1, 8, 8, 1 , &
     8, 9, 9, 8 , &
     1, 5, 5, 1 , &
     2, 6, 6, 2 , & ! 2
     5, 9, 9, 5 , &
     6, 9, 9, 6 , &
     2, 5, 5, 2 , &
     3, 6, 6, 3 , & ! 3
     7, 9, 9, 7 , &
     3, 7, 7, 3 , &
     6, 9, 9, 6 , &
     7, 9, 9, 7 , & ! 4
     4, 8, 8, 4 , &
     4, 7, 7, 4 , &
     8, 9, 9, 8   &
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(3),FE_Nips(3)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(4),1:FE_Nips(4),4) = &  ! element 27
    reshape(int([&
     9,17,17, 9 , & ! 1
     1,16,16, 1 , &
    16,17,17,16 , &
     1, 9, 9, 1 , &
    10,18,18,10 , & ! 2
     9,17,17, 9 , &
    17,18,18,17 , &
     9,10,10, 9 , &
     2,11,11, 2 , & ! 3
    10,18,18,10 , &
    11,18,18,11 , &
     2,10,10, 2 , &
    17,20,20,17 , & ! 4
    15,16,16,15 , &
    15,20,20,15 , &
    16,17,17,16 , &
    18,19,19,18 , & ! 5
    17,20,20,17 , &
    19,20,20,19 , &
    17,18,18,17 , &
    11,12,12,11 , & ! 6
    18,19,19,18 , &
    12,19,19,12 , &
    11,18,18,11 , &
    14,20,20,14 , & ! 7
     4,15,15, 4 , &
     4,14,14, 4 , &
    15,20,20,15 , &
    13,19,19,13 , & ! 8
    14,20,20,14 , &
    13,14,14,13 , &
    19,20,20,19 , &
     3,12,12, 3 , & ! 9
    13,19,19,13 , &
     3,13,13, 3 , &
    12,19,19,12   &
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(4),FE_Nips(4)])

 !FE_subNodeOnIPFace(:FE_NipFaceNodes,:FE_NipNeighbors(5),:FE_Nips(5),5) = &  ! element 157
 !   reshape((/&
 !   *still to be defined*
 !   /),(/FE_NipFaceNodes,FE_NipNeighbors(5),FE_Nips(5)/))

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(6),1:FE_Nips(6),6) = &  ! element 136
    reshape(int([&
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
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(6),FE_Nips(6)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(7),1:FE_Nips(7),7) = &  ! element 21
    reshape(int([&
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
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(7),FE_Nips(7)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(8),1:FE_Nips(8),8) = &  ! element 117
    reshape(int([&
     2, 3, 7, 6, & ! 1
     1, 5, 8, 4, &
     3, 4, 8, 7, &
     1, 2, 6, 5, &
     5, 6, 7, 8, &
     1, 4, 3, 2  &
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(8),FE_Nips(8)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(9),1:FE_Nips(9),9) = &  ! element 57 (c3d20r == c3d8 --> copy of 7)
    reshape(int([&
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
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(9),FE_Nips(9)])

 FE_subNodeOnIPFace(1:FE_NipFaceNodes,1:FE_NipNeighbors(10),1:FE_Nips(10),10) = &  ! element 155, 125, 128
    reshape(int([&
     4, 7, 7, 4 , & ! 1
     1, 6, 6, 1 , &
     6, 7, 7, 6 , &
     1, 4, 4, 1 , &
     2, 5, 5, 2 , & ! 2
     4, 7, 7, 4 , &
     5, 7, 7, 5 , &
     2, 4, 4, 2 , &
     5, 7, 7, 5 , & ! 3
     3, 6, 6, 3 , &
     3, 5, 5, 3 , &
     6, 7, 7, 6   &
    ],pInt),[FE_NipFaceNodes,FE_NipNeighbors(10),FE_Nips(10)])

end subroutine mesh_build_FEdata


!********************************************************************
! figure out table styles (Marc only)
!
! initialcondTableStyle, hypoelasticTableStyle
!********************************************************************
subroutine mesh_marc_get_tableStyles(myUnit)

 use IO,   only: IO_lc, &
                 IO_intValue, &
                 IO_stringValue, &
                 IO_stringPos
 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 6_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line

 initialcondTableStyle = 0_pInt
 hypoelasticTableStyle = 0_pInt
 
610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'table' .and. myPos(1_pInt) .GT. 5) then
     initialcondTableStyle = IO_intValue(line,myPos,4_pInt)
     hypoelasticTableStyle = IO_intValue(line,myPos,5_pInt)
     exit
   endif
 enddo

620 end subroutine mesh_marc_get_tableStyles


!********************************************************************
! get any additional damask options from input file
!
! mesh_periodicSurface
!********************************************************************
subroutine mesh_get_damaskOptions(myUnit)

use DAMASK_interface, only: FEsolver
use IO,  only: IO_lc, &
                IO_stringValue, &
                IO_stringPos

implicit none
integer(pInt), intent(in) :: myUnit

integer(pInt), parameter :: maxNchunks = 5_pInt
integer(pInt), dimension (1+2*maxNchunks) :: myPos
integer(pInt) chunk, Nchunks
character(len=300) line, keyword, damaskOption, v

mesh_periodicSurface = (/.false., .false., .false./)

610 FORMAT(A300)

select case (FEsolver)
  case ('Spectral') ! no special keyword needed, the damask option directly goes into the header
  case ('Marc')
    keyword = '$damask'
  case ('Abaqus')
    keyword = '**damask'
end select

rewind(myUnit)
do 
  read (myUnit,610,END=620) line
  myPos = IO_stringPos(line,maxNchunks)
  Nchunks = myPos(1)
  select case (FEsolver)
    case ('Marc','Abaqus')
      if (IO_lc(IO_stringValue(line,myPos,1_pInt)) == keyword .and. Nchunks > 1_pInt) then  ! found keyword for damask option and there is at least one more chunk to read
        damaskOption = IO_lc(IO_stringValue(line,myPos,2_pInt))
        select case(damaskOption)
          case('periodic')                                                        ! damask Option that allows to specify periodic fluxes
            do chunk = 3_pInt,Nchunks                                                  ! loop through chunks (skipping the keyword)
              v = IO_lc(IO_stringValue(line,myPos,chunk))                       ! chunk matches keyvalues x,y, or z?
              mesh_periodicSurface(1) = mesh_periodicSurface(1) .or. v == 'x'
              mesh_periodicSurface(2) = mesh_periodicSurface(2) .or. v == 'y'
              mesh_periodicSurface(3) = mesh_periodicSurface(3) .or. v == 'z'
            enddo
        endselect
      endif
    case('Spectral')
      damaskOption = IO_lc(IO_stringValue(line,myPos,1_pInt))
      select case(damaskOption)
        case('periodic')                                                        ! damask Option that allows to specify periodic fluxes
          do chunk = 2_pInt,Nchunks                                                  ! loop through chunks (skipping the keyword)
            v = IO_lc(IO_stringValue(line,myPos,chunk))                       ! chunk matches keyvalues x,y, or z?
            mesh_periodicSurface(1) = mesh_periodicSurface(1) .or. v == 'x'
            mesh_periodicSurface(2) = mesh_periodicSurface(2) .or. v == 'y'
            mesh_periodicSurface(3) = mesh_periodicSurface(3) .or. v == 'z'
          enddo
      endselect
  endselect
enddo

620 end subroutine mesh_get_damaskOptions

 
!********************************************************************
! count overall number of nodes and elements in mesh
!
! mesh_Nelems, mesh_Nnodes
!********************************************************************
subroutine mesh_spectral_count_nodesAndElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_intValue, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error
 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 7_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 integer(pInt) :: a = 0_pInt, &
                  b = 0_pInt, &
                  c = 0_pInt, &
                  headerLength = 0_pInt, &
                  i,j
 character(len=1024) line,keyword

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt
 
 rewind(myUnit)
 read(myUnit,'(a1024)') line
 myPos = IO_stringPos(line,2_pInt)
 keyword = IO_lc(IO_StringValue(line,myPos,2_pInt))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,myPos,1_pInt) + 1_pInt 
 else
   call IO_error(error_ID=842_pInt)
 endif
 
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   myPos = IO_stringPos(line,maxNchunks)             
   if ( IO_lc(IO_StringValue(line,myPos,1_pInt)) == 'resolution') then
     do j = 2_pInt,6_pInt,2_pInt
       select case (IO_lc(IO_stringValue(line,myPos,j)))
         case('a')
           a = IO_intValue(line,myPos,j+1_pInt)
         case('b')
           b = IO_intValue(line,myPos,j+1_pInt)
         case('c')
           c = IO_intValue(line,myPos,j+1_pInt)
       end select
     enddo
     mesh_Nelems = a * b * c
     mesh_Nnodes = (1_pInt + a)*(1_pInt + b)*(1_pInt + c)
   endif
 enddo

end subroutine mesh_spectral_count_nodesAndElements

!********************************************************************
! count overall number of nodes and elements in mesh
!
! mesh_Nelems, mesh_Nnodes
!********************************************************************
subroutine mesh_marc_count_nodesAndElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_IntValue
 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,myPos,1_pInt)) == 'sizing') then
       mesh_Nelems = IO_IntValue (line,myPos,3_pInt)
       mesh_Nnodes = IO_IntValue (line,myPos,4_pInt)
     exit
   endif
 enddo

620 end subroutine mesh_marc_count_nodesAndElements

!********************************************************************
! count overall number of nodes and elements in mesh
!
! mesh_Nelems, mesh_Nnodes
!********************************************************************
subroutine mesh_abaqus_count_nodesAndElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_error
                 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) :: line
 logical :: inPart

 mesh_Nnodes = 0_pInt
 mesh_Nelems = 0_pInt
 
610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.
   
   if (inPart .or. noPart) then
     select case ( IO_lc(IO_stringValue(line,myPos,1_pInt)))
       case('*node')
          if( &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'print'    .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'file'     .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' &
             ) &
            mesh_Nnodes = mesh_Nnodes + IO_countDataLines(myUnit)
       case('*element')
          if( &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'matrix'   .and. &
              IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' &
             ) then
            mesh_Nelems = mesh_Nelems + IO_countDataLines(myUnit)
          endif
     endselect
   endif
 enddo
 
620 if (mesh_Nnodes < 2_pInt)  call IO_error(error_ID=900_pInt)
 if (mesh_Nelems == 0_pInt) call IO_error(error_ID=901_pInt)
 
end subroutine mesh_abaqus_count_nodesAndElements

 
!********************************************************************
! count overall number of element sets in mesh
!
! mesh_NelemSets, mesh_maxNelemInSet
!********************************************************************
 subroutine mesh_marc_count_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countContinousIntValues
                 
 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line

 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_StringValue(line,myPos,1_pInt)) == 'define' .and. &
        IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'element' ) then
     mesh_NelemSets = mesh_NelemSets + 1_pInt
     mesh_maxNelemInSet = max(mesh_maxNelemInSet, &
                              IO_countContinousIntValues(myUnit))
   endif
 enddo

620 end subroutine mesh_marc_count_elementSets


!********************************************************************
! count overall number of element sets in mesh
!
! mesh_NelemSets, mesh_maxNelemInSet
!********************************************************************
subroutine mesh_abaqus_count_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) :: line
 logical :: inPart
 
 mesh_NelemSets     = 0_pInt
 mesh_maxNelemInSet = mesh_Nelems               ! have to be conservative, since Abaqus allows for recursive definitons
 
610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.
   
   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*elset' ) &
     mesh_NelemSets = mesh_NelemSets + 1_pInt
 enddo

620 continue
 if (mesh_NelemSets == 0) call IO_error(error_ID=902_pInt)

end subroutine mesh_abaqus_count_elementSets


!********************************************************************
! count overall number of solid sections sets in mesh (Abaqus only)
!
! mesh_Nmaterials
!********************************************************************
subroutine mesh_abaqus_count_materials(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 logical inPart
 
 mesh_Nmaterials = 0_pInt
 
610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,myPos,1_pInt)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'section' ) &
     mesh_Nmaterials = mesh_Nmaterials + 1_pInt
 enddo

620 if (mesh_Nmaterials == 0_pInt) call IO_error(error_ID=903_pInt)
 
end subroutine mesh_abaqus_count_materials


!********************************************************************
! count overall number of cpElements in mesh
!
! mesh_NcpElems
!********************************************************************
subroutine mesh_spectral_count_cpElements

 implicit none

 mesh_NcpElems = mesh_Nelems
 
end subroutine mesh_spectral_count_cpElements
 

!********************************************************************
! count overall number of cpElements in mesh
!
! mesh_NcpElems
!********************************************************************
subroutine mesh_marc_count_cpElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countContinousIntValues
                 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 integer(pInt) :: i
 character(len=300):: line

 mesh_NcpElems = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)

   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'hypoelastic') then
       do i=1_pInt,3_pInt+hypoelasticTableStyle  ! Skip 3 or 4 lines
         read (myUnit,610,END=620) line
       enddo
       mesh_NcpElems = mesh_NcpElems + IO_countContinousIntValues(myUnit)
     exit
   endif
 enddo

620 end subroutine mesh_marc_count_cpElements
 

!********************************************************************
! count overall number of cpElements in mesh
!
! mesh_NcpElems
!********************************************************************
subroutine mesh_abaqus_count_cpElements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_extractValue
                 
 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) line
 integer(pInt) :: i,k
 logical :: materialFound = .false.
 character(len=64) ::materialName,elemSetName
 
 mesh_NcpElems = 0_pInt
 
610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_stringValue(line,myPos,1_pInt)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'name'))     ! extract name=value
       materialFound = materialName /= ''                                           ! valid name?
     case('*user')
       if (IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'material' .and. materialFound) then
         do i = 1_pInt,mesh_Nmaterials                                                   ! look thru material names
           if (materialName == mesh_nameMaterial(i)) then                           ! found one
             elemSetName = mesh_mapMaterial(i)                                      ! take corresponding elemSet
             do k = 1_pInt,mesh_NelemSets                                                ! look thru all elemSet definitions
               if (elemSetName == mesh_nameElemSet(k)) &                            ! matched?
                 mesh_NcpElems = mesh_NcpElems + mesh_mapElemSet(1,k)               ! add those elem count
             enddo
           endif
         enddo
         materialFound = .false.
       endif
   endselect
 enddo
 
620 if (mesh_NcpElems == 0_pInt) call IO_error(error_ID=906_pInt)

end subroutine mesh_abaqus_count_cpElements


!********************************************************************
! map element sets
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!********************************************************************
subroutine mesh_marc_map_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_continousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1+2*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: elemSet = 0_pInt

 allocate (mesh_nameElemSet(mesh_NelemSets))                     ; mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1_pInt+mesh_maxNelemInSet,mesh_NelemSets)) ; mesh_mapElemSet = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=640) line
   myPos = IO_stringPos(line,maxNchunks)
   if( (IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'define' ) .and. &
       (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'element' ) ) then
      elemSet = elemSet+1_pInt
      mesh_nameElemSet(elemSet) = trim(IO_stringValue(line,myPos,4_pInt))
      mesh_mapElemSet(:,elemSet) = IO_continousIntValues(myUnit,mesh_maxNelemInSet,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
   endif
 enddo
 
640 end subroutine mesh_marc_map_elementSets


!********************************************************************
! Build element set mapping 
!
! allocate globals: mesh_nameElemSet, mesh_mapElemSet
!********************************************************************
subroutine mesh_abaqus_map_elementSets(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_continousIntValues, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: elemSet = 0_pInt,i
 logical :: inPart = .false.

 allocate (mesh_nameElemSet(mesh_NelemSets))                          ; mesh_nameElemSet = ''
 allocate (mesh_mapElemSet(1_pInt+mesh_maxNelemInSet,mesh_NelemSets)) ; mesh_mapElemSet  = 0_pInt

610 FORMAT(A300)


 rewind(myUnit)
 do
   read (myUnit,610,END=640) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.
   
   if ( (inPart .or. noPart) .and. IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*elset' ) then
     elemSet = elemSet + 1_pInt
     mesh_nameElemSet(elemSet)  = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'elset'))
     mesh_mapElemSet(:,elemSet) = IO_continousIntValues(myUnit,mesh_Nelems,mesh_nameElemSet,&
                                          mesh_mapElemSet,elemSet-1_pInt)
   endif
 enddo

640 do i = 1_pInt,elemSet
!   write(6,*)'elemSetName: ',mesh_nameElemSet(i)
!   write(6,*)'elems in Elset',mesh_mapElemSet(:,i)
   if (mesh_mapElemSet(1,i) == 0_pInt) call IO_error(error_ID=904_pInt,ext_msg=mesh_nameElemSet(i))
 enddo

end subroutine mesh_abaqus_map_elementSets


!********************************************************************
! map solid section (Abaqus only)
!
! allocate globals: mesh_nameMaterial, mesh_mapMaterial
!********************************************************************
subroutine mesh_abaqus_map_materials(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 20_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt) :: i,c = 0_pInt
 logical :: inPart = .false.
 character(len=64) :: elemSetName,materialName
 
 allocate (mesh_nameMaterial(mesh_Nmaterials)) ; mesh_nameMaterial = ''
 allocate (mesh_mapMaterial(mesh_Nmaterials)) ;  mesh_mapMaterial = ''

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if ( (inPart .or. noPart) .and. &
        IO_lc(IO_StringValue(line,myPos,1_pInt)) == '*solid' .and. &
        IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'section' ) then

     elemSetName = ''
     materialName = ''

     do i = 3_pInt,myPos(1_pInt)
       if (IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'elset') /= '') &
         elemSetName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'elset'))
       if (IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'material') /= '') &
         materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,i)),'material'))
     enddo

     if (elemSetName /= '' .and. materialName /= '') then
       c = c + 1_pInt
       mesh_nameMaterial(c) = materialName         ! name of material used for this section
       mesh_mapMaterial(c)  = elemSetName          ! mapped to respective element set
     endif       
   endif
 enddo

620 if (c==0_pInt) call IO_error(error_ID=905_pInt)
 do i=1_pInt,c
!   write(6,*)'name of materials: ',i,mesh_nameMaterial(i)
!   write(6,*)'name of elemSets:  ',i,mesh_mapMaterial(i)
   if (mesh_nameMaterial(i)=='' .or. mesh_mapMaterial(i)=='') call IO_error(error_ID=905_pInt)
 enddo

 end subroutine mesh_abaqus_map_materials



!********************************************************************
! map nodes from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPnode
!********************************************************************
subroutine mesh_spectral_map_nodes

 implicit none
 integer(pInt) :: i

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

 forall (i = 1_pInt:mesh_Nnodes) &
   mesh_mapFEtoCPnode(1:2,i) = i
 
end subroutine mesh_spectral_map_nodes



!********************************************************************
! map nodes from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPnode
!********************************************************************
subroutine mesh_marc_map_nodes(myUnit)

 use math, only: qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_fixedIntValue
                 
 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt), dimension (mesh_Nnodes) :: node_count
 integer(pInt) :: i

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

610 FORMAT(A300)

 node_count = 0_pInt

 rewind(myUnit)
 do
   read (myUnit,610,END=650) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'coordinates' ) then
     read (myUnit,610,END=650) line                                         ! skip crap line
     do i = 1_pInt,mesh_Nnodes
       read (myUnit,610,END=650) line
       mesh_mapFEtoCPnode(1_pInt,i) = IO_fixedIntValue (line,[ 0_pInt,10_pInt],1_pInt)
       mesh_mapFEtoCPnode(2_pInt,i) = i
     enddo
     exit
   endif
 enddo

650 call qsort(mesh_mapFEtoCPnode,1_pInt,int(size(mesh_mapFEtoCPnode,2_pInt),pInt))
 
end subroutine mesh_marc_map_nodes



!********************************************************************
! map nodes from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPnode
!********************************************************************
subroutine mesh_abaqus_map_nodes(myUnit)

 use math, only: qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_countDataLines, &
                 IO_intValue, &
                 IO_error

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt) :: i,c,cpNode = 0_pInt
 logical :: inPart = .false.

 allocate (mesh_mapFEtoCPnode(2_pInt,mesh_Nnodes)) ; mesh_mapFEtoCPnode = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=650) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
   ) then
     c = IO_countDataLines(myUnit)
     do i = 1_pInt,c
       backspace(myUnit)
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=650) line
       myPos = IO_stringPos(line,maxNchunks)
       cpNode = cpNode + 1_pInt
       mesh_mapFEtoCPnode(1_pInt,cpNode) = IO_intValue(line,myPos,1_pInt)
       mesh_mapFEtoCPnode(2_pInt,cpNode) = cpNode
     enddo
   endif
 enddo

650 call qsort(mesh_mapFEtoCPnode,1_pInt,int(size(mesh_mapFEtoCPnode,2_pInt),pInt))

 if (int(size(mesh_mapFEtoCPnode),pInt) == 0_pInt) call IO_error(error_ID=908_pInt)

end subroutine mesh_abaqus_map_nodes


!********************************************************************
! map elements from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPelem
!********************************************************************
subroutine mesh_spectral_map_elements

 implicit none
 integer(pInt) :: i

 allocate (mesh_mapFEtoCPelem(2_pInt,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

 forall (i = 1_pInt:mesh_NcpElems) &
   mesh_mapFEtoCPelem(1:2,i) = i

end subroutine mesh_spectral_map_elements



!********************************************************************
! map elements from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPelem
!********************************************************************
subroutine mesh_marc_map_elements(myUnit)

 use math, only: qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_continousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt), dimension (1_pInt+mesh_NcpElems) :: contInts
 integer(pInt) :: i,cpElem = 0_pInt

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=660) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'hypoelastic' ) then
     do i=1_pInt,3_pInt+hypoelasticTableStyle                                       ! skip three (or four if new table style!) lines
       read (myUnit,610,END=660) line 
     enddo
     contInts = IO_continousIntValues(myUnit,mesh_NcpElems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
     do i = 1_pInt,contInts(1)
       cpElem = cpElem+1_pInt
       mesh_mapFEtoCPelem(1,cpElem) = contInts(1_pInt+i)
       mesh_mapFEtoCPelem(2,cpElem) = cpElem
     enddo
   endif
 enddo

660 call qsort(mesh_mapFEtoCPelem,1_pInt,int(size(mesh_mapFEtoCPelem,2_pInt),pInt))           ! should be mesh_NcpElems

end subroutine mesh_marc_map_elements


!********************************************************************
! map elements from FE id to internal (consecutive) representation
!
! allocate globals: mesh_mapFEtoCPelem
!********************************************************************
subroutine mesh_abaqus_map_elements(myUnit)

 use math, only: qsort
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue, &
                 IO_error
                 
 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) ::i,j,k,cpElem = 0_pInt
 logical :: materialFound = .false.
 character (len=64) materialName,elemSetName ! why limited to 64? ABAQUS?

 allocate (mesh_mapFEtoCPelem(2,mesh_NcpElems)) ; mesh_mapFEtoCPelem = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do 
   read (myUnit,610,END=660) line
   myPos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_stringValue(line,myPos,1_pInt)) )
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'name'))     ! extract name=value
       materialFound = materialName /= ''                                           ! valid name?
     case('*user')
       if (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'material' .and. materialFound) then
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

660 call qsort(mesh_mapFEtoCPelem,1_pInt,int(size(mesh_mapFEtoCPelem,2_pInt),pInt))            ! should be mesh_NcpElems

 if (int(size(mesh_mapFEtoCPelem),pInt) < 2_pInt) call IO_error(error_ID=907_pInt)

end subroutine mesh_abaqus_map_elements


!********************************************************************
! get maximum count of nodes, IPs, IP neighbors, and subNodes
! among cpElements
!
! _maxNnodes, _maxNips, _maxNipNeighbors, _maxNsubNodes
!********************************************************************
subroutine mesh_spectral_count_cpSizes
 
 implicit none
 integer(pInt) :: t
 
 t = FE_mapElemtype('C3D8R')                   ! fake 3D hexahedral 8 node 1 IP element

 mesh_maxNnodes =       FE_Nnodes(t)
 mesh_maxNips =         FE_Nips(t)
 mesh_maxNipNeighbors = FE_NipNeighbors(t)
 mesh_maxNsubNodes =    FE_NsubNodes(t)

end subroutine mesh_spectral_count_cpSizes


!********************************************************************
! get maximum count of nodes, IPs, IP neighbors, and subNodes
! among cpElements
!
! _maxNnodes, _maxNips, _maxNipNeighbors, _maxNsubNodes
!********************************************************************
subroutine mesh_marc_count_cpSizes(myUnit)
 
 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_skipChunks

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,t,e

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNsubNodes    = 0_pInt
 
610 FORMAT(A300)
 rewind(myUnit)
 do
   read (myUnit,610,END=630) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'connectivity' ) then
     read (myUnit,610,END=630) line                                                                 ! Garbage line
     do i=1_pInt,mesh_Nelems                                                                        ! read all elements
       read (myUnit,610,END=630) line
       myPos = IO_stringPos(line,maxNchunks)                                                        ! limit to id and type
       e = mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt))
       if (e /= 0_pInt) then
         t = FE_mapElemtype(IO_stringValue(line,myPos,2_pInt))
         mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(t))
         mesh_maxNips =         max(mesh_maxNips,FE_Nips(t))
         mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(t))
         mesh_maxNsubNodes =    max(mesh_maxNsubNodes,FE_NsubNodes(t))
         call IO_skipChunks(myUnit,FE_NoriginalNodes(t)-(myPos(1_pInt)-2_pInt))        ! read on if FE_Nnodes exceeds node count present on current line
       endif
     enddo
     exit
   endif
 enddo
 
630 end subroutine mesh_marc_count_cpSizes


!********************************************************************
! get maximum count of nodes, IPs, IP neighbors, and subNodes
! among cpElements
!
! _maxNnodes, _maxNips, _maxNipNeighbors, _maxNsubNodes
!********************************************************************
subroutine mesh_abaqus_count_cpSizes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_extractValue ,&
                 IO_error, &
                 IO_countDataLines, &
                 IO_intValue

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,c,t
 logical :: inPart

 mesh_maxNnodes       = 0_pInt
 mesh_maxNips         = 0_pInt
 mesh_maxNipNeighbors = 0_pInt
 mesh_maxNsubNodes    = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do
   read (myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'type'))  ! remember elem type
     if (t==0_pInt) call IO_error(error_ID=910_pInt,ext_msg='mesh_abaqus_count_cpSizes')
     c = IO_countDataLines(myUnit)
     do i = 1_pInt,c
       backspace(myUnit)
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=620) line
       myPos = IO_stringPos(line,maxNchunks)                                  ! limit to 64 nodes max
       if (mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt)) /= 0_pInt) then                                                     ! disregard non CP elems
         mesh_maxNnodes =       max(mesh_maxNnodes,FE_Nnodes(t))
         mesh_maxNips =         max(mesh_maxNips,FE_Nips(t))
         mesh_maxNipNeighbors = max(mesh_maxNipNeighbors,FE_NipNeighbors(t))
         mesh_maxNsubNodes =    max(mesh_maxNsubNodes,FE_NsubNodes(t))
       endif
     enddo
   endif
 enddo
 
620 end subroutine mesh_abaqus_count_cpSizes


!********************************************************************
! store x,y,z coordinates of all nodes in mesh
!
! allocate globals:
! _node
!********************************************************************
subroutine mesh_spectral_build_nodes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_floatValue, &
                 IO_intValue

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 7_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 integer(pInt) :: a = 1_pInt, &
                  b = 1_pInt, &
                  c = 1_pInt, & 
                  headerLength = 0_pInt,i,j,n
 real(pReal) ::   x = 1.0_pReal, &
                  y = 1.0_pReal, &
                  z = 1.0_pReal
 logical ::  gotResolution = .false. ,gotDimension = .false.
 character(len=1024) :: line, keyword

 allocate ( mesh_node0 (3,mesh_Nnodes) ); mesh_node0 = 0.0_pReal
 allocate ( mesh_node  (3,mesh_Nnodes) ); mesh_node  = 0.0_pReal
 
 rewind(myUnit)
 read(myUnit,'(a1024)') line
 myPos = IO_stringPos(line,2_pInt)
 keyword = IO_lc(IO_StringValue(line,myPos,2_pInt))
 if (keyword(1:4) == 'head') then 
   headerLength = IO_intValue(line,myPos,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=842_pInt)
 endif
 
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   myPos = IO_stringPos(line,maxNchunks)             
   select case ( IO_lc(IO_StringValue(line,myPos,1_pInt)) )
     case ('dimension')
       gotDimension = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,myPos,j)))
           case('x')
              x = IO_floatValue(line,myPos,j+1_pInt)
           case('y')
              y = IO_floatValue(line,myPos,j+1_pInt)
           case('z')
              z = IO_floatValue(line,myPos,j+1_pInt)
         end select
       enddo
     case ('resolution')
       gotResolution = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,myPos,j)))
           case('a')
             a = 1_pInt + IO_intValue(line,myPos,j+1_pInt)
           case('b')
             b = 1_pInt + IO_intValue(line,myPos,j+1_pInt)
           case('c')
             c = 1_pInt + IO_intValue(line,myPos,j+1_pInt)
         end select
       enddo
   end select
 enddo

! --- sanity checks ---

 if ((.not. gotDimension) .or. (.not. gotResolution)) call IO_error(error_ID=842_pInt)
 if ((a < 1_pInt) .or. (b < 1_pInt) .or. (c < 0_pInt)) call IO_error(error_ID=843_pInt)           ! 1_pInt is already added
 if ((x <= 0.0_pReal) .or. (y <= 0.0_pReal) .or. (z <= 0.0_pReal)) call IO_error(error_ID=844_pInt)
 
 forall (n = 0_pInt:mesh_Nnodes-1_pInt)
   mesh_node0(1,n+1_pInt) = x * real(mod(n,a),pReal) / real(a-1_pInt,pReal)
   mesh_node0(2,n+1_pInt) = y * real(mod(n/a,b),pReal) / real(b-1_pInt,pReal)
   mesh_node0(3,n+1_pInt) = z * real(mod(n/a/b,c),pReal) / real(c-1_pInt,pReal)
 end forall 

 mesh_node = mesh_node0                                         !why?

end subroutine mesh_spectral_build_nodes


!********************************************************************
! store x,y,z coordinates of all nodes in mesh
!
! allocate globals:
! _node
!********************************************************************
subroutine mesh_marc_build_nodes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_stringPos, &
                 IO_fixedIntValue, &
                 IO_fixedNoEFloatValue

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), dimension(5), parameter :: node_ends = int([0,10,30,50,70],pInt)
 integer(pInt), parameter :: maxNchunks = 1_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,j,m

 allocate ( mesh_node0 (3,mesh_Nnodes) ); mesh_node0 = 0.0_pReal
 allocate ( mesh_node  (3,mesh_Nnodes) ); mesh_node  = 0.0_pReal

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=670) line
   myPos = IO_stringPos(line,maxNchunks)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'coordinates' ) then
     read (myUnit,610,END=670) line                                         ! skip crap line
     do i=1_pInt,mesh_Nnodes
       read (myUnit,610,END=670) line
       m = mesh_FEasCP('node',IO_fixedIntValue(line,node_ends,1_pInt))
       forall (j = 1_pInt:3_pInt) mesh_node0(j,m) = IO_fixedNoEFloatValue(line,node_ends,j+1_pInt)
     enddo
     exit
   endif
 enddo

670 mesh_node = mesh_node0

end subroutine mesh_marc_build_nodes


!********************************************************************
! store x,y,z coordinates of all nodes in mesh
!
! allocate globals:
! _node
!********************************************************************
subroutine mesh_abaqus_build_nodes(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_floatValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_countDataLines, &
                 IO_intValue

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 4_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) :: line
 integer(pInt) :: i,j,m,c
 logical :: inPart
 
 allocate ( mesh_node0 (3,mesh_Nnodes) ); mesh_node0 = 0.0_pReal
 allocate ( mesh_node  (3,mesh_Nnodes) ); mesh_node  = 0.0_pReal

610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do
   read (myUnit,610,END=670) line
   myPos = IO_stringPos(line,maxNchunks)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*node' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'print'    .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'file'     .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
   ) then
     c = IO_countDataLines(myUnit)                                                                  ! how many nodes are defined here?
     do i = 1_pInt,c
       backspace(myUnit)                                                                            ! rewind to first entry
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=670) line
       myPos = IO_stringPos(line,maxNchunks)
       m = mesh_FEasCP('node',IO_intValue(line,myPos,1_pInt))
       forall (j=1_pInt:3_pInt) mesh_node0(j,m) = IO_floatValue(line,myPos,j+1_pInt)
     enddo
   endif
 enddo

670 if (int(size(mesh_node0,2_pInt),pInt) /= mesh_Nnodes) call IO_error(error_ID=909_pInt)
 mesh_node = mesh_node0

end subroutine mesh_abaqus_build_nodes


!********************************************************************
! store FEid, type, mat, tex, and node list per element
!
! allocate globals:
! _element
!********************************************************************
subroutine mesh_spectral_build_elements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_floatValue, &
                 IO_stringPos, &
                 IO_error, &
                 IO_continousIntValues, &
                 IO_intValue, &
                 IO_countContinousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 7_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 integer(pInt) :: a = 1_pInt, b = 1_pInt, c = 1_pInt
 integer(pInt) :: e, i, j, homog = 0_pInt, headerLength = 0_pInt, maxIntCount
 integer(pInt), dimension(:), allocatable :: microstructures
 integer(pInt), dimension(1,1) :: dummySet = 0_pInt
 character(len=65536) :: line,keyword
 character(len=64), dimension(1) :: dummyName = ''

 rewind(myUnit)
 read(myUnit,'(a65536)') line
 myPos = IO_stringPos(line,2_pInt)
 keyword = IO_lc(IO_StringValue(line,myPos,2_pInt))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,myPos,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=842_pInt)
 endif
 
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a65536)') line
   myPos = IO_stringPos(line,maxNchunks)             
   select case ( IO_lc(IO_StringValue(line,myPos,1_pInt)) )
     case ('resolution')
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,myPos,j)))
           case('a')
             a = 1_pInt + IO_intValue(line,myPos,j+1_pInt)
           case('b')
             b = 1_pInt + IO_intValue(line,myPos,j+1_pInt)
           case('c')
             c = 1_pInt + IO_intValue(line,myPos,j+1_pInt)
         end select
       enddo
     case ('homogenization')
       homog = IO_intValue(line,myPos,2_pInt)
   end select
 enddo

 maxIntCount = 0_pInt
 i = 1_pInt

 do while (i > 0_pInt)
   i = IO_countContinousIntValues(myUnit)
   maxIntCount = max(maxIntCount, i)
 enddo

 rewind (myUnit)
 do i=1_pInt,headerLength                                         ! skip header
   read(myUnit,'(a65536)') line
 enddo

 allocate (mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt
 allocate (microstructures (1_pInt+maxIntCount))             ; microstructures = 2_pInt
 
 e = 0_pInt
 do while (e < mesh_NcpElems .and. microstructures(1) > 0_pInt)                                     ! fill expected number of elements, stop at end of data (or blank line!)
   microstructures = IO_continousIntValues(myUnit,maxIntCount,dummyName,dummySet,0_pInt)            ! get affected elements
   do i = 1_pInt,microstructures(1_pInt)
     e = e+1_pInt                                                                                   ! valid element entry
     mesh_element( 1,e) = e                                                                         ! FE id
     mesh_element( 2,e) = FE_mapElemtype('C3D8R')                                                   ! elem type
     mesh_element( 3,e) = homog                                                                     ! homogenization
     mesh_element( 4,e) = microstructures(1_pInt+i)                                                 ! microstructure
     mesh_element( 5,e) = e + (e-1_pInt)/(a-1_pInt) + ((e-1_pInt)/((a-1_pInt)*(b-1_pInt)))*a        ! base node
     mesh_element( 6,e) = mesh_element(5,e) + 1_pInt
     mesh_element( 7,e) = mesh_element(5,e) + a + 1_pInt
     mesh_element( 8,e) = mesh_element(5,e) + a
     mesh_element( 9,e) = mesh_element(5,e) + a * b                                                 ! second floor base node
     mesh_element(10,e) = mesh_element(9,e) + 1_pInt
     mesh_element(11,e) = mesh_element(9,e) + a + 1_pInt
     mesh_element(12,e) = mesh_element(9,e) + a
     mesh_maxValStateVar(1) = max(mesh_maxValStateVar(1),mesh_element(3,e))                         !needed for statistics
     mesh_maxValStateVar(2) = max(mesh_maxValStateVar(2),mesh_element(4,e))              
   enddo
 enddo

 deallocate(microstructures)
 if (e /= mesh_NcpElems) call IO_error(880_pInt,e)

end subroutine mesh_spectral_build_elements


!********************************************************************
! store FEid, type, mat, tex, and node list per element
!
! allocate globals:
! _element
!********************************************************************
subroutine mesh_marc_build_elements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_fixedNoEFloatValue, &
                 IO_skipChunks, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_continousIntValues

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 66_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos
 character(len=300) line

 integer(pInt), dimension(1_pInt+mesh_NcpElems) :: contInts
 integer(pInt) :: i,j,sv,myVal,e

 allocate (mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

610 FORMAT(A300)

 rewind(myUnit)
 do
   read (myUnit,610,END=620) line
   myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
   if( IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'connectivity' ) then
     read (myUnit,610,END=620) line                                                                 ! Garbage line
     do i = 1_pInt,mesh_Nelems
       read (myUnit,610,END=620) line
       myPos = IO_stringPos(line,maxNchunks)                                                        ! limit to 64 nodes max (plus ID, type)
       e = mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt))
       if (e /= 0_pInt) then                                                                        ! disregard non CP elems
         mesh_element(1,e) = IO_IntValue (line,myPos,1_pInt)                                        ! FE id
         mesh_element(2,e) = FE_mapElemtype(IO_StringValue(line,myPos,2_pInt))                      ! elem type
           forall (j = 1_pInt:FE_Nnodes(mesh_element(2,e))) &
             mesh_element(j+4_pInt,e) = IO_IntValue(line,myPos,j+2_pInt)                            ! copy FE ids of nodes
           call IO_skipChunks(myUnit,FE_NoriginalNodes(mesh_element(2_pInt,e))-(myPos(1_pInt)-2_pInt))        ! read on if FE_Nnodes exceeds node count present on current line
       endif
     enddo
     exit
   endif
 enddo
 
620 rewind(myUnit)                                                                                  ! just in case "initial state" apears before "connectivity"
 read (myUnit,610,END=620) line
 do
   myPos(1:1+2*2) = IO_stringPos(line,2_pInt)
   if( (IO_lc(IO_stringValue(line,myPos,1_pInt)) == 'initial') .and. &
       (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'state') ) then
     if (initialcondTableStyle == 2_pInt) read (myUnit,610,END=620) line                            ! read extra line for new style     
     read (myUnit,610,END=630) line                                                                 ! read line with index of state var
     myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
     sv = IO_IntValue(line,myPos,1_pInt)                                                            ! figure state variable index
     if( (sv == 2_pInt).or.(sv == 3_pInt) ) then                                                    ! only state vars 2 and 3 of interest
       read (myUnit,610,END=620) line                                                               ! read line with value of state var
       myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
       do while (scan(IO_stringValue(line,myPos,1_pInt),'+-',back=.true.)>1)                        ! is noEfloat value?
         myVal = nint(IO_fixedNoEFloatValue(line,[0_pInt,20_pInt],1_pInt),pInt)                     ! state var's value
         mesh_maxValStateVar(sv-1_pInt) = max(myVal,mesh_maxValStateVar(sv-1_pInt))                 ! remember max val of homogenization and microstructure index
         if (initialcondTableStyle == 2_pInt) then
           read (myUnit,610,END=630) line                                                           ! read extra line     
           read (myUnit,610,END=630) line                                                           ! read extra line     
         endif
         contInts = IO_continousIntValues&                                                          ! get affected elements
                   (myUnit,mesh_Nelems,mesh_nameElemSet,mesh_mapElemSet,mesh_NelemSets)
         do i = 1_pInt,contInts(1)
           e = mesh_FEasCP('elem',contInts(1_pInt+i))
           mesh_element(1_pInt+sv,e) = myVal
         enddo
         if (initialcondTableStyle == 0_pInt) read (myUnit,610,END=620) line                        ! ignore IP range for old table style
         read (myUnit,610,END=630) line
         myPos(1:1+2*1) = IO_stringPos(line,1_pInt)
       enddo
     endif
   else   
     read (myUnit,610,END=630) line
   endif
 enddo

630 end subroutine mesh_marc_build_elements

!********************************************************************
! store FEid, type, mat, tex, and node list per element
!
! allocate globals:
! _element
!********************************************************************
subroutine mesh_abaqus_build_elements(myUnit)

 use IO,   only: IO_lc, &
                 IO_stringValue, &
                 IO_skipChunks, &
                 IO_stringPos, &
                 IO_intValue, &
                 IO_extractValue, &
                 IO_floatValue, &
                 IO_error, &
                 IO_countDataLines

 implicit none
 integer(pInt), intent(in) :: myUnit

 integer(pInt), parameter :: maxNchunks = 65_pInt
 integer(pInt), dimension (1_pInt+2_pInt*maxNchunks) :: myPos

 integer(pInt) :: i,j,k,c,e,t,homog,micro
 logical inPart,materialFound
 character (len=64) :: materialName,elemSetName
 character(len=300) :: line

 allocate (mesh_element (4_pInt+mesh_maxNnodes,mesh_NcpElems)) ; mesh_element = 0_pInt

610 FORMAT(A300)

 inPart = .false.
 rewind(myUnit)
 do
   read (myUnit,610,END=620) line
   myPos(1:1+2*2) = IO_stringPos(line,2_pInt)
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) inPart = .true.
   if ( IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*end' .and. &
        IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'part' ) inPart = .false.

   if( (inPart .or. noPart) .and. &
       IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*element' .and. &
       ( IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'output'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'matrix'   .and. &
         IO_lc(IO_stringValue(line,myPos,2_pInt)) /= 'response' ) &
     ) then
     t = FE_mapElemtype(IO_extractValue(IO_lc(IO_stringValue(line,myPos,2_pInt)),'type'))  ! remember elem type
     if (t==0_pInt) call IO_error(error_ID=910_pInt,ext_msg='mesh_abaqus_build_elements')
     c = IO_countDataLines(myUnit)
     do i = 1_pInt,c
       backspace(myUnit)
     enddo
     do i = 1_pInt,c
       read (myUnit,610,END=620) line
       myPos = IO_stringPos(line,maxNchunks)                                  ! limit to 64 nodes max
       e = mesh_FEasCP('elem',IO_intValue(line,myPos,1_pInt))
       if (e /= 0_pInt) then                                                     ! disregard non CP elems
         mesh_element(1,e) = IO_intValue(line,myPos,1_pInt)                        ! FE id
         mesh_element(2,e) = t                                              ! elem type
         forall (j=1_pInt:FE_Nnodes(t)) &
           mesh_element(4_pInt+j,e) = IO_intValue(line,myPos,1_pInt+j)                  ! copy FE ids of nodes to position 5:
         call IO_skipChunks(myUnit,FE_NoriginalNodes(t)-(myPos(1_pInt)-1_pInt))           ! read on (even multiple lines) if FE_NoriginalNodes exceeds required node count
       endif
     enddo
   endif
 enddo

 
620 rewind(myUnit)                                                                                  ! just in case "*material" definitions apear before "*element"

 materialFound = .false.
 do 
   read (myUnit,610,END=630) line
   myPos = IO_stringPos(line,maxNchunks)
   select case ( IO_lc(IO_StringValue(line,myPos,1_pInt)))
     case('*material')
       materialName = trim(IO_extractValue(IO_lc(IO_StringValue(line,myPos,2_pInt)),'name'))        ! extract name=value
       materialFound = materialName /= ''                                                           ! valid name?
     case('*user')
       if ( IO_lc(IO_StringValue(line,myPos,2_pInt)) == 'material' .and. &
            materialFound ) then
         read (myUnit,610,END=630) line                                                             ! read homogenization and microstructure
         myPos(1:1+2*2) = IO_stringPos(line,2_pInt)
         homog = nint(IO_floatValue(line,myPos,1_pInt),pInt)
         micro = nint(IO_floatValue(line,myPos,2_pInt),pInt)
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


!********************************************************************
! get maximum count of shared elements among cpElements and
! build list of elements shared by each node in mesh
!
! _maxNsharedElems
! _sharedElem
!********************************************************************
subroutine mesh_build_sharedElems

implicit none
integer(pint)   e, &                                                                                ! element index
                t, &                                                                                ! element type
                node, &                                                                             ! CP node index
                j, &                                                                                ! node index per element 
                myDim, &                                                                            ! dimension index 
                nodeTwin                                                                            ! node twin in the specified dimension
integer(pInt), dimension (mesh_Nnodes) :: node_count
integer(pInt), dimension (:), allocatable :: node_seen

allocate(node_seen(maxval(FE_Nnodes)))


node_count = 0_pInt

do e = 1_pInt,mesh_NcpElems
  t = mesh_element(2,e)                                                                             ! get element type

  node_seen = 0_pInt                                                                                ! reset node duplicates
  do j = 1_pInt,FE_Nnodes(t)                                                                             ! check each node of element
    node = mesh_FEasCP('node',mesh_element(4+j,e))                                                  ! translate to internal (consecutive) numbering
    if (all(node_seen /= node)) then
      node_count(node) = node_count(node) + 1_pInt                                                  ! if FE node not yet encountered -> count it
      do myDim = 1_pInt,3_pInt                                                                                  ! check in each dimension...
        nodeTwin = mesh_nodeTwins(myDim,node)
        if (nodeTwin > 0_pInt) &                                                                    ! if I am a twin of some node...
          node_count(nodeTwin) = node_count(nodeTwin) + 1_pInt                                      ! -> count me again for the twin node
      enddo
    endif
    node_seen(j) = node                                                                             ! remember this node to be counted already
  enddo
enddo

mesh_maxNsharedElems = int(maxval(node_count),pInt)                                                 ! most shared node

allocate(mesh_sharedElem(1+mesh_maxNsharedElems,mesh_Nnodes))
mesh_sharedElem = 0_pInt

do e = 1_pInt,mesh_NcpElems
  t = mesh_element(2,e)
  node_seen = 0_pInt
  do j = 1_pInt,FE_Nnodes(t)
    node = mesh_FEasCP('node',mesh_element(4_pInt+j,e))
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
    node_seen(j) = node
  enddo
enddo

deallocate(node_seen)

end subroutine mesh_build_sharedElems


!***********************************************************
! build up of IP neighborhood
!
! allocate globals
! _ipNeighborhood
!***********************************************************
subroutine mesh_build_ipNeighborhood

implicit none
integer(pInt)                   myElem, &                                                           ! my CP element index
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
                                a, anchor
integer(pInt), dimension(FE_maxmaxNnodesAtIP) :: &
                                linkedNodes = 0_pInt, &
                                matchingNodes
logical checkTwins

allocate(mesh_ipNeighborhood(2,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems))
mesh_ipNeighborhood = 0_pInt


do myElem = 1_pInt,mesh_NcpElems                                                                    ! loop over cpElems
  myType = mesh_element(2,myElem)                                                                   ! get elemType
  do myIP = 1_pInt,FE_Nips(myType)                                                                  ! loop over IPs of elem

    do neighbor = 1_pInt,FE_NipNeighbors(myType)                                                    ! loop over neighbors of IP
      neighboringIPkey = FE_ipNeighbor(neighbor,myIP,myType)

      !*** if the key is positive, the neighbor is inside the element
      !*** that means, we have already found our neighboring IP
      
      if (neighboringIPkey > 0_pInt) then
        mesh_ipNeighborhood(1,neighbor,myIP,myElem) = myElem
        mesh_ipNeighborhood(2,neighbor,myIP,myElem) = neighboringIPkey


      !*** if the key is negative, the neighbor resides in a neighboring element
      !*** that means, we have to look through the face indicated by the key and see which element is behind that face
      
      elseif (neighboringIPkey < 0_pInt) then                                            ! neighboring element's IP
        myFace = -neighboringIPkey
        call mesh_faceMatch(myElem, myFace, matchingElem, matchingFace)             ! get face and CP elem id of face match
        if (matchingElem > 0_pInt) then                                             ! found match?
          neighboringType = mesh_element(2,matchingElem)

          !*** trivial solution if neighbor has only one IP
          
          if (FE_Nips(neighboringType) == 1_pInt) then            
            mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
            mesh_ipNeighborhood(2,neighbor,myIP,myElem) = 1_pInt
            cycle
          endif

          !*** find those nodes which build the link to the neighbor
          
          NlinkedNodes = 0_pInt
          linkedNodes = 0_pInt
          do a = 1_pInt,FE_maxNnodesAtIP(myType)                                         ! figure my anchor nodes on connecting face
            anchor = FE_nodesAtIP(a,myIP,myType)
            if (anchor /= 0_pInt) then                                              ! valid anchor node
              if (any(FE_nodeOnFace(:,myFace,myType) == anchor)) then               ! ip anchor sits on face?
                NlinkedNodes = NlinkedNodes + 1_pInt
                linkedNodes(NlinkedNodes) = &
                   mesh_FEasCP('node',mesh_element(4_pInt+anchor,myElem))                ! CP id of anchor node
              else                                                                  ! something went wrong with the linkage, since not all anchors sit on my face
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
            do a = 1_pInt,FE_maxNnodesAtIP(neighboringType)                              ! check each anchor node of that ip
              anchor = FE_nodesAtIP(a,candidateIP,neighboringType)
              if (anchor /= 0_pInt) then                                            ! valid anchor node
                if (any(FE_nodeOnFace(:,matchingFace,neighboringType) == anchor)) then  ! sits on matching face?
                  NmatchingNodes = NmatchingNodes + 1_pInt
                  matchingNodes(NmatchingNodes) = &
                     mesh_FEasCP('node',mesh_element(4+anchor,matchingElem))        ! CP id of neighbor's anchor node
                else                                                                ! no matching, because not all nodes sit on the matching face
                  NmatchingNodes = 0_pInt
                  matchingNodes = 0_pInt
                  exit
                endif
              endif
            enddo

            if (NmatchingNodes /= NlinkedNodes) &                                   ! this ip has wrong count of anchors on face
              cycle checkCandidateIP
            
            !*** check "normal" nodes whether they match or not
            
            checkTwins = .false.
            do a = 1_pInt,NlinkedNodes
              if (all(matchingNodes /= linkedNodes(a))) then                        ! this linkedNode does not match any matchingNode
                checkTwins = .true.
                exit                                                                ! no need to search further
              endif
            enddo
            
            !*** if no match found, then also check node twins
            
            if(checkTwins) then
              dir = int(maxloc(abs(mesh_ipAreaNormal(1:3,neighbor,myIP,myElem)),1),pInt)      ! check for twins only in direction of the surface normal
              do a = 1_pInt,NlinkedNodes
                twin_of_linkedNode = mesh_nodeTwins(dir,linkedNodes(a))
                if (twin_of_linkedNode == 0_pInt .or. &                             ! twin of linkedNode does not exist...
                    all(matchingNodes /= twin_of_linkedNode)) then                  ! ... or it does not match any matchingNode
                  cycle checkCandidateIP                                            ! ... then check next candidateIP
                endif
              enddo
            endif

            !*** we found a match !!!

            mesh_ipNeighborhood(1,neighbor,myIP,myElem) = matchingElem
            mesh_ipNeighborhood(2,neighbor,myIP,myElem) = candidateIP
            exit checkCandidateIP            
          enddo checkCandidateIP
        endif                                                                       ! end of valid external matching
      endif                                                                         ! end of internal/external matching
    enddo
  enddo
enddo

end subroutine mesh_build_ipNeighborhood



!***********************************************************
! assignment of coordinates for subnodes in each cp element
!
! allocate globals
! _subNodeCoord
!***********************************************************
subroutine mesh_build_subNodeCoords
 
 implicit none
 integer(pInt) e,t,n,p
 
 if (.not. allocated(mesh_subNodeCoord)) then
   allocate(mesh_subNodeCoord(3,mesh_maxNnodes+mesh_maxNsubNodes,mesh_NcpElems))
 endif
 mesh_subNodeCoord = 0.0_pReal
 
 do e = 1_pInt,mesh_NcpElems                   ! loop over cpElems
   t = mesh_element(2,e)                  ! get elemType
   do n = 1_pInt,FE_Nnodes(t)
     mesh_subNodeCoord(1:3,n,e) = mesh_node(1:3,mesh_FEasCP('node',mesh_element(4_pInt+n,e))) ! loop over nodes of this element type
   enddo
   do n = 1_pInt,FE_NsubNodes(t)               ! now for the true subnodes
     do p = 1_pInt,FE_Nips(t)                  ! loop through possible parent nodes
       if (FE_subNodeParent(p,n,t) > 0_pInt) & ! valid parent node
         mesh_subNodeCoord(1:3,FE_Nnodes(t)+n,e) = mesh_subNodeCoord(1:3,FE_Nnodes(t)+n,e) &
                                                 + mesh_node(1:3,mesh_FEasCP('node',mesh_element(4_pInt+FE_subNodeParent(p,n,t),e))) ! add up parents
     enddo
     mesh_subNodeCoord(1:3,n+FE_Nnodes(t),e) = mesh_subNodeCoord(1:3,n+FE_Nnodes(t),e) &
                                             /real(count(FE_subNodeParent(:,n,t) > 0_pInt),pReal)
   enddo
 enddo 
 
end subroutine mesh_build_subNodeCoords


!***********************************************************
! calculation of IP coordinates
!
! allocate globals
! _ipCenterOfGravity
!***********************************************************
subroutine mesh_build_ipCoordinates
 
 use prec, only: tol_gravityNodePos
 
 implicit none
 integer(pInt) :: e,f,t,i,j,k,n
 logical, dimension(mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNode            ! flagList to find subnodes determining center of grav
 real(pReal), dimension(3,mesh_maxNnodes+mesh_maxNsubNodes) :: gravityNodePos   ! coordinates of subnodes determining center of grav
 real(pReal), dimension(3) :: centerOfGravity

 if (.not. allocated(mesh_ipCenterOfGravity)) then
   allocate(mesh_ipCenterOfGravity(3,mesh_maxNips,mesh_NcpElems))
 endif
 
 do e = 1_pInt,mesh_NcpElems                                    ! loop over cpElems
   t = mesh_element(2,e)                                   ! get elemType
   do i = 1_pInt,FE_Nips(t)                                     ! loop over IPs of elem
     gravityNode = .false.                                 ! reset flagList
     gravityNodePos = 0.0_pReal                            ! reset coordinates
     do f = 1_pInt,FE_NipNeighbors(t)                           ! loop over interfaces of IP
       do n = 1_pInt,FE_NipFaceNodes                            ! loop over nodes on interface
         gravityNode(FE_subNodeOnIPFace(n,f,i,t)) = .true.
         gravityNodePos(:,FE_subNodeOnIPFace(n,f,i,t)) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       enddo
     enddo
     
     do j = 1_pInt,mesh_maxNnodes+mesh_maxNsubNodes-1_pInt           ! walk through entire flagList except last
       if (gravityNode(j)) then                            ! valid node index
         do k = j+1_pInt,mesh_maxNnodes+mesh_maxNsubNodes       ! walk through remainder of list
           if (gravityNode(k) .and. all(abs(gravityNodePos(:,j) - gravityNodePos(:,k)) < tol_gravityNodePos)) then   ! found duplicate
             gravityNode(j) = .false.                      ! delete first instance
             gravityNodePos(:,j) = 0.0_pReal
             exit                                          ! continue with next suspect
           endif
         enddo
       endif
     enddo
     centerOfGravity = sum(gravityNodePos,2)/real(count(gravityNode),pReal)
     mesh_ipCenterOfGravity(:,i,e) = centerOfGravity
   enddo
 enddo

end subroutine mesh_build_ipCoordinates


!***********************************************************
! calculation of IP volume
!
! allocate globals
! _ipVolume
!***********************************************************
subroutine mesh_build_ipVolumes
 
 use math, only: math_volTetrahedron
 implicit none
 
 integer(pInt) :: e,f,t,i,j,n
 integer(pInt), parameter :: Ntriangles = FE_NipFaceNodes-2_pInt                ! each interface is made up of this many triangles
 real(pReal), dimension(3,FE_NipFaceNodes) :: nPos                              ! coordinates of nodes on IP face
 real(pReal), dimension(Ntriangles,FE_NipFaceNodes) :: volume                   ! volumes of possible tetrahedra

 if (.not. allocated(mesh_ipVolume)) then
   allocate(mesh_ipVolume(mesh_maxNips,mesh_NcpElems))
 endif
 
 mesh_ipVolume = 0.0_pReal 
 do e = 1_pInt,mesh_NcpElems                                    ! loop over cpElems
   t = mesh_element(2_pInt,e)                                   ! get elemType
   do i = 1_pInt,FE_Nips(t)                                     ! loop over IPs of elem
     do f = 1_pInt,FE_NipNeighbors(t)         ! loop over interfaces of IP and add tetrahedra which connect to CoG
       forall (n = 1_pInt:FE_NipFaceNodes) &
         nPos(:,n) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       forall (n = 1_pInt:FE_NipFaceNodes, j = 1_pInt:Ntriangles) &  ! start at each interface node and build valid triangles to cover interface
         volume(j,n) = math_volTetrahedron(nPos(:,n), &    ! calc volume of respective tetrahedron to CoG
                                           nPos(:,1_pInt+mod(n-1_pInt +j  ,FE_NipFaceNodes)), & ! start at offset j
                                           nPos(:,1_pInt+mod(n-1_pInt +j+1_pInt,FE_NipFaceNodes)), & ! and take j's neighbor
                                           mesh_ipCenterOfGravity(:,i,e))
       mesh_ipVolume(i,e) = mesh_ipVolume(i,e) + sum(volume)    ! add contribution from this interface
     enddo
     mesh_ipVolume(i,e) = mesh_ipVolume(i,e) / FE_NipFaceNodes  ! renormalize with interfaceNodeNum due to loop over them
   enddo
 enddo

end subroutine mesh_build_ipVolumes


!***********************************************************
! calculation of IP interface areas
!
! allocate globals
! _ipArea, _ipAreaNormal
!***********************************************************
subroutine mesh_build_ipAreas
 
 use math, only: math_vectorproduct
 
 implicit none
 integer(pInt) :: e,f,t,i,j,n
 integer(pInt), parameter :: Ntriangles = FE_NipFaceNodes-2_pInt     ! each interface is made up of this many triangles
 real(pReal), dimension (3,FE_NipFaceNodes) :: nPos             ! coordinates of nodes on IP face
 real(pReal), dimension(3,Ntriangles,FE_NipFaceNodes) :: normal
 real(pReal), dimension(Ntriangles,FE_NipFaceNodes)   :: area

 allocate(mesh_ipArea(mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ;         mesh_ipArea       = 0.0_pReal
 allocate(mesh_ipAreaNormal(3_pInt,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems)) ; mesh_ipAreaNormal = 0.0_pReal
 do e = 1_pInt,mesh_NcpElems                                         ! loop over cpElems
   t = mesh_element(2,e)                                        ! get elemType
   do i = 1_pInt,FE_Nips(t)                                          ! loop over IPs of elem
     do f = 1_pInt,FE_NipNeighbors(t)                                ! loop over interfaces of IP 
       forall (n = 1_pInt:FE_NipFaceNodes) nPos(:,n) = mesh_subNodeCoord(:,FE_subNodeOnIPFace(n,f,i,t),e)
       forall (n = 1_pInt:FE_NipFaceNodes, j = 1_pInt:Ntriangles)         ! start at each interface node and build valid triangles to cover interface
         normal(:,j,n) = math_vectorproduct(nPos(:,1_pInt+mod(n+j-1_pInt,FE_NipFaceNodes)) - nPos(:,n), &    ! calc their normal vectors
                                            nPos(:,1_pInt+mod(n+j-0_pInt,FE_NipFaceNodes)) - nPos(:,n))
         area(j,n) = sqrt(sum(normal(:,j,n)*normal(:,j,n)))                                       ! and area
       end forall
       forall (n = 1_pInt:FE_NipFaceNodes, j = 1_pInt:Ntriangles, area(j,n) > 0.0_pReal) &
         normal(1:3,j,n) = normal(1:3,j,n) / area(j,n)            ! make myUnit normal
       
       mesh_ipArea(f,i,e) = sum(area) / (FE_NipFaceNodes*2.0_pReal)                   ! area of parallelograms instead of triangles
       mesh_ipAreaNormal(:,f,i,e) = sum(sum(normal,3),2_pInt)/&                            ! average of all valid normals
                                        real(count(area > 0.0_pReal),pReal)
     enddo
   enddo
 enddo
 
 end subroutine mesh_build_ipAreas


!***********************************************************
! assignment of twin nodes for each cp node
!
! allocate globals
! _nodeTwins
!***********************************************************
subroutine mesh_build_nodeTwins

implicit none
integer(pInt) dir, &      ! direction of periodicity
              node, &
              minimumNode, &
              maximumNode, &
              n1, &
              n2
integer(pInt), dimension(mesh_Nnodes+1) :: minimumNodes, maximumNodes ! list of surface nodes (minimum and maximum coordinate value) with first entry giving the number of nodes
real(pReal)   minCoord, maxCoord, &     ! extreme positions in one dimension
              tolerance                 ! tolerance below which positions are assumed identical
real(pReal), dimension(3) ::  distance  ! distance between two nodes in all three coordinates
logical, dimension(mesh_Nnodes) :: unpaired

allocate(mesh_nodeTwins(3,mesh_Nnodes))
mesh_nodeTwins = 0_pInt

tolerance = 0.001_pReal * minval(mesh_ipVolume) ** 0.333_pReal

do dir = 1_pInt,3_pInt                                    ! check periodicity in directions of x,y,z
  if (mesh_periodicSurface(dir)) then           ! only if periodicity is requested

    
    !*** find out which nodes sit on the surface 
    !*** and have a minimum or maximum position in this dimension
    
    minimumNodes = 0_pInt
    maximumNodes = 0_pInt
    minCoord = minval(mesh_node0(dir,:))
    maxCoord = maxval(mesh_node0(dir,:))
    do node = 1_pInt,mesh_Nnodes                     ! loop through all nodes and find surface nodes
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
          if (sum(distance) - distance(dir) <= tolerance) then        ! minimum possible distance (within tolerance)
            mesh_nodeTwins(dir,minimumNode) = maximumNode
            mesh_nodeTwins(dir,maximumNode) = minimumNode
            unpaired(maximumNode) = .false.                           ! remember this node, we don't have to look for his partner again
            exit
          endif
        enddo
      endif
    enddo

  endif
enddo

end subroutine mesh_build_nodeTwins



!***********************************************************
! write statistics regarding input file parsing
! to the output file
! 
!***********************************************************
subroutine mesh_tell_statistics

 use math,  only: math_range
 use IO,    only: IO_error
 use debug, only: debug_what, &
                  debug_mesh, &
                  debug_levelBasic, &
                  debug_levelExtensive, &
                  debug_levelSelective, &
                  debug_e, &
                  debug_i

 implicit none
 integer(pInt), dimension (:,:), allocatable :: mesh_HomogMicro
 character(len=64) :: myFmt
 integer(pInt) :: i,e,n,f,t, myDebug
 
 myDebug = debug_what(debug_mesh)

 if (mesh_maxValStateVar(1) < 1_pInt) call IO_error(error_ID=170_pInt) ! no homogenization specified
 if (mesh_maxValStateVar(2) < 1_pInt) call IO_error(error_ID=180_pInt) ! no microstructure specified
 
 allocate (mesh_HomogMicro(mesh_maxValStateVar(1),mesh_maxValStateVar(2))); mesh_HomogMicro = 0_pInt
do e = 1_pInt,mesh_NcpElems
  if (mesh_element(3,e) < 1_pInt) call IO_error(error_ID=170_pInt,e=e) ! no homogenization specified
  if (mesh_element(4,e) < 1_pInt) call IO_error(error_ID=180_pInt,e=e) ! no microstructure specified
  mesh_HomogMicro(mesh_element(3,e),mesh_element(4,e)) = &
  mesh_HomogMicro(mesh_element(3,e),mesh_element(4,e)) + 1_pInt ! count combinations of homogenization and microstructure
enddo
!$OMP CRITICAL (write2out)
  if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
    write (6,*)
    write (6,*) 'Input Parser: STATISTICS'
    write (6,*)
    write (6,*) mesh_Nelems,           ' : total number of elements in mesh'
    write (6,*) mesh_NcpElems,         ' : total number of CP elements in mesh'
    write (6,*) mesh_Nnodes,           ' : total number of nodes in mesh'
    write (6,*) mesh_maxNnodes,        ' : max number of nodes in any CP element'
    write (6,*) mesh_maxNips,          ' : max number of IPs in any CP element'
    write (6,*) mesh_maxNipNeighbors,  ' : max number of IP neighbors in any CP element'
    write (6,*) mesh_maxNsubNodes,     ' : max number of (additional) subnodes in any CP element'
    write (6,*) mesh_maxNsharedElems,  ' : max number of CP elements sharing a node'
    write (6,*)
    write (6,*) 'Input Parser: HOMOGENIZATION/MICROSTRUCTURE'
    write (6,*)
    write (6,*) mesh_maxValStateVar(1), ' : maximum homogenization index'
    write (6,*) mesh_maxValStateVar(2), ' : maximum microstructure index'
    write (6,*)
    write (myFmt,'(a,i32.32,a)') '(9x,a2,1x,',mesh_maxValStateVar(2),'(i8))'
    write (6,myFmt) '+-',math_range(mesh_maxValStateVar(2))
    write (myFmt,'(a,i32.32,a)') '(i8,1x,a2,1x,',mesh_maxValStateVar(2),'(i8))'
    do i=1_pInt,mesh_maxValStateVar(1)      ! loop over all (possibly assigned) homogenizations
      write (6,myFmt) i,'| ',mesh_HomogMicro(i,:) ! loop over all (possibly assigned) microstructures
    enddo
    write(6,*)
    write(6,*) 'Input Parser: ADDITIONAL MPIE OPTIONS'
    write(6,*)
    write(6,*) 'periodic surface : ', mesh_periodicSurface
    write(6,*)
    call flush(6)
  endif

  if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
    write (6,*)
    write (6,*) 'Input Parser: SUBNODE COORDINATES'
    write (6,*)
    write(6,'(a8,1x,a5,1x,2(a15,1x),a20,3(1x,a12))')&
                              'elem','IP','IP neighbor','IPFaceNodes','subNodeOnIPFace','x','y','z'
    do e = 1_pInt,mesh_NcpElems                  ! loop over cpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = mesh_element(2,e)                 ! get elemType
      do i = 1_pInt,FE_Nips(t)                   ! loop over IPs of elem
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        do f = 1_pInt,FE_NipNeighbors(t)         ! loop over interfaces of IP
          do n = 1_pInt,FE_NipFaceNodes          ! loop over nodes on interface
            write(6,'(i8,1x,i5,2(1x,i15),1x,i20,3(1x,f12.8))') e,i,f,n,FE_subNodeOnIPFace(n,f,i,t),&
                                               mesh_subNodeCoord(1,FE_subNodeOnIPFace(n,f,i,t),e),&
                                               mesh_subNodeCoord(2,FE_subNodeOnIPFace(n,f,i,t),e),&
                                               mesh_subNodeCoord(3,FE_subNodeOnIPFace(n,f,i,t),e)
          enddo
        enddo
      enddo
    enddo
    write(6,*)
    write(6,*) 'Input Parser: IP COORDINATES'
    write(6,'(a8,1x,a5,3(1x,a12))') 'elem','IP','x','y','z'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      do i = 1_pInt,FE_Nips(mesh_element(2,e))
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write (6,'(i8,1x,i5,3(1x,f12.8))') e, i, mesh_ipCenterOfGravity(:,i,e)
      enddo
    enddo 
    write (6,*)
    write (6,*) 'Input Parser: ELEMENT VOLUME'
    write (6,*)
    write (6,'(a13,1x,e15.8)') 'total volume', sum(mesh_ipVolume)
    write (6,*)
    write (6,'(a8,1x,a5,1x,a15,1x,a5,1x,a15,1x,a16)') 'elem','IP','volume','face','area','-- normal --'
    do e = 1_pInt,mesh_NcpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      do i = 1_pInt,FE_Nips(mesh_element(2,e))
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        write (6,'(i8,1x,i5,1x,e15.8)') e,i,mesh_IPvolume(i,e)
        do f = 1_pInt,FE_NipNeighbors(mesh_element(2,e))
          write (6,'(i33,1x,e15.8,1x,3(f6.3,1x))') f,mesh_ipArea(f,i,e),mesh_ipAreaNormal(:,f,i,e)
        enddo
      enddo
    enddo
    write (6,*)
    write (6,*) 'Input Parser: NODE TWINS'
    write (6,*)
    write(6,'(a6,3(3x,a6))') '  node','twin_x','twin_y','twin_z'
    do n = 1_pInt,mesh_Nnodes                    ! loop over cpNodes
      if (debug_e <= mesh_NcpElems) then
        if (any(mesh_element(5:,debug_e) == n)) then
          write(6,'(i6,3(3x,i6))') n, mesh_nodeTwins(1:3,n)
        endif
      endif
    enddo
    write(6,*)
    write(6,*) 'Input Parser: IP NEIGHBORHOOD'
    write(6,*)
    write(6,'(a8,1x,a10,1x,a10,1x,a3,1x,a13,1x,a13)') 'elem','IP','neighbor','','elemNeighbor','ipNeighbor'
    do e = 1_pInt,mesh_NcpElems                  ! loop over cpElems
      if (iand(myDebug,debug_levelSelective)   /= 0_pInt .and. debug_e /= e) cycle
      t = mesh_element(2,e)                 ! get elemType
      do i = 1_pInt,FE_Nips(t)                   ! loop over IPs of elem
        if (iand(myDebug,debug_levelSelective) /= 0_pInt .and. debug_i /= i) cycle
        do n = 1_pInt,FE_NipNeighbors(t)         ! loop over neighbors of IP
          write (6,'(i8,1x,i10,1x,i10,1x,a3,1x,i13,1x,i13)') e,i,n,'-->',mesh_ipNeighborhood(1,n,i,e),mesh_ipNeighborhood(2,n,i,e)
        enddo
      enddo
    enddo
  endif
!$OMP END CRITICAL (write2out)

deallocate(mesh_HomogMicro)
 
end subroutine mesh_tell_statistics

subroutine mesh_regrid(res,resNew)           !use new_res=0.0 for automatic determination of new grid
 use prec, only pInt, pReal
 use DAMASK_interface, only : getSolverJobName
 use IO, only : IO_read_jobBinaryFile

 integer(pInt), dimension(3), intent(in) :: res
 integer(pInt), dimension(3), intent(inout) :: resNew
 real(pReal), dimension(res(1),res(2),res(3),3,3) :: F

 real(pReal), dimension(:,:,:,:,:),      allocatable  :: crystallite_F0, &
                                                         CPFEM_dcsdE, &
                                                         crystallite_Fp0, &
                                                         crystallite_Lp0
 real(pReal), dimension (:,:,:,:,:,:,:), allocatable  :: crystallite_dPdF0
 real(pReal), dimension (:,:,:,:),       allocatable  :: crystallite_Tstar0_v, &
                                                         convergedStateConst 
 integer(pInt), dimension (:,:),         allocatable  :: convergedSizeConst 

 call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',trim(getSolverJobName()),size(F))
 read (777,rec=1) F
 close (777)
  
end subroutine mesh_regrid


end module mesh
