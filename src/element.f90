!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module element
 use prec, only: &
   pInt, &
   pReal

 implicit none
 private 

!---------------------------------------------------------------------------------------------------
!> Properties of a single element (the element used in the mesh)
!---------------------------------------------------------------------------------------------------
 type, public :: tElement
   integer(pInt) :: &
     elemType, &
     geomType, &                                                                                    ! geometry type (same for same dimension and same number of integration points)
     cellType, &
     Nnodes, &
     Ncellnodes, &
     NcellnodesPerCell, &
     nIPs, &
     nIPneighbors, &                                                                                   ! ToDo: MD: Do all IPs in one element type have the same number of neighbors?
     maxNnodeAtIP
   integer(pInt), dimension(:,:), allocatable :: &
     Cell, &                                                                                        ! intra-element (cell) nodes that constitute a cell
     NnodeAtIP, &
     IPneighbor, &
     cellFace
   real(pReal), dimension(:,:), allocatable :: &
 ! center of gravity of the weighted nodes gives the position of the cell node.
 ! example: face-centered cell node with face nodes 1,2,5,6 to be used in,
 !          e.g., an 8 node element, would be encoded:
 !          1, 1, 0, 0, 1, 1, 0, 0
     cellNodeParentNodeWeights
   contains
     procedure :: init => tElement_init
 end type

 integer(pInt), parameter, private :: &
   NELEMTYPE = 13_pInt

 integer(pInt), dimension(NelemType), parameter, private :: NNODE = &
 int([ &
      3, & ! 2D 3node 1ip
      6, & ! 2D 6node 3ip
      4, & ! 2D 4node 4ip
      8, & ! 2D 8node 9ip
      8, & ! 2D 8node 4ip
      !--------------------
      4, & ! 3D 4node 1ip
      5, & ! 3D 5node 4ip
     10, & ! 3D 10node 4ip
      6, & ! 3D 6node 6ip
      8, & ! 3D 8node 1ip
      8, & ! 3D 8node 8ip
     20, & ! 3D 20node 8ip
     20  & ! 3D 20node 27ip
  ],pInt)                                                                                           !< number of nodes that constitute a specific type of element

 integer(pInt), dimension(NelemType), parameter, public :: GEOMTYPE = &
 int([ &
      1, & ! 2D 3node 1ip
      2, & ! 2D 6node 3ip
      3, & ! 2D 4node 4ip
      4, & ! 2D 8node 9ip
      3, & ! 2D 8node 4ip
      !--------------------
      5, & ! 3D 4node 1ip
      6, & ! 3D 5node 4ip
      6, & ! 3D 10node 4ip
      7, & ! 3D 6node 6ip
      8, & ! 3D 8node 1ip
      9, & ! 3D 8node 8ip
      9, & ! 3D 20node 8ip
     10  & ! 3D 20node 27ip
  ],pInt)                                                                                           !< geometry type of particular element type

 !integer(pInt), dimension(maxval(geomType)), parameter, private :: NCELLNODE = &                   ! Intel 16.0 complains
 integer(pInt), dimension(10), parameter, private :: NCELLNODE = &
 int([ &
      3, &
      7, &
      9, &
     16, &
      4, &
     15, &
     21, &
      8, &
     27, &
     64  &
  ],pInt)                                                                                           !< number of cell nodes in a specific geometry type

 !integer(pInt), dimension(maxval(geomType)), parameter, private :: NIP = &                         ! Intel 16.0 complains
 integer(pInt), dimension(10), parameter, private :: NIP = &
 int([ &
      1, & 
      3, & 
      4, & 
      9, & 
      1, & 
      4, & 
      6, & 
      1, & 
      8, & 
     27  & 
  ],pInt)                                                                                           !< number of IPs in a specific geometry type

 !integer(pInt), dimension(maxval(geomType)), parameter, private  :: CELLTYPE = &                   ! Intel 16.0 complains
 integer(pInt), dimension(10), parameter, private  :: CELLTYPE = &                                  !< cell type that is used by each geometry type
 int([ &
      1, & ! 2D 3node
      2, & ! 2D 4node
      2, & ! 2D 4node
      2, & ! 2D 4node
      3, & ! 3D 4node
      4, & ! 3D 8node
      4, & ! 3D 8node
      4, & ! 3D 8node
      4, & ! 3D 8node
      4  & ! 3D 8node
  ],pInt)
 
 !integer(pInt), dimension(maxval(cellType)), parameter, private :: nIPNeighbor = &                 ! causes problem with Intel 16.0
 integer(pInt), dimension(4), parameter, private :: NIPNEIGHBOR = &                                 !< number of ip neighbors / cell faces in a specific cell type
 int([&
      3, & ! 2D 3node
      4, & ! 2D 4node
      4, & ! 3D 4node
      6  & ! 3D 8node
  ],pInt)

 !integer(pInt), dimension(maxval(cellType)), parameter, private :: NCELLNODESPERCELLFACE = &            
 integer(pInt), dimension(4), parameter, private :: NCELLNODEPERCELLFACE = &                        !< number of cell nodes in a specific cell type
 int([ &
      2, & ! 2D 3node
      2, & ! 2D 4node
      3, & ! 3D 4node
      4  & ! 3D 8node
  ],pInt)

 !integer(pInt), dimension(maxval(geomType)), parameter, private :: maxNodeAtIP = &                 ! causes problem with Intel 16.0
 integer(pInt), dimension(10), parameter, private :: maxNnodeAtIP = &                               !< maximum number of parent nodes that belong to an IP for a specific type of element
 int([ &
      3, &
      1, &
      1, &
      2, &
      4, &
      1, &
      1, &
      8, &
      1, &
      4  &
  ],pInt)


 !integer(pInt), dimension(maxval(CELLTYPE)), parameter, private  :: NCELLNODEPERCELL = &           ! Intel 16.0 complains
 integer(pInt), dimension(4), parameter, private :: NCELLNODEPERCELL = &                            !< number of cell nodes in a specific cell type
 int([ &
      3, & ! 2D 3node
      4, & ! 2D 4node
      4, & ! 3D 4node
      8  & ! 3D 8node
  ],pInt)

 integer(pInt), dimension(maxNnodeAtIP(1),nIP(1)), parameter, private :: NnodeAtIP1 = &
    reshape(int([&
    1,2,3   &
    ],pInt),[maxNnodeAtIP(1),nIP(1)])

 integer(pInt), dimension(maxNnodeAtIP(2),nIP(2)), parameter, private :: NnodeAtIP2 = &
    reshape(int([&
    1,  &
    2,  &
    3   &
    ],pInt),[maxNnodeAtIP(2),nIP(2)])

 integer(pInt), dimension(maxNnodeAtIP(3),nIP(3)), parameter, private :: NnodeAtIP3 = &
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3   &
    ],pInt),[maxNnodeAtIP(3),nIP(3)])

 integer(pInt), dimension(maxNnodeAtIP(4),nIP(4)), parameter, private :: NnodeAtIP4 = &
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
    ],pInt),[maxNnodeAtIP(4),nIP(4)])

 integer(pInt), dimension(maxNnodeAtIP(5),nIP(5)), parameter, private :: NnodeAtIP5 = &
    reshape(int([&
    1,2,3,4   &
    ],pInt),[maxNnodeAtIP(5),nIP(5)])

 integer(pInt), dimension(maxNnodeAtIP(6),nIP(6)), parameter, private :: NnodeAtIP6 = &
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4   &
    ],pInt),[maxNnodeAtIP(6),nIP(6)])

 integer(pInt), dimension(maxNnodeAtIP(7),nIP(7)), parameter, private :: NnodeAtIP7 = &
    reshape(int([&
    1,  &
    2,  &
    3,  &
    4,  &
    5,  &
    6   &
    ],pInt),[maxNnodeAtIP(7),nIP(7)])

 integer(pInt), dimension(maxNnodeAtIP(8),nIP(8)), parameter, private :: NnodeAtIP8 = &
    reshape(int([&
    1,2,3,4,5,6,7,8   &
    ],pInt),[maxNnodeAtIP(8),nIP(8)])

 integer(pInt), dimension(maxNnodeAtIP(9),nIP(9)), parameter, private :: NnodeAtIP9 = &
    reshape(int([&
    1,  &
    2,  &
    4,  &
    3,  &
    5,  &
    6,  &
    8,  &
    7   &
    ],pInt),[maxNnodeAtIP(9),nIP(9)])

 integer(pInt), dimension(maxNnodeAtIP(10),nIP(10)), parameter, private :: NnodeAtIP10 = &
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
    ],pInt),[maxNnodeAtIP(10),nIP(10)])

 ! *** FE_ipNeighbor ***
 ! is a list of the neighborhood of each IP.
 ! It is sorted in (local) +x,-x, +y,-y, +z,-z direction.
 ! Positive integers denote an intra-FE IP identifier.
 ! Negative integers denote the interface behind which the neighboring (extra-FE) IP will be located.


 integer(pInt), dimension(nIPneighbor(cellType(1)),nIP(1)), parameter, private :: IPneighbor1 = &
    reshape(int([&
    -2,-3,-1   &
    ],pInt),[nIPneighbor(cellType(1)),nIP(1)])

 integer(pInt), dimension(nIPneighbor(cellType(2)),nIP(2)), parameter, private :: IPneighbor2 = &
     reshape(int([&
      2,-3, 3,-1,  &
     -2, 1, 3,-1,  &
      2,-3,-2, 1   &
    ],pInt),[nIPneighbor(cellType(2)),nIP(2)])
 
 integer(pInt), dimension(nIPneighbor(cellType(3)),nIP(3)), parameter, private :: IPneighbor3 = &
     reshape(int([&
      2,-4, 3,-1,  &
     -2, 1, 4,-1,  &
      4,-4,-3, 1,  &
     -2, 3,-3, 2   &
    ],pInt),[nIPneighbor(cellType(3)),nIP(3)])
 
 integer(pInt), dimension(nIPneighbor(cellType(4)),nIP(4)), parameter, private :: IPneighbor4 = &
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
    ],pInt),[nIPneighbor(cellType(4)),nIP(4)])
 
 integer(pInt), dimension(nIPneighbor(cellType(5)),nIP(5)), parameter, private :: IPneighbor5 = &
     reshape(int([&
     -1,-2,-3,-4   &
    ],pInt),[nIPneighbor(cellType(5)),nIP(5)])
 
 integer(pInt), dimension(nIPneighbor(cellType(6)),nIP(6)), parameter, private :: IPneighbor6 = &
     reshape(int([&
      2,-4, 3,-2, 4,-1,  &
     -2, 1, 3,-2, 4,-1,  &
      2,-4,-3, 1, 4,-1,  &
      2,-4, 3,-2,-3, 1   &
    ],pInt),[nIPneighbor(cellType(6)),nIP(6)])
 
 integer(pInt), dimension(nIPneighbor(cellType(7)),nIP(7)), parameter, private :: IPneighbor7 = &
     reshape(int([&
      2,-4, 3,-2, 4,-1,  &
     -3, 1, 3,-2, 5,-1,  &
      2,-4,-3, 1, 6,-1,  &
      5,-4, 6,-2,-5, 1,  &
     -3, 4, 6,-2,-5, 2,  &
      5,-4,-3, 4,-5, 3   &
    ],pInt),[nIPneighbor(cellType(7)),nIP(7)])
 
 integer(pInt), dimension(nIPneighbor(cellType(8)),nIP(8)), parameter, private :: IPneighbor8 = &
     reshape(int([&
     -3,-5,-4,-2,-6,-1   &
    ],pInt),[nIPneighbor(cellType(8)),nIP(8)])
 
 integer(pInt), dimension(nIPneighbor(cellType(9)),nIP(9)), parameter, private :: IPneighbor9 = &
     reshape(int([&
      2,-5, 3,-2, 5,-1,  &
     -3, 1, 4,-2, 6,-1,  &
      4,-5,-4, 1, 7,-1,  &
     -3, 3,-4, 2, 8,-1,  &
      6,-5, 7,-2,-6, 1,  &
     -3, 5, 8,-2,-6, 2,  &
      8,-5,-4, 5,-6, 3,  &
     -3, 7,-4, 6,-6, 4   &
    ],pInt),[nIPneighbor(cellType(9)),nIP(9)])
 
 integer(pInt), dimension(nIPneighbor(cellType(10)),nIP(10)), parameter, private :: IPneighbor10 = &
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
    ],pInt),[nIPneighbor(cellType(10)),nIP(10)])


 real(pReal), dimension(nNode(1),NcellNode(geomType(1))), parameter :: cellNodeParentNodeWeights1 = &
    reshape(real([&
    1, 0, 0,  &
    0, 1, 0,  &
    0, 0, 1   &
    ],pReal),[nNode(1),NcellNode(geomType(1))]) ! 2D 3node 1ip

 real(pReal), dimension(nNode(2),NcellNode(geomType(2))), parameter :: cellNodeParentNodeWeights2 = &
     reshape(real([&
     1, 0, 0, 0, 0, 0,  &
     0, 1, 0, 0, 0, 0,  &
     0, 0, 1, 0, 0, 0,  &
     0, 0, 0, 1, 0, 0,  &
     0, 0, 0, 0, 1, 0,  &
     0, 0, 0, 0, 0, 1,  &
     1, 1, 1, 2, 2, 2   &
    ],pReal),[nNode(2),NcellNode(geomType(2))]) ! 2D 6node 3ip

 real(pReal), dimension(nNode(3),NcellNode(geomType(3))), parameter :: cellNodeParentNodeWeights3 = &
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
    ],pReal),[nNode(3),NcellNode(geomType(3))]) ! 2D 6node 3ip

 real(pReal), dimension(nNode(4),NcellNode(geomType(4))), parameter :: cellNodeParentNodeWeights4 = &
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
    ],pReal),[nNode(4),NcellNode(geomType(4))])  ! 2D 8node 9ip

 real(pReal), dimension(nNode(5),NcellNode(geomType(5))), parameter :: cellNodeParentNodeWeights5 = &
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
    ],pReal),[nNode(5),NcellNode(geomType(5))])  ! 2D 8node 4ip

 real(pReal), dimension(nNode(6),NcellNode(geomType(6))), parameter :: cellNodeParentNodeWeights6 = &
    reshape(real([&
    1, 0, 0, 0,  &
    0, 1, 0, 0,  &
    0, 0, 1, 0,  &
    0, 0, 0, 1   &
    ],pReal),[nNode(6),NcellNode(geomType(6))]) ! 3D 4node 1ip

 real(pReal), dimension(nNode(7),NcellNode(geomType(7))), parameter :: cellNodeParentNodeWeights7 = &
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
    ],pReal),[nNode(7),NcellNode(geomType(7))])  ! 3D 5node 4ip

 real(pReal), dimension(nNode(8),NcellNode(geomType(8))), parameter :: cellNodeParentNodeWeights8 = &
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
    ],pReal),[nNode(8),NcellNode(geomType(8))])  ! 3D 10node 4ip

 real(pReal), dimension(nNode(9),NcellNode(geomType(9))), parameter :: cellNodeParentNodeWeights9 = &
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
    ],pReal),[nNode(9),NcellNode(geomType(9))])  ! 3D 6node 6ip

 real(pReal), dimension(nNode(10),NcellNode(geomType(10))), parameter :: cellNodeParentNodeWeights10 = &
    reshape(real([&
    1, 0, 0, 0, 0, 0, 0, 0,  &
    0, 1, 0, 0, 0, 0, 0, 0,  &
    0, 0, 1, 0, 0, 0, 0, 0,  &
    0, 0, 0, 1, 0, 0, 0, 0,  &
    0, 0, 0, 0, 1, 0, 0, 0,  &
    0, 0, 0, 0, 0, 1, 0, 0,  &
    0, 0, 0, 0, 0, 0, 1, 0,  &
    0, 0, 0, 0, 0, 0, 0, 1   &
    ],pReal),[nNode(10),NcellNode(geomType(10))])  ! 3D 8node 1ip
 
 real(pReal), dimension(nNode(11),NcellNode(geomType(11))), parameter :: cellNodeParentNodeWeights11 = &
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
    ],pReal),[nNode(11),NcellNode(geomType(11))])  ! 3D 8node 8ip
 
 real(pReal), dimension(nNode(12),NcellNode(geomType(12))), parameter :: cellNodeParentNodeWeights12 = &
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
    ],pReal),[nNode(12),NcellNode(geomType(12))])  ! 3D 20node 8ip
 
 real(pReal), dimension(nNode(13),NcellNode(geomType(13))), parameter :: cellNodeParentNodeWeights13 = &
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
    ],pReal),[nNode(13),NcellNode(geomType(13))])  ! 3D 20node 27ip


 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(1)),NIP(1)), parameter :: CELL1 = &
    reshape(int([&
    1,2,3   &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(1)),NIP(1)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(2)),NIP(2)), parameter :: CELL2 = &
    reshape(int([&
    1, 4, 7, 6,   &
    2, 5, 7, 4,   &
    3, 6, 7, 5    &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(2)),NIP(2)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(3)),NIP(3)), parameter :: CELL3 = &
    reshape(int([&
    1, 5, 9, 8,   &
    5, 2, 6, 9,   &
    8, 9, 7, 4,   &
    9, 6, 3, 7    &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(3)),NIP(3)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(4)),NIP(4)), parameter :: CELL4 = &
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
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(4)),NIP(4)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(5)),NIP(5)), parameter :: CELL5 = &
    reshape(int([&
    1, 2, 3, 4   &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(5)),NIP(5)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(6)),NIP(6)), parameter :: CELL6 = &
    reshape(int([&
    1, 5,11, 7, 8,12,15,14,  &
    5, 2, 6,11,12, 9,13,15,  &
    7,11, 6, 3,14,15,13,10,  &
    8,12,15, 4, 4, 9,13,10   &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(6)),NIP(6)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(7)),NIP(7)), parameter :: CELL7 = &
    reshape(int([&
    1, 7,16, 9,10,17,21,19,  &
    7, 2, 8,16,17,11,18,21,  &
    9,16, 8, 3,19,21,18,12,  &
   10,17,21,19, 4,13,20,15,  &
   17,11,18,21,13, 5,14,20,  &
   19,21,18,12,15,20,14, 6   &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(7)),NIP(7)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(8)),NIP(8)), parameter :: CELL8 = &
    reshape(int([&
    1, 2, 3, 4, 5, 6, 7, 8   &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(8)),NIP(8)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(9)),NIP(9)), parameter :: CELL9 = &
    reshape(int([&
    1, 9,21,12,13,22,27,25,  &
    9, 2,10,21,22,14,23,27,  &
   12,21,11, 4,25,27,24,16,  &
   21,10, 3,11,27,23,15,24,  &
   13,22,27,25, 5,17,26,20,  &
   22,14,23,27,17, 6,18,26,  &
   25,27,24,16,20,26,19, 8,  &
   27,23,15,24,26,18, 7,19   &
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(9)),NIP(9)])

 integer(pInt), dimension(NCELLNODEPERCELL(CELLTYPE(10)),NIP(10)), parameter :: CELL10 = &
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
    ],pInt),[NCELLNODEPERCELL(CELLTYPE(10)),NIP(10)])


  integer(pInt), dimension(NCELLNODEPERCELLFACE(1),NIPNEIGHBOR(1)), parameter :: CELLFACE1 = &
    reshape(int([&
    2,3,  &
    3,1,  &
    1,2   &
    ],pInt),[NCELLNODEPERCELLFACE(1),NIPNEIGHBOR(1)])                   ! 2D 3node, VTK_TRIANGLE (5)

 integer(pInt), dimension(NCELLNODEPERCELLFACE(2),NIPNEIGHBOR(2)), parameter :: CELLFACE2 = &
    reshape(int([&
    2,3,  &
    4,1,  &
    3,4,  &
    1,2   &
    ],pInt),[NCELLNODEPERCELLFACE(2),NIPNEIGHBOR(2)])              ! 2D 4node, VTK_QUAD (9)

 integer(pInt), dimension(NCELLNODEPERCELLFACE(3),NIPNEIGHBOR(3)), parameter :: CELLFACE3 = &
    reshape(int([&
    1,3,2,  &
    1,2,4,  &
    2,3,4,  &
    1,4,3   &
    ],pInt),[NCELLNODEPERCELLFACE(3),NIPNEIGHBOR(3)])                     ! 3D 4node, VTK_TETRA (10)

 integer(pInt), dimension(NCELLNODEPERCELLFACE(4),NIPNEIGHBOR(4)), parameter :: CELLFACE4 = &
    reshape(int([&
    2,3,7,6,  &
    4,1,5,8,  &
    3,4,8,7,  &
    1,2,6,5,  &
    5,6,7,8,  &
    1,4,3,2   &
    ],pInt),[NCELLNODEPERCELLFACE(4),NIPNEIGHBOR(4)])              ! 3D 8node, VTK_HEXAHEDRON (12)


contains

 subroutine tElement_init(self,elemType)
   implicit none
   class(tElement) :: self
   integer(pInt), intent(in) :: elemType
   self%elemType = elemType

   self%Nnodes     = Nnode    (self%elemType)
   self%geomType   = geomType (self%elemType)
   select case (self%elemType)
     case(1_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights1
     case(2_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights2
     case(3_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights3
     case(4_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights4
     case(5_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights5
     case(6_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights6
     case(7_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights7
     case(8_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights8
     case(9_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights9
     case(10_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights10
     case(11_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights11
     case(12_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights12
     case(13_pInt)
       self%cellNodeParentNodeWeights  = cellNodeParentNodeWeights13
     case default
       print*, 'Mist'
  end select
  

   self%NcellNodes   = NcellNode    (self%geomType)
   self%maxNnodeAtIP = maxNnodeAtIP (self%geomType) 
   self%nIPs         = nIP          (self%geomType)
   self%cellType     = cellType     (self%geomType)


   select case (self%geomType)
     case(1_pInt)
       self%NnodeAtIP   = NnodeAtIP1
       self%IPneighbor  = IPneighbor1
       self%cell        = CELL1
     case(2_pInt)
       self%NnodeAtIP   = NnodeAtIP2
       self%IPneighbor  = IPneighbor2
       self%cell        = CELL2
     case(3_pInt)
       self%NnodeAtIP   = NnodeAtIP3
       self%IPneighbor  = IPneighbor3
       self%cell        = CELL3
     case(4_pInt)
       self%NnodeAtIP   = NnodeAtIP4
       self%IPneighbor  = IPneighbor4
       self%cell        = CELL4
     case(5_pInt)
       self%NnodeAtIP   = NnodeAtIP5
       self%IPneighbor  = IPneighbor5
       self%cell        = CELL5
     case(6_pInt)
       self%NnodeAtIP   = NnodeAtIP6
       self%IPneighbor  = IPneighbor6
       self%cell        = CELL6
     case(7_pInt)
       self%NnodeAtIP   = NnodeAtIP7
       self%IPneighbor  = IPneighbor7
       self%cell        = CELL7
     case(8_pInt)
       self%NnodeAtIP   = NnodeAtIP8
       self%IPneighbor  = IPneighbor8
       self%cell        = CELL8
     case(9_pInt)
       self%NnodeAtIP   = NnodeAtIP9
       self%IPneighbor  = IPneighbor9
       self%cell        = CELL9
     case(10_pInt)
       self%NnodeAtIP   = NnodeAtIP10
       self%IPneighbor  = IPneighbor10
       self%cell        = CELL10
  end select
   self%NcellNodesPerCell = NCELLNODEPERCELL(self%cellType)
 
  select case(self%cellType)
    case(1_pInt)
      self%cellFace = CELLFACE1
    case(2_pInt)
      self%cellFace = CELLFACE2
    case(3_pInt)
      self%cellFace = CELLFACE3
    case(4_pInt)
      self%cellFace = CELLFACE4
   end select 
 end subroutine tElement_init



end module element
