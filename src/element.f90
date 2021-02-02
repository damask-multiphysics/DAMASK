!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Koords, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module element
  use IO

  implicit none
  private

!---------------------------------------------------------------------------------------------------
!> Properties of a single element
!---------------------------------------------------------------------------------------------------
  type, public :: tElement
    integer :: &
      elemType, &
      geomType, &                                                                                   !< geometry type (same for same dimension and same number of integration points)
      cellType, &
      Nnodes, &
      Ncellnodes, &
      NcellnodesPerCell, &
      nIPs, &
      nIPneighbors
    character(len=:),        allocatable :: &
      vtkType
    integer, dimension(:,:), allocatable :: &
      Cell, &                                                                                       !< intra-element (cell) nodes that constitute a cell
      IPneighbor, &
      cellFace
    integer, dimension(:,:), allocatable :: &
      ! center of gravity of the weighted nodes gives the position of the cell node.
      ! example: face-centered cell node with face nodes 1,2,5,6 to be used in,
      !          e.g., an 8 node element, would be encoded: 1, 1, 0, 0, 1, 1, 0, 0
      cellNodeParentNodeWeights
    contains
      procedure :: init => tElement_init
  end type tElement


  integer, parameter :: &
    NELEMTYPE = 13

  integer, dimension(NELEMTYPE), parameter :: NNODE = &
    [ &
        3, & ! 2D, 1 IP
        6, & ! 2D, 3 IP
        4, & ! 2D, 4 IP
        8, & ! 2D, 9 IP
        8, & ! 2D, 4 IP
        !----------------------
        4, & ! 3D, 1 IP
        5, & ! 3D, 4 IP
       10, & ! 3D, 4 IP
        6, & ! 3D, 6 IP
        8, & ! 3D, 1 IP
        8, & ! 3D, 8 IP
       20, & ! 3D, 8 IP
       20  & ! 3D, 27 IP
    ]                                                                                               !< number of nodes that constitute a specific type of element

  integer, dimension(NELEMTYPE), parameter :: GEOMTYPE = &
    [ &
        1, & ! 1 triangle
        2, & ! 3 quadrilaterals
        3, & ! 4 quadrilaterals
        4, & ! 9 quadrilaterals
        3, & ! 4 quadrilaterals
        !----------------------
        5, & ! 1 tetrahedron
        6, & ! 4 hexahedrons
        6, & ! 4 hexahedrons
        7, & ! 6 hexahedrons
        8, & ! 1 hexahedron
        9, & ! 8 hexahedrons
        9, & ! 8 hexahedrons
       10  & ! 27 hexahedrons
    ]                                                                                               !< geometry type (same number of cell nodes and IPs)

  integer, dimension(maxval(GEOMTYPE)), parameter :: NCELLNODE = &
    [ &
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
    ]                                                                                               !< number of cell nodes

  integer, dimension(maxval(GEOMTYPE)), parameter :: NIP = &
    [ &
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
    ]                                                                                               !< number of IPs

  integer, dimension(maxval(GEOMTYPE)), parameter :: CELLTYPE = &
    [ &
       1, & ! 2D, 3 node (Triangle)
       2, & ! 2D, 4 node (Quadrilateral)
       2, & !   - " -
       2, & !   - " -
       3, & ! 3D, 4 node (Tetrahedron)
       4, & ! 3D, 4 node (Hexahedron)
       4, & !   - " -
       4, & !   - " -
       4, & !   - " -
       4  & !   - " -
    ]                                                                                               !< cell type

  integer, dimension(maxval(CELLTYPE)), parameter :: NIPNEIGHBOR = &
    [ &
       3, &
       4, &
       4, &
       6  &
    ]                                                                                               !< number of ip neighbors / cell faces

  integer, dimension(maxval(CELLTYPE)), parameter :: NCELLNODEPERCELLFACE = &
    [ &
       2, &
       2, &
       3, &
       4  &
    ]                                                                                               !< number of cell nodes per face

  integer, dimension(maxval(CELLTYPE)), parameter  :: NCELLNODEPERCELL = &
    [ &
       3, &
       4, &
       4, &
       8  &
    ]                                                                                               !< number of total cell nodes

  ! *** IPneighbor ***
  ! list of the neighborhood of each IP.
  ! It is sorted in (local) +x,-x, +y,-y, +z,-z direction.
  ! Positive integers denote an intra-element IP identifier.
  ! Negative integers denote the interface behind which the neighboring (extra-element) IP will be located.

  integer, dimension(NIPNEIGHBOR(CELLTYPE(1)),NIP(1)),   parameter :: IPNEIGHBOR1 = &
    reshape([&
      -2,-3,-1  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR1))
#else
    ],[NIPNEIGHBOR(CELLTYPE(1)),NIP(1)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(2)),NIP(2)),   parameter :: IPNEIGHBOR2 = &
    reshape([&
       2,-3, 3,-1, &
      -2, 1, 3,-1, &
       2,-3,-2, 1  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR2))
#else
    ],[NIPNEIGHBOR(CELLTYPE(2)),NIP(2)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(3)),NIP(3)),   parameter :: IPNEIGHBOR3 = &
    reshape([&
       2,-4, 3,-1, &
      -2, 1, 4,-1, &
       4,-4,-3, 1, &
      -2, 3,-3, 2  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR3))
#else
    ],[NIPNEIGHBOR(CELLTYPE(3)),NIP(3)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(4)),NIP(4)),   parameter :: IPNEIGHBOR4 = &
    reshape([&
       2,-4, 4,-1, &
       3, 1, 5,-1, &
      -2, 2, 6,-1, &
       5,-4, 7, 1, &
       6, 4, 8, 2, &
      -2, 5, 9, 3, &
       8,-4,-3, 4, &
       9, 7,-3, 5, &
      -2, 8,-3, 6  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR4))
#else
    ],[NIPNEIGHBOR(CELLTYPE(4)),NIP(4)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(5)),NIP(5)),   parameter :: IPNEIGHBOR5 = &
    reshape([&
      -1,-2,-3,-4 &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR5))
#else
    ],[NIPNEIGHBOR(CELLTYPE(5)),NIP(5)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(6)),NIP(6)),   parameter :: IPNEIGHBOR6 = &
    reshape([&
       2,-4, 3,-2, 4,-1, &
      -2, 1, 3,-2, 4,-1, &
       2,-4,-3, 1, 4,-1, &
       2,-4, 3,-2,-3, 1  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR6))
#else
    ],[NIPNEIGHBOR(CELLTYPE(6)),NIP(6)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(7)),NIP(7)),   parameter :: IPNEIGHBOR7 = &
    reshape([&
       2,-4, 3,-2, 4,-1, &
      -3, 1, 3,-2, 5,-1, &
       2,-4,-3, 1, 6,-1, &
       5,-4, 6,-2,-5, 1, &
      -3, 4, 6,-2,-5, 2, &
       5,-4,-3, 4,-5, 3  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR7))
#else
    ],[NIPNEIGHBOR(CELLTYPE(7)),NIP(7)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(8)),NIP(8)),   parameter :: IPNEIGHBOR8 = &
    reshape([&
      -3,-5,-4,-2,-6,-1  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR8))
#else
    ],[NIPNEIGHBOR(CELLTYPE(8)),NIP(8)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(9)),NIP(9)),   parameter :: IPNEIGHBOR9 = &
    reshape([&
       2,-5, 3,-2, 5,-1, &
      -3, 1, 4,-2, 6,-1, &
       4,-5,-4, 1, 7,-1, &
      -3, 3,-4, 2, 8,-1, &
       6,-5, 7,-2,-6, 1, &
      -3, 5, 8,-2,-6, 2, &
       8,-5,-4, 5,-6, 3, &
      -3, 7,-4, 6,-6, 4  &
#if !defined(__GFORTRAN__)
    ],shape(IPNEIGHBOR9))
#else
    ],[NIPNEIGHBOR(CELLTYPE(9)),NIP(9)])
#endif

  integer, dimension(NIPNEIGHBOR(CELLTYPE(10)),NIP(10)), parameter :: IPNEIGHBOR10 = &
    reshape([&
       2,-5, 4,-2,10,-1, &
       3, 1, 5,-2,11,-1, &
      -3, 2, 6,-2,12,-1, &
       5,-5, 7, 1,13,-1, &
       6, 4, 8, 2,14,-1, &
      -3, 5, 9, 3,15,-1, &
       8,-5,-4, 4,16,-1, &
       9, 7,-4, 5,17,-1, &
      -3, 8,-4, 6,18,-1, &
      11,-5,13,-2,19, 1, &
      12,10,14,-2,20, 2, &
      -3,11,15,-2,21, 3, &
      14,-5,16,10,22, 4, &
      15,13,17,11,23, 5, &
      -3,14,18,12,24, 6, &
      17,-5,-4,13,25, 7, &
      18,16,-4,14,26, 8, &
      -3,17,-4,15,27, 9, &
      20,-5,22,-2,-6,10, &
      21,19,23,-2,-6,11, &
      -3,20,24,-2,-6,12, &
      23,-5,25,19,-6,13, &
      24,22,26,20,-6,14, &
      -3,23,27,21,-6,15, &
      26,-5,-4,22,-6,16, &
      27,25,-4,23,-6,17, &
      -3,26,-4,24,-6,18  &
#if !defined(__GFORTRAN__)
     ],shape(IPNEIGHBOR10))
#else
    ],[NIPNEIGHBOR(CELLTYPE(10)),NIP(10)])
#endif


  integer, dimension(NNODE(1),NCELLNODE(GEOMTYPE(1))),   parameter :: CELLNODEPARENTNODEWEIGHTS1 = &
    reshape([&
      1, 0, 0, &
      0, 1, 0, &
      0, 0, 1  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS1))                                                            !< 2D 3node 1ip
#else
    ],[NNODE(1),NCELLNODE(GEOMTYPE(1))])
#endif

  integer, dimension(NNODE(2),NCELLNODE(GEOMTYPE(2))),   parameter :: CELLNODEPARENTNODEWEIGHTS2 = &
    reshape([&
      1, 0, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, 0, &
      0, 0, 1, 0, 0, 0, &
      0, 0, 0, 1, 0, 0, &
      0, 0, 0, 0, 1, 0, &
      0, 0, 0, 0, 0, 1, &
      1, 1, 1, 2, 2, 2  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS2))                                                            !< 2D 6node 3ip
#else
    ],[NNODE(2),NCELLNODE(GEOMTYPE(2))])
#endif

  integer, dimension(NNODE(3),NCELLNODE(GEOMTYPE(3))),   parameter :: CELLNODEPARENTNODEWEIGHTS3 = &
    reshape([&
      1, 0, 0, 0, &
      0, 1, 0, 0, &
      0, 0, 1, 0, &
      0, 0, 0, 1, &
      1, 1, 0, 0, &
      0, 1, 1, 0, &
      0, 0, 1, 1, &
      1, 0, 0, 1, &
      1, 1, 1, 1  &
#if !defined(__GFORTRAN__)
     ],shape(CELLNODEPARENTNODEWEIGHTS3))                                                            !< 2D 6node 3ip
#else
    ],[NNODE(3),NCELLNODE(GEOMTYPE(3))])
#endif

  integer, dimension(NNODE(4),NCELLNODE(GEOMTYPE(4))),   parameter :: CELLNODEPARENTNODEWEIGHTS4 = &
    reshape([&
      1, 0, 0, 0, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, 0, 0, 0, &
      0, 0, 1, 0, 0, 0, 0, 0, &
      0, 0, 0, 1, 0, 0, 0, 0, &
      1, 0, 0, 0, 2, 0, 0, 0, &
      0, 1, 0, 0, 2, 0, 0, 0, &
      0, 1, 0, 0, 0, 2, 0, 0, &
      0, 0, 1, 0, 0, 2, 0, 0, &
      0, 0, 1, 0, 0, 0, 2, 0, &
      0, 0, 0, 1, 0, 0, 2, 0, &
      0, 0, 0, 1, 0, 0, 0, 2, &
      1, 0, 0, 0, 0, 0, 0, 2, &
      4, 1, 1, 1, 8, 2, 2, 8, &
      1, 4, 1, 1, 8, 8, 2, 2, &
      1, 1, 4, 1, 2, 8, 8, 2, &
      1, 1, 1, 4, 2, 2, 8, 8  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS4))                                                            !< 2D 8node 9ip
#else
    ],[NNODE(4),NCELLNODE(GEOMTYPE(4))])
#endif

  integer, dimension(NNODE(5),NCELLNODE(GEOMTYPE(5))),   parameter :: CELLNODEPARENTNODEWEIGHTS5 = &
    reshape([&
      1, 0, 0, 0, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, 0, 0, 0, &
      0, 0, 1, 0, 0, 0, 0, 0, &
      0, 0, 0, 1, 0, 0, 0, 0, &
      0, 0, 0, 0, 1, 0, 0, 0, &
      0, 0, 0, 0, 0, 1, 0, 0, &
      0, 0, 0, 0, 0, 0, 1, 0, &
      0, 0, 0, 0, 0, 0, 0, 1, &
      1, 1, 1, 1, 2, 2, 2, 2  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS5))                                                            !< 2D 8node 4ip
#else
    ],[NNODE(5),NCELLNODE(GEOMTYPE(5))])
#endif

  integer, dimension(NNODE(6),NcellNode(GEOMTYPE(6))),   parameter :: CELLNODEPARENTNODEWEIGHTS6 = &
    reshape([&
      1, 0, 0, 0, &
      0, 1, 0, 0, &
      0, 0, 1, 0, &
      0, 0, 0, 1  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS6))                                                            !< 3D 4node 1ip
#else
    ],[NNODE(6),NcellNode(GEOMTYPE(6))])
#endif

  integer, dimension(NNODE(7),NCELLNODE(GEOMTYPE(7))),   parameter :: CELLNODEPARENTNODEWEIGHTS7 = &
    reshape([&
      1, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, &
      0, 0, 1, 0, 0, &
      0, 0, 0, 1, 0, &
      1, 1, 0, 0, 0, &
      0, 1, 1, 0, 0, &
      1, 0, 1, 0, 0, &
      1, 0, 0, 1, 0, &
      0, 1, 0, 1, 0, &
      0, 0, 1, 1, 0, &
      1, 1, 1, 0, 0, &
      1, 1, 0, 1, 0, &
      0, 1, 1, 1, 0, &
      1, 0, 1, 1, 0, &
      0, 0, 0, 0, 1  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS7))                                                            !< 3D 5node 4ip
#else
    ],[NNODE(7),NCELLNODE(GEOMTYPE(7))])
#endif

  integer, dimension(NNODE(8),NCELLNODE(GEOMTYPE(8))),   parameter :: CELLNODEPARENTNODEWEIGHTS8 = &
    reshape([&
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, &
      0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
      1, 1, 1, 0, 2, 2, 2, 0, 0, 0, &
      1, 1, 0, 1, 2, 0, 0, 2, 2, 0, &
      0, 1, 1, 1, 0, 2, 0, 0, 2, 2, &
      1, 0, 1, 1, 0, 0, 2, 2, 0, 2, &
      3, 3, 3, 3, 4, 4, 4, 4, 4, 4  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS8))                                                            !< 3D 10node 4ip
#else
    ],[NNODE(8),NCELLNODE(GEOMTYPE(8))])
#endif

  integer, dimension(NNODE(9),NCELLNODE(GEOMTYPE(9))),   parameter :: CELLNODEPARENTNODEWEIGHTS9 = &
    reshape([&
      1, 0, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, 0, &
      0, 0, 1, 0, 0, 0, &
      0, 0, 0, 1, 0, 0, &
      0, 0, 0, 0, 1, 0, &
      0, 0, 0, 0, 0, 1, &
      1, 1, 0, 0, 0, 0, &
      0, 1, 1, 0, 0, 0, &
      1, 0, 1, 0, 0, 0, &
      1, 0, 0, 1, 0, 0, &
      0, 1, 0, 0, 1, 0, &
      0, 0, 1, 0, 0, 1, &
      0, 0, 0, 1, 1, 0, &
      0, 0, 0, 0, 1, 1, &
      0, 0, 0, 1, 0, 1, &
      1, 1, 1, 0, 0, 0, &
      1, 1, 0, 1, 1, 0, &
      0, 1, 1, 0, 1, 1, &
      1, 0, 1, 1, 0, 1, &
      0, 0, 0, 1, 1, 1, &
      1, 1, 1, 1, 1, 1  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS9))                                                            !< 3D 6node 6ip
#else
    ],[NNODE(9),NCELLNODE(GEOMTYPE(9))])
#endif

  integer, dimension(NNODE(10),NCELLNODE(GEOMTYPE(10))), parameter :: CELLNODEPARENTNODEWEIGHTS10 = &
    reshape([&
      1, 0, 0, 0, 0, 0, 0, 0, &
      0, 1, 0, 0, 0, 0, 0, 0, &
      0, 0, 1, 0, 0, 0, 0, 0, &
      0, 0, 0, 1, 0, 0, 0, 0, &
      0, 0, 0, 0, 1, 0, 0, 0, &
      0, 0, 0, 0, 0, 1, 0, 0, &
      0, 0, 0, 0, 0, 0, 1, 0, &
      0, 0, 0, 0, 0, 0, 0, 1  &
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS10))                                                           !< 3D 8node 1ip
#else
    ],[NNODE(10),NCELLNODE(GEOMTYPE(10))])
#endif

  integer, dimension(NNODE(11),NCELLNODE(GEOMTYPE(11))), parameter :: CELLNODEPARENTNODEWEIGHTS11 = &
    reshape([&
      1, 0, 0, 0,  0, 0, 0, 0, &   !
      0, 1, 0, 0,  0, 0, 0, 0, &   !
      0, 0, 1, 0,  0, 0, 0, 0, &   !
      0, 0, 0, 1,  0, 0, 0, 0, &   !
      0, 0, 0, 0,  1, 0, 0, 0, &   !  5
      0, 0, 0, 0,  0, 1, 0, 0, &   !
      0, 0, 0, 0,  0, 0, 1, 0, &   !
      0, 0, 0, 0,  0, 0, 0, 1, &   !
      1, 1, 0, 0,  0, 0, 0, 0, &   !
      0, 1, 1, 0,  0, 0, 0, 0, &   ! 10
      0, 0, 1, 1,  0, 0, 0, 0, &   !
      1, 0, 0, 1,  0, 0, 0, 0, &   !
      0, 0, 0, 0,  1, 1, 0, 0, &   !
      0, 0, 0, 0,  0, 1, 1, 0, &   !
      0, 0, 0, 0,  0, 0, 1, 1, &   ! 15
      0, 0, 0, 0,  1, 0, 0, 1, &   !
      1, 0, 0, 0,  1, 0, 0, 0, &   !
      0, 1, 0, 0,  0, 1, 0, 0, &   !
      0, 0, 1, 0,  0, 0, 1, 0, &   !
      0, 0, 0, 1,  0, 0, 0, 1, &   ! 20
      1, 1, 1, 1,  0, 0, 0, 0, &   !
      1, 1, 0, 0,  1, 1, 0, 0, &   !
      0, 1, 1, 0,  0, 1, 1, 0, &   !
      0, 0, 1, 1,  0, 0, 1, 1, &   !
      1, 0, 0, 1,  1, 0, 0, 1, &   ! 25
      0, 0, 0, 0,  1, 1, 1, 1, &   !
      1, 1, 1, 1,  1, 1, 1, 1  &   !
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS11))                                                           !< 3D 8node 8ip
#else
    ],[NNODE(11),NCELLNODE(GEOMTYPE(11))])
#endif

  integer, dimension(NNODE(12),NCELLNODE(GEOMTYPE(12))), parameter :: CELLNODEPARENTNODEWEIGHTS12 = &
    reshape([&
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
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS12))                                                           !< 3D 20node 8ip
#else
    ],[NNODE(12),NCELLNODE(GEOMTYPE(12))])
#endif

  integer, dimension(NNODE(13),NCELLNODE(GEOMTYPE(13))), parameter :: CELLNODEPARENTNODEWEIGHTS13 = &
    reshape([&
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
#if !defined(__GFORTRAN__)
    ],shape(CELLNODEPARENTNODEWEIGHTS13))                                                           !< 3D 20node 27ip
#else
    ],[NNODE(13),NCELLNODE(GEOMTYPE(13))])
#endif


  integer, dimension(NCELLNODEPERCELL(CELLTYPE(1)),NIP(1)),   parameter :: CELL1 = &
    reshape([&
      1,2,3 &
#if !defined(__GFORTRAN__)
    ],shape(CELL1))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(1)),NIP(1)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(2)),NIP(2)),   parameter :: CELL2 = &
    reshape([&
      1, 4, 7, 6, &
      2, 5, 7, 4, &
      3, 6, 7, 5  &
#if !defined(__GFORTRAN__)
    ],shape(CELL2))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(2)),NIP(2)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(3)),NIP(3)),   parameter :: CELL3 = &
    reshape([&
      1, 5, 9, 8, &
      5, 2, 6, 9, &
      8, 9, 7, 4, &
      9, 6, 3, 7  &
#if !defined(__GFORTRAN__)
    ],shape(CELL3))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(3)),NIP(3)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(4)),NIP(4)),   parameter :: CELL4 = &
    reshape([&
       1, 5,13,12, &
       5, 6,14,13, &
       6, 2, 7,14, &
      12,13,16,11, &
      13,14,15,16, &
      14, 7, 8,15, &
      11,16,10, 4, &
      16,15, 9,10, &
      15, 8, 3, 9  &
#if !defined(__GFORTRAN__)
    ],shape(CELL4))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(4)),NIP(4)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(5)),NIP(5)),   parameter :: CELL5 = &
    reshape([&
      1, 2, 3, 4 &
#if !defined(__GFORTRAN__)
    ],shape(CELL5))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(5)),NIP(5)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(6)),NIP(6)),   parameter :: CELL6 = &
    reshape([&
      1, 5,11, 7, 8,12,15,14, &
      5, 2, 6,11,12, 9,13,15, &
      7,11, 6, 3,14,15,13,10, &
      8,12,15,14, 4, 9,13,10  &
#if !defined(__GFORTRAN__)
    ],shape(CELL6))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(6)),NIP(6)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(7)),NIP(7)),   parameter :: CELL7 = &
    reshape([&
      1, 7,16, 9,10,17,21,19, &
      7, 2, 8,16,17,11,18,21, &
      9,16, 8, 3,19,21,18,12, &
     10,17,21,19, 4,13,20,15, &
     17,11,18,21,13, 5,14,20, &
     19,21,18,12,15,20,14, 6  &
#if !defined(__GFORTRAN__)
    ],shape(CELL7))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(7)),NIP(7)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(8)),NIP(8)),   parameter :: CELL8 = &
    reshape([&
      1, 2, 3, 4, 5, 6, 7, 8 &
#if !defined(__GFORTRAN__)
    ],shape(CELL8))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(8)),NIP(8)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(9)),NIP(9)),   parameter :: CELL9 = &
    reshape([&
       1, 9,21,12,17,22,27,25, &
       9, 2,10,21,22,18,23,27, &
      12,21,11, 4,25,27,24,20, &
      21,10, 3,11,27,23,19,24, &
      17,22,27,25, 5,13,26,16, &
      22,18,23,27,13, 6,14,26, &
      25,27,24,20,16,26,15, 8, &
      27,23,19,24,26,14, 7,15  &
#if !defined(__GFORTRAN__)
    ],shape(CELL9))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(9)),NIP(9)])
#endif

  integer, dimension(NCELLNODEPERCELL(CELLTYPE(10)),NIP(10)), parameter :: CELL10 = &
    reshape([&
       1, 9,33,16,17,37,57,44, &
       9,10,34,33,37,38,58,57, &
      10, 2,11,34,38,18,39,58, &
      16,33,36,15,44,57,60,43, &
      33,34,35,36,57,58,59,60, &
      34,11,12,35,58,39,40,59, &
      15,36,14, 4,43,60,42,20, &
      36,35,13,14,60,59,41,42, &
      35,12, 3,13,59,40,19,41, &
      17,37,57,44,21,45,61,52, &
      37,38,58,57,45,46,62,61, &
      38,18,39,58,46,22,47,62, &
      44,57,60,43,52,61,64,51, &
      57,58,59,60,61,62,63,64, &
      58,39,40,59,62,47,48,63, &
      43,60,42,20,51,64,50,24, &
      60,59,41,42,64,63,49,50, &
      59,40,19,41,63,48,23,49, &
      21,45,61,52, 5,25,53,32, &
      45,46,62,61,25,26,54,53, &
      46,22,47,62,26, 6,27,54, &
      52,61,64,51,32,53,56,31, &
      61,62,63,64,53,54,55,56, &
      62,47,48,63,54,27,28,55, &
      51,64,50,24,31,56,30, 8, &
      64,63,49,50,56,55,29,30, &
      63,48,23,49,55,28, 7,29  &
#if !defined(__GFORTRAN__)
    ],shape(CELL10))
#else
    ],[NCELLNODEPERCELL(CELLTYPE(10)),NIP(10)])
#endif


  integer, dimension(NCELLNODEPERCELLFACE(1),NIPNEIGHBOR(1)), parameter :: CELLFACE1 = &
    reshape([&
      2,3, &
      3,1, &
      1,2  &
#if !defined(__GFORTRAN__)
    ],shape(CELLFACE1))                                                                             !< 2D 3node, VTK_TRIANGLE (5)
#else
    ],[NCELLNODEPERCELLFACE(1),NIPNEIGHBOR(1)])
#endif

  integer, dimension(NCELLNODEPERCELLFACE(2),NIPNEIGHBOR(2)), parameter :: CELLFACE2 = &
    reshape([&
      2,3, &
      4,1, &
      3,4, &
      1,2  &
#if !defined(__GFORTRAN__)
    ],shape(CELLFACE2))                                                                             !< 2D 4node, VTK_QUAD (9)
#else
    ],[NCELLNODEPERCELLFACE(2),NIPNEIGHBOR(2)])
#endif

  integer, dimension(NCELLNODEPERCELLFACE(3),NIPNEIGHBOR(3)), parameter :: CELLFACE3 = &
    reshape([&
      1,3,2, &
      1,2,4, &
      2,3,4, &
      1,4,3  &
#if !defined(__GFORTRAN__)
    ],shape(CELLFACE3))                                                                             !< 3D 4node, VTK_TETRA (10)
#else
    ],[NCELLNODEPERCELLFACE(3),NIPNEIGHBOR(3)])
#endif

  integer, dimension(NCELLNODEPERCELLFACE(4),NIPNEIGHBOR(4)), parameter :: CELLFACE4 = &
    reshape([&
      2,3,7,6, &
      4,1,5,8, &
      3,4,8,7, &
      1,2,6,5, &
      5,6,7,8, &
      1,4,3,2  &
#if !defined(__GFORTRAN__)
    ],shape(CELLFACE4))                                                                             !< 3D 8node, VTK_HEXAHEDRON (12)
#else
    ],[NCELLNODEPERCELLFACE(4),NIPNEIGHBOR(4)])
#endif


contains


!---------------------------------------------------------------------------------------------------
!> define properties of an element
!---------------------------------------------------------------------------------------------------
subroutine tElement_init(self,elemType)

  class(tElement)     :: self
  integer, intent(in) :: elemType

  self%elemType = elemType

  self%Nnodes     = NNODE   (self%elemType)
  self%geomType   = GEOMTYPE(self%elemType)

  select case (self%elemType)
    case(1)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS1
    case(2)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS2
    case(3)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS3
    case(4)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS4
    case(5)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS5
    case(6)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS6
    case(7)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS7
    case(8)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS8
    case(9)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS9
    case(10)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS10
    case(11)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS11
    case(12)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS12
    case(13)
      self%cellNodeParentNodeWeights = CELLNODEPARENTNODEWEIGHTS13
    case default
      call IO_error(0,ext_msg='invalid element type')
  end select


  self%NcellNodes = NCELLNODE(self%geomType)
  self%nIPs       = NIP      (self%geomType)
  self%cellType   = CELLTYPE (self%geomType)

  select case (self%geomType)
    case(1)
      self%IPneighbor = IPNEIGHBOR1
      self%cell       = CELL1
    case(2)
      self%IPneighbor = IPNEIGHBOR2
      self%cell       = CELL2
    case(3)
      self%IPneighbor = IPNEIGHBOR3
      self%cell       = CELL3
    case(4)
      self%IPneighbor = IPNEIGHBOR4
      self%cell       = CELL4
    case(5)
      self%IPneighbor = IPNEIGHBOR5
      self%cell       = CELL5
    case(6)
      self%IPneighbor = IPNEIGHBOR6
      self%cell       = CELL6
    case(7)
      self%IPneighbor = IPNEIGHBOR7
      self%cell       = CELL7
    case(8)
      self%IPneighbor = IPNEIGHBOR8
      self%cell       = CELL8
    case(9)
      self%IPneighbor = IPNEIGHBOR9
      self%cell       = CELL9
    case(10)
      self%IPneighbor = IPNEIGHBOR10
      self%cell       = CELL10
  end select

  self%NcellnodesPerCell = NCELLNODEPERCELL(self%cellType)

  select case(self%cellType)
    case(1)
      self%cellFace = CELLFACE1
      self%vtkType  = 'TRIANGLE'
    case(2)
      self%cellFace = CELLFACE2
      self%vtkType  = 'QUAD'
    case(3)
      self%cellFace = CELLFACE3
      self%vtkType  = 'TETRA'
    case(4)
      self%cellFace = CELLFACE4
      self%vtkType  = 'HEXAHEDRON'
  end select

  self%nIPneighbors = size(self%IPneighbor,1)

  print'(/,a)', ' <<<+-  element_init  -+>>>'; flush(IO_STDOUT)

  print*, 'element type:      ',self%elemType
  print*, '  geom type:       ',self%geomType
  print*, '  cell type:       ',self%cellType
  print*, '  # node:          ',self%Nnodes
  print*, '  # IP:            ',self%nIPs
  print*, '  # cellnode:      ',self%Ncellnodes
  print*, '  # cellnode/cell: ',self%NcellnodesPerCell
  print*, '  # IP neighbor:   ',self%nIPneighbors

end subroutine tElement_init

end module element
