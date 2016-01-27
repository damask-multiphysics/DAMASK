!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  defines lattice structure definitions, slip and twin system definitions, Schimd matrix
!>         calculation and non-Schmid behavior
!--------------------------------------------------------------------------------------------------
module lattice
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt), parameter, public :: &
   LATTICE_maxNslipFamily     =  13_pInt, &                                                         !< max # of slip system families over lattice structures
   LATTICE_maxNtwinFamily     =  4_pInt, &                                                          !< max # of twin system families over lattice structures
   LATTICE_maxNtransFamily    =  2_pInt, &                                                          !< max # of transformation system families over lattice structures
   LATTICE_maxNcleavageFamily =  3_pInt, &                                                          !< max # of transformation system families over lattice structures
   LATTICE_maxNslip           = 52_pInt, &                                                          !< max # of slip systems over lattice structures
   LATTICE_maxNtwin           = 24_pInt, &                                                          !< max # of twin systems over lattice structures
   LATTICE_maxNinteraction    = 182_pInt, &                                                         !< max # of interaction types (in hardening matrix part)
   LATTICE_maxNnonSchmid      = 6_pInt, &                                                           !< max # of non schmid contributions over lattice structures
   LATTICE_maxNtrans          = 12_pInt, &                                                          !< max # of transformations over lattice structures
   LATTICE_maxNcleavage       = 9_pInt                                                              !< max # of cleavage over lattice structures
 
 integer(pInt), allocatable, dimension(:,:), protected, public :: &
   lattice_NslipSystem, &                                                                           !< total # of slip systems in each family
   lattice_NtwinSystem, &                                                                           !< total # of twin systems in each family
   lattice_NtransSystem, &                                                                          !< total # of transformation systems in each family
   lattice_NcleavageSystem                                                                          !< total # of transformation systems in each family

 integer(pInt), allocatable, dimension(:,:,:), protected, public :: &
   lattice_interactionSlipSlip, &                                                                   !< Slip--slip interaction type 
   lattice_interactionSlipTwin, &                                                                   !< Slip--twin interaction type 
   lattice_interactionTwinSlip, &                                                                   !< Twin--slip interaction type 
   lattice_interactionTwinTwin, &                                                                   !< Twin--twin interaction type 
   lattice_interactionSlipTrans, &                                                                  !< Slip--trans interaction type 
   lattice_interactionTransSlip, &                                                                  !< Trans--slip interaction type 
   lattice_interactionTransTrans                                                                    !< Trans--trans interaction type

 real(pReal), allocatable, dimension(:,:,:,:,:), protected, public :: &
   lattice_Sslip, &                                                                                 !< Schmid and non-Schmid matrices
   lattice_Scleavage                                                                                !< Schmid matrices for cleavage systems
  
 real(pReal), allocatable, dimension(:,:,:,:), protected, public :: &
   lattice_Sslip_v, &                                                                               !< Mandel notation of lattice_Sslip
   lattice_Scleavage_v                                                                              !< Mandel notation of lattice_Scleavege
  
 real(pReal), allocatable, dimension(:,:,:), protected, public :: &
   lattice_sn, &                                                                                    !< normal direction of slip system
   lattice_sd, &                                                                                    !< slip direction of slip system
   lattice_st                                                                                       !< sd x sn

! rotation and Schmid matrices, normal, shear direction and d x n of twin systems
 real(pReal), allocatable, dimension(:,:,:,:), protected, public  :: &
   lattice_Stwin, &
   lattice_Qtwin

 real(pReal), allocatable, dimension(:,:,:), protected, public :: &
   lattice_Stwin_v, &
   lattice_tn, &
   lattice_td, &
   lattice_tt

 real(pReal), allocatable, dimension(:,:,:), protected, public :: &
   lattice_Strans_v, &                                                                              !< Eigendeformation tensor in vector form
   lattice_projectionTrans                                                                          !< Matrix for projection of slip to fault-band (twin) systems for strain-induced martensite nucleation

 real(pReal), allocatable, dimension(:,:,:,:), protected, public :: &
   lattice_Qtrans, &                                                                                !< Total rotation: Q = R*B
   lattice_Strans                                                                                   !< Eigendeformation tensor for phase transformation

 real(pReal), allocatable, dimension(:,:), protected, public :: &
   lattice_shearTwin, &                                                                             !< characteristic twin shear
   lattice_shearTrans                                                                               !< characteristic transformation shear

 integer(pInt), allocatable, dimension(:), protected, public :: &
   lattice_NnonSchmid                                                                               !< total # of non-Schmid contributions for each structure

!--------------------------------------------------------------------------------------------------
! fcc
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, public :: & 
   LATTICE_fcc_NslipSystem = int([12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],pInt)                     !< total # of slip systems per family for fcc
   
 integer(pInt), dimension(LATTICE_maxNtwinFamily), parameter, public :: &
   LATTICE_fcc_NtwinSystem = int([12, 0, 0, 0],pInt)                                                !< total # of twin systems per family for fcc

 integer(pInt), dimension(LATTICE_maxNtransFamily), parameter, public :: &
   LATTICE_fcc_NtransSystem = int([12, 0],pInt)                                                     !< total # of transformation systems per family for fcc
   
 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, public :: &
   LATTICE_fcc_NcleavageSystem = int([3, 4, 0],pInt)                                                !< total # of cleavage systems per family for fcc
   
 integer(pInt), parameter, private  :: &
   LATTICE_fcc_Nslip = 12_pInt, & ! sum(lattice_fcc_NslipSystem), &                                 !< total # of slip systems for fcc
   LATTICE_fcc_Ntwin = 12_pInt, & ! sum(lattice_fcc_NtwinSystem)                                    !< total # of twin systems for fcc
   LATTICE_fcc_NnonSchmid = 0_pInt, &                                                               !< total # of non-Schmid contributions for fcc
   LATTICE_fcc_Ntrans     = 12_pInt, &                                                              !< total # of transformations for fcc
   LATTICE_fcc_Ncleavage  = 7_pInt                                                                  !< total # of cleavage systems for fcc
 
 real(pReal), dimension(3+3,LATTICE_fcc_Nslip), parameter, private :: &
   LATTICE_fcc_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
      0, 1,-1,     1, 1, 1, &
     -1, 0, 1,     1, 1, 1, &
      1,-1, 0,     1, 1, 1, &
      0,-1,-1,    -1,-1, 1, &
      1, 0, 1,    -1,-1, 1, &
     -1, 1, 0,    -1,-1, 1, &
      0,-1, 1,     1,-1,-1, &
     -1, 0,-1,     1,-1,-1, &
      1, 1, 0,     1,-1,-1, &
      0, 1, 1,    -1, 1,-1, &
      1, 0,-1,    -1, 1,-1, &
     -1,-1, 0,    -1, 1,-1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_fcc_Nslip])                                                 !< Slip system <110>{111} directions. Sorted according to Eisenlohr & Hantcherli

 real(pReal), dimension(3+3,LATTICE_fcc_Ntwin), parameter, private :: &
   LATTICE_fcc_systemTwin = reshape(real( [&
     -2, 1, 1,     1, 1, 1, &
      1,-2, 1,     1, 1, 1, &
      1, 1,-2,     1, 1, 1, &
      2,-1, 1,    -1,-1, 1, &
     -1, 2, 1,    -1,-1, 1, &
     -1,-1,-2,    -1,-1, 1, &
     -2,-1,-1,     1,-1,-1, &
      1, 2,-1,     1,-1,-1, &
      1,-1, 2,     1,-1,-1, &
      2, 1,-1,    -1, 1,-1, &
     -1,-2,-1,    -1, 1,-1, &
     -1, 1, 2,    -1, 1,-1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_fcc_Ntwin])                                                 !< Twin system <112>{111} directions. Sorted according to Eisenlohr & Hantcherli

 real(pReal), dimension(3+3,LATTICE_fcc_Ntrans), parameter, private :: &
   LATTICE_fccTohex_systemTrans = reshape(real( [&
     -2, 1, 1,     1, 1, 1, &
      1,-2, 1,     1, 1, 1, &
      1, 1,-2,     1, 1, 1, &
      2,-1, 1,    -1,-1, 1, &
     -1, 2, 1,    -1,-1, 1, &
     -1,-1,-2,    -1,-1, 1, &
     -2,-1,-1,     1,-1,-1, &
      1, 2,-1,     1,-1,-1, &
      1,-1, 2,     1,-1,-1, &
      2, 1,-1,    -1, 1,-1, &
     -1,-2,-1,    -1, 1,-1, &
     -1, 1, 2,    -1, 1,-1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_fcc_Ntrans])

 real(pReal), dimension(LATTICE_fcc_Ntwin), parameter, private :: &
   LATTICE_fcc_shearTwin = 0.5_pReal*sqrt(2.0_pReal)                                                !< Twin system <112>{111} ??? Sorted according to Eisenlohr & Hantcherli

 integer(pInt), dimension(2_pInt,LATTICE_fcc_Ntwin), parameter, public :: &
   LATTICE_fcc_twinNucleationSlipPair = reshape(int( [&
     2,3, &
     1,3, &
     1,2, &
     5,6, &
     4,6, &
     4,5, &
     8,9, &
     7,9, &
     7,8, &
     11,12, &
     10,12, &
     10,11 &
     ],pInt),[2_pInt,LATTICE_fcc_Ntwin])

 integer(pInt), dimension(LATTICE_fcc_Nslip,lattice_fcc_Nslip), parameter, public :: &
   LATTICE_fcc_interactionSlipSlip = reshape(int( [&
     1,2,2,4,6,5,3,5,5,4,5,6, &  ! ---> slip
     2,1,2,6,4,5,5,4,6,5,3,5, &  ! |
     2,2,1,5,5,3,5,6,4,6,5,4, &  ! |
     4,6,5,1,2,2,4,5,6,3,5,5, &  ! v slip
     6,4,5,2,1,2,5,3,5,5,4,6, &
     5,5,3,2,2,1,6,5,4,5,6,4, &
     3,5,5,4,5,6,1,2,2,4,6,5, &
     5,4,6,5,3,5,2,1,2,6,4,5, &
     5,6,4,6,5,4,2,2,1,5,5,3, &
     4,5,6,3,5,5,4,6,5,1,2,2, &
     5,3,5,5,4,6,6,4,5,2,1,2, &
     6,5,4,5,6,4,5,5,3,2,2,1  &
     ],pInt),[LATTICE_fcc_Nslip,LATTICE_fcc_Nslip],order=[2,1])                                     !< Slip--slip interaction types for fcc
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: coplanar interaction
                                                                                                    !< 3: collinear interaction
                                                                                                    !< 4: Hirth locks
                                                                                                    !< 5: glissile junctions
                                                                                                    !< 6: Lomer locks
 integer(pInt), dimension(LATTICE_fcc_Nslip,LATTICE_fcc_Ntwin), parameter, public :: &
   LATTICE_fcc_interactionSlipTwin = reshape(int( [&
     1,1,1,3,3,3,2,2,2,3,3,3, & ! ---> twin
     1,1,1,3,3,3,3,3,3,2,2,2, & ! |
     1,1,1,2,2,2,3,3,3,3,3,3, & ! |
     3,3,3,1,1,1,3,3,3,2,2,2, & ! v slip
     3,3,3,1,1,1,2,2,2,3,3,3, &
     2,2,2,1,1,1,3,3,3,3,3,3, &
     2,2,2,3,3,3,1,1,1,3,3,3, &
     3,3,3,2,2,2,1,1,1,3,3,3, &
     3,3,3,3,3,3,1,1,1,2,2,2, &
     3,3,3,2,2,2,3,3,3,1,1,1, &
     2,2,2,3,3,3,3,3,3,1,1,1, &
     3,3,3,3,3,3,2,2,2,1,1,1  &
     ],pInt),[LATTICE_fcc_Nslip,LATTICE_fcc_Ntwin],order=[2,1])                                     !< Slip--twin interaction types for fcc
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
 integer(pInt), dimension(LATTICE_fcc_Ntwin,LATTICE_fcc_Nslip), parameter, public :: &
   LATTICE_fcc_interactionTwinSlip = 1_pInt                                                         !< Twin--Slip interaction types for fcc

 integer(pInt), dimension(LATTICE_fcc_Ntwin,LATTICE_fcc_Ntwin), parameter,public :: &
   LATTICE_fcc_interactionTwinTwin = reshape(int( [&
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! ---> twin
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     2,2,2,1,1,1,2,2,2,2,2,2, &  ! v twin
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1  &
     ],pInt),[lattice_fcc_Ntwin,lattice_fcc_Ntwin],order=[2,1])                                     !< Twin--twin interaction types for fcc

 integer(pInt), dimension(LATTICE_fcc_Nslip,LATTICE_fcc_Ntrans), parameter, public :: &
   LATTICE_fccTohex_interactionSlipTrans = reshape(int( [&
     1,1,1,3,3,3,2,2,2,3,3,3, & ! ---> trans
     1,1,1,3,3,3,3,3,3,2,2,2, & ! |
     1,1,1,2,2,2,3,3,3,3,3,3, & ! |
     3,3,3,1,1,1,3,3,3,2,2,2, & ! v slip
     3,3,3,1,1,1,2,2,2,3,3,3, &
     2,2,2,1,1,1,3,3,3,3,3,3, &
     2,2,2,3,3,3,1,1,1,3,3,3, &
     3,3,3,2,2,2,1,1,1,3,3,3, &
     3,3,3,3,3,3,1,1,1,2,2,2, &
     3,3,3,2,2,2,3,3,3,1,1,1, &
     2,2,2,3,3,3,3,3,3,1,1,1, &
     3,3,3,3,3,3,2,2,2,1,1,1  &
     ],pInt),[LATTICE_fcc_Nslip,LATTICE_fcc_Ntrans],order=[2,1])                                    !< Slip--trans interaction types for fcc

 integer(pInt), dimension(LATTICE_fcc_Ntrans,LATTICE_fcc_Nslip), parameter, public :: &
   LATTICE_fccTohex_interactionTransSlip = 1_pInt                                                   !< Trans--Slip interaction types for fcc

 integer(pInt), dimension(LATTICE_fcc_Ntrans,LATTICE_fcc_Ntrans), parameter,public :: &
   LATTICE_fccTohex_interactionTransTrans = reshape(int( [&
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! ---> trans
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     2,2,2,1,1,1,2,2,2,2,2,2, &  ! v trans
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1  &
     ],pInt),[LATTICE_fcc_Ntrans,LATTICE_fcc_Ntrans],order=[2,1])                                   !< Trans--trans interaction types for fcc
 
 real(pReal), dimension(LATTICE_fcc_Ntrans), parameter, private :: &
   LATTICE_fccTohex_shearTrans = sqrt(2.0_pReal)/4.0_pReal
    
 real(pReal), dimension(4,LATTICE_fcc_Ntrans), parameter, private :: &
   LATTICE_fccTobcc_systemTrans = reshape([&
     0.0, 1.0, 0.0,     10.26, &                                                                    ! Pitsch OR (Ma & Hartmaier 2014, Table 3)
     0.0, 1.0, 0.0,    -10.26, &
     0.0, 0.0, 1.0,     10.26, &
     0.0, 0.0, 1.0,    -10.26, &
     1.0, 0.0, 0.0,     10.26, &
     1.0, 0.0, 0.0,    -10.26, &
     0.0, 0.0, 1.0,     10.26, &
     0.0, 0.0, 1.0,    -10.26, &
     1.0, 0.0, 0.0,     10.26, &
     1.0, 0.0, 0.0,    -10.26, &
     0.0, 1.0, 0.0,     10.26, &
     0.0, 1.0, 0.0,    -10.26  &
     ],[ 4_pInt,LATTICE_fcc_Ntrans])

 integer(pInt), dimension(9,LATTICE_fcc_Ntrans), parameter, private :: &
   LATTICE_fccTobcc_bainVariant = reshape(int( [&
     1, 0, 0, 0, 1, 0, 0, 0, 1, &                                                                   ! Pitsch OR (Ma & Hartmaier 2014, Table 3)
     1, 0, 0, 0, 1, 0, 0, 0, 1, &
     1, 0, 0, 0, 1, 0, 0, 0, 1, &
     1, 0, 0, 0, 1, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 0, 1, 1, 0, 0, 0, 1, 0, &
     0, 0, 1, 1, 0, 0, 0, 1, 0, &
     0, 0, 1, 1, 0, 0, 0, 1, 0, &
     0, 0, 1, 1, 0, 0, 0, 1, 0  & 
     ],pInt),[ 9_pInt, LATTICE_fcc_Ntrans])

 real(pReal), dimension(4,LATTICE_fcc_Ntrans), parameter, private :: &
   LATTICE_fccTobcc_bainRot = reshape([&
     1.0, 0.0, 0.0,     45.0, &                                                                    ! Rotate fcc austensite to bain variant
     1.0, 0.0, 0.0,     45.0, &
     1.0, 0.0, 0.0,     45.0, &
     1.0, 0.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 0.0, 1.0,     45.0, &
     0.0, 0.0, 1.0,     45.0, &
     0.0, 0.0, 1.0,     45.0, &
     0.0, 0.0, 1.0,     45.0  &
     ],[ 4_pInt,LATTICE_fcc_Ntrans])

 real(pReal), dimension(LATTICE_fcc_Ntrans,LATTICE_fcc_Ntrans), parameter, private :: &            ! Matrix for projection of shear from slip system to fault-band (twin) systems
   LATTICE_fccTobcc_projectionTrans = reshape(real([&                                                        ! For ns = nt = nr
     0, 1,-1,  0, 0, 0,  0, 0, 0,  0, 0, 0, &                                                                    
    -1, 0, 1,  0, 0, 0,  0, 0, 0,  0, 0, 0, &
     1,-1, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0, &
     0, 0, 0,  0, 1,-1,  0, 0, 0,  0, 0, 0, &
     0, 0, 0, -1, 0, 1,  0, 0, 0,  0, 0, 0, &
     0, 0, 0,  1,-1, 0,  0, 0, 0,  0, 0, 0, &
     0, 0, 0,  0, 0, 0,  0, 1,-1,  0, 0, 0, &
     0, 0, 0,  0, 0, 0, -1, 0, 1,  0, 0, 0, &
     0, 0, 0,  0, 0, 0,  1,-1, 0,  0, 0, 0, &
     0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 1,-1, &
     0, 0, 0,  0, 0, 0,  0, 0, 0, -1, 0, 1, &
     0, 0, 0,  0, 0, 0,  0, 0, 0,  1,-1, 0  &
     ],pReal),[LATTICE_fcc_Ntrans,LATTICE_fcc_Ntrans],order=[2,1])

 real(pReal), parameter, private  :: &
   LATTICE_fccTobcc_projectionTransFactor = sqrt(3.0_pReal/4.0_pReal)

 real(pReal), parameter, public  :: &
   LATTICE_fccTobcc_shearCritTrans = 0.0224

 integer(pInt), dimension(2_pInt,LATTICE_fcc_Ntrans), parameter, public :: &
   LATTICE_fccTobcc_transNucleationTwinPair = reshape(int( [&
     4,  7, &
     1, 10, &
     1,  4, &
     7, 10, &
     2,  8, &
     5, 11, &
     8, 11, &
     2,  5, &
     6, 12, &
     3,  9, &
     3, 12, &
     6,  9  &
     ],pInt),[2_pInt,LATTICE_fcc_Ntrans])

 real(pReal), dimension(3+3,LATTICE_fcc_Ncleavage), parameter, private :: &
   LATTICE_fcc_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1, &
      0, 1,-1,     1, 1, 1, &
      0,-1,-1,    -1,-1, 1, &
     -1, 0,-1,     1,-1,-1, &
      0, 1, 1,    -1, 1,-1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_fcc_Ncleavage])

!--------------------------------------------------------------------------------------------------
! bcc
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, public :: &
   LATTICE_bcc_NslipSystem = int([ 12, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], pInt)                      !< total # of slip systems per family for bcc
   
 integer(pInt), dimension(LATTICE_maxNtwinFamily), parameter, public :: &
   LATTICE_bcc_NtwinSystem = int([ 12, 0, 0, 0], pInt)                                              !< total # of twin systems per family for bcc

 integer(pInt), dimension(LATTICE_maxNtransFamily), parameter, public :: &
   LATTICE_bcc_NtransSystem = int([0,0],pInt)                                                       !< total # of transformation systems per family for bcc

 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, public :: &
   LATTICE_bcc_NcleavageSystem = int([3,6,0],pInt)                                                    !< total # of cleavage systems per family for bcc
   
 integer(pInt), parameter, private  :: &
   LATTICE_bcc_Nslip = 24_pInt, & ! sum(lattice_bcc_NslipSystem), &                                 !< total # of slip systems for bcc
   LATTICE_bcc_Ntwin = 12_pInt, & ! sum(lattice_bcc_NtwinSystem)                                    !< total # of twin systems for bcc
   LATTICE_bcc_NnonSchmid = 6_pInt, &                                                               !< # of non-Schmid contributions for bcc. 6 known non schmid contributions for BCC (A. Koester, A. Ma, A. Hartmaier 2012)
   LATTICE_bcc_Ntrans  = 0_pInt, &                                                                  !< total # of transformations for bcc
   LATTICE_bcc_Ncleavage  = 9_pInt                                                                  !< total # of cleavage systems for bcc

   
 real(pReal), dimension(3+3,LATTICE_bcc_Nslip), parameter, private :: &
   LATTICE_bcc_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
    ! Slip system <111>{110} 
      1,-1, 1,     0, 1, 1, &
     -1,-1, 1,     0, 1, 1, &
      1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     1, 0, 1, &
     -1,-1, 1,     1, 0, 1, &
      1, 1, 1,    -1, 0, 1, &
      1,-1, 1,    -1, 0, 1, &
     -1, 1, 1,     1, 1, 0, &
     -1, 1,-1,     1, 1, 0, &
      1, 1, 1,    -1, 1, 0, &
      1, 1,-1,    -1, 1, 0, &
    ! Slip system <111>{112}
     -1, 1, 1,     2, 1, 1, &
      1, 1, 1,    -2, 1, 1, &
      1, 1,-1,     2,-1, 1, &
      1,-1, 1,     2, 1,-1, &
      1,-1, 1,     1, 2, 1, &
      1, 1,-1,    -1, 2, 1, &
      1, 1, 1,     1,-2, 1, &
     -1, 1, 1,     1, 2,-1, &
      1, 1,-1,     1, 1, 2, &
      1,-1, 1,    -1, 1, 2, &
     -1, 1, 1,     1,-1, 2, &
      1, 1, 1,     1, 1,-2  &
   !  Slip system <111>{123}
   !  1, 1,-1,     1, 2, 3, &
   !  1,-1, 1,    -1, 2, 3, &
   ! -1, 1, 1,     1,-2, 3, &
   !  1, 1, 1,     1, 2,-3, &
   !  1,-1, 1,     1, 3, 2, &
   !  1, 1,-1,    -1, 3, 2, &
   !  1, 1, 1,     1,-3, 2, &
   ! -1, 1, 1,     1, 3,-2, &
   !  1, 1,-1,     2, 1, 3, &
   !  1,-1, 1,    -2, 1, 3, &
   ! -1, 1, 1,     2,-1, 3, &
   !  1, 1, 1,     2, 1,-3, &
   !  1,-1, 1,     2, 3, 1, &
   !  1, 1,-1,    -2, 3, 1, &
   !  1, 1, 1,     2,-3, 1, &
   ! -1, 1, 1,     2, 3,-1, &
   ! -1, 1, 1,     3, 1, 2, &
   !  1, 1, 1,    -3, 1, 2, &
   !  1, 1,-1,     3,-1, 2, &
   !  1,-1, 1,     3, 1,-2, &
   ! -1, 1, 1,     3, 2, 1, &
   !  1, 1, 1,    -3, 2, 1, &
   !  1, 1,-1,     3,-2, 1, &
   !  1,-1, 1,     3, 2,-1  &
     ],pReal),[ 3_pInt + 3_pInt ,LATTICE_bcc_Nslip])

 real(pReal), dimension(3+3,LATTICE_bcc_Ntwin), parameter, private :: &
   LATTICE_bcc_systemTwin = reshape(real([&
    ! Twin system <111>{112}
     -1, 1, 1,     2, 1, 1, &
      1, 1, 1,    -2, 1, 1, &
      1, 1,-1,     2,-1, 1, &
      1,-1, 1,     2, 1,-1, &
      1,-1, 1,     1, 2, 1, &
      1, 1,-1,    -1, 2, 1, &
      1, 1, 1,     1,-2, 1, &
     -1, 1, 1,     1, 2,-1, &
      1, 1,-1,     1, 1, 2, &
      1,-1, 1,    -1, 1, 2, &
     -1, 1, 1,     1,-1, 2, &
      1, 1, 1,     1, 1,-2  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_bcc_Ntwin])

 real(pReal), dimension(LATTICE_bcc_Ntwin), parameter, private :: &
   LATTICE_bcc_shearTwin = 0.5_pReal*sqrt(2.0_pReal)

 integer(pInt), dimension(LATTICE_bcc_Nslip,LATTICE_bcc_Nslip), parameter, public :: &
   LATTICE_bcc_interactionSlipSlip = reshape(int( [&
     1,2,6,6,5,4,4,3,4,3,5,4, 6,6,4,3,3,4,6,6,4,3,6,6, &  ! ---> slip  
     2,1,6,6,4,3,5,4,5,4,4,3, 6,6,3,4,4,3,6,6,3,4,6,6, &  ! |
     6,6,1,2,4,5,3,4,4,5,3,4, 4,3,6,6,6,6,3,4,6,6,4,3, &  ! |
     6,6,2,1,3,4,4,5,3,4,4,5, 3,4,6,6,6,6,4,3,6,6,3,4, &  ! v slip
     5,4,4,3,1,2,6,6,3,4,5,4, 3,6,4,6,6,4,6,3,4,6,3,6, &
     4,3,5,4,2,1,6,6,4,5,4,3, 4,6,3,6,6,3,6,4,3,6,4,6, &
     4,5,3,4,6,6,1,2,5,4,3,4, 6,3,6,4,4,6,3,6,6,4,6,3, &
     3,4,4,5,6,6,2,1,4,3,4,5, 6,4,6,3,3,6,4,6,6,3,6,4, &
     4,5,4,3,3,4,5,4,1,2,6,6, 3,6,6,4,4,6,6,3,6,4,3,6, &
     3,4,5,4,4,5,4,3,2,1,6,6, 4,6,6,3,3,6,6,4,6,3,4,6, &
     5,4,3,4,5,4,3,4,6,6,1,2, 6,3,4,6,6,4,3,6,4,6,6,3, &
     4,3,4,5,4,3,4,5,6,6,2,1, 6,4,3,6,6,3,4,6,3,6,6,4, &
     !
     6,6,4,3,3,4,6,6,3,4,6,6, 1,5,6,6,5,6,6,3,5,6,3,6, &
     6,6,3,4,6,6,3,4,6,6,3,4, 5,1,6,6,6,5,3,6,6,5,6,3, &
     4,3,6,6,4,3,6,6,6,6,4,3, 6,6,1,5,6,3,5,6,3,6,5,6, &
     3,4,6,6,6,6,4,3,4,3,6,6, 6,6,5,1,3,6,6,5,6,3,6,5, &
     3,4,6,6,6,6,4,3,4,3,6,6, 5,6,6,3,1,6,5,6,5,3,6,6, &
     4,3,6,6,4,3,6,6,6,6,4,3, 6,5,3,6,6,1,6,5,3,5,6,6, &
     6,6,3,4,6,6,3,4,6,6,3,4, 6,3,5,6,5,6,1,6,6,6,5,3, &
     6,6,4,3,3,4,6,6,3,4,6,6, 3,6,6,5,6,5,6,1,6,6,3,5, &
     4,3,6,6,4,3,6,6,6,6,4,3, 5,6,3,6,5,3,6,6,1,6,6,5, &
     3,4,6,6,6,6,4,3,4,3,6,6, 6,5,6,3,3,5,6,6,6,1,5,6, &
     6,6,4,3,3,4,6,6,3,4,6,6, 3,6,5,6,6,6,5,3,6,5,1,6, &
     6,6,3,4,6,6,3,4,6,6,3,4, 6,3,6,5,6,6,3,5,5,6,6,1  &
     ],pInt),[lattice_bcc_Nslip,lattice_bcc_Nslip],order=[2,1])                                     !< Slip--slip interaction types for bcc from Queyreau et al. Int J Plast 25 (2009) 361–377
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: coplanar interaction
                                                                                                    !< 3: collinear interaction
                                                                                                    !< 4: mixed-asymmetrical junction
                                                                                                    !< 5: mixed-symmetrical junction
                                                                                                    !< 6: edge junction
 integer(pInt), dimension(LATTICE_bcc_Nslip,LATTICE_bcc_Ntwin), parameter, public :: &
  LATTICE_bcc_interactionSlipTwin = reshape(int( [&
     3,3,3,2,2,3,3,3,3,2,3,3, &  ! ---> twin
     3,3,2,3,3,2,3,3,2,3,3,3, &  ! |
     3,2,3,3,3,3,2,3,3,3,3,2, &  ! |
     2,3,3,3,3,3,3,2,3,3,2,3, &  ! v slip 
     2,3,3,3,3,3,3,2,3,3,2,3, &
     3,3,2,3,3,2,3,3,2,3,3,3, &
     3,2,3,3,3,3,2,3,3,3,3,2, &
     3,3,3,2,2,3,3,3,3,2,3,3, &
     2,3,3,3,3,3,3,2,3,3,2,3, &
     3,3,3,2,2,3,3,3,3,2,3,3, &
     3,2,3,3,3,3,2,3,3,3,3,2, &
     3,3,2,3,3,2,3,3,2,3,3,3, &
     !
     1,3,3,3,3,3,3,2,3,3,2,3, &
     3,1,3,3,3,3,2,3,3,3,3,2, &
     3,3,1,3,3,2,3,3,2,3,3,3, &
     3,3,3,1,2,3,3,3,3,2,3,3, &
     3,3,3,2,1,3,3,3,3,2,3,3, &
     3,3,2,3,3,1,3,3,2,3,3,3, &
     3,2,3,3,3,3,1,3,3,3,3,2, &
     2,3,3,3,3,3,3,1,3,3,2,3, &
     3,3,2,3,3,2,3,3,1,3,3,3, &
     3,3,3,2,2,3,3,3,3,1,3,3, &
     2,3,3,3,3,3,3,2,3,3,1,3, &
     3,2,3,3,3,3,2,3,3,3,3,1  &
     ],pInt),[LATTICE_bcc_Nslip,LATTICE_bcc_Ntwin],order=[2,1])                                     !< Slip--twin interaction types for bcc
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
 integer(pInt), dimension(LATTICE_bcc_Ntwin,LATTICE_bcc_Nslip), parameter, public :: &
   LATTICE_bcc_interactionTwinSlip = 1_pInt                                                         !< Twin--slip interaction types for bcc @todo not implemented yet

 integer(pInt), dimension(LATTICE_bcc_Ntwin,LATTICE_bcc_Ntwin), parameter, public :: &
   LATTICE_bcc_interactionTwinTwin = reshape(int( [&
     1,3,3,3,3,3,3,2,3,3,2,3, &  ! ---> twin
     3,1,3,3,3,3,2,3,3,3,3,2, &  ! |
     3,3,1,3,3,2,3,3,2,3,3,3, &  ! |
     3,3,3,1,2,3,3,3,3,2,3,3, &  ! v twin
     3,3,3,2,1,3,3,3,3,2,3,3, &
     3,3,2,3,3,1,3,3,2,3,3,3, &
     3,2,3,3,3,3,1,3,3,3,3,2, &
     2,3,3,3,3,3,3,1,3,3,2,3, &
     3,3,2,3,3,2,3,3,1,3,3,3, &
     3,3,3,2,2,3,3,3,3,1,3,3, &
     2,3,3,3,3,3,3,2,3,3,1,3, &
     3,2,3,3,3,3,2,3,3,3,3,1  &
     ],pInt),[LATTICE_bcc_Ntwin,LATTICE_bcc_Ntwin],order=[2,1])                                     !< Twin--twin interaction types for bcc
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: collinear interaction
                                                                                                    !< 3: other interaction
 real(pReal), dimension(3+3,LATTICE_bcc_Ncleavage), parameter, private :: &
   LATTICE_bcc_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1, &
      1,-1, 1,     0, 1, 1, &
      1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     1, 0, 1, &
      1, 1, 1,    -1, 0, 1, &
     -1, 1, 1,     1, 1, 0, &
      1, 1, 1,    -1, 1, 0  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_bcc_Ncleavage])

!--------------------------------------------------------------------------------------------------
! hex
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, public :: &
   lattice_hex_NslipSystem = int([ 3, 3, 3, 6, 12, 6, 0, 0, 0, 0, 0, 0, 0],pInt)                             !< # of slip systems per family for hex
   
 integer(pInt), dimension(LATTICE_maxNtwinFamily), parameter, public :: &
   lattice_hex_NtwinSystem = int([ 6, 6, 6, 6],pInt)                                                !< # of slip systems per family for hex 

 integer(pInt), dimension(LATTICE_maxNtransFamily), parameter, public :: &
   LATTICE_hex_NtransSystem = int([0,0],pInt)                                                       !< total # of transformation systems per family for hex
   
 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, public :: &
   LATTICE_hex_NcleavageSystem = int([3,0,0],pInt)                                                    !< total # of cleavage systems per family for hex
   
 integer(pInt), parameter , private :: &
   LATTICE_hex_Nslip = 33_pInt, & ! sum(lattice_hex_NslipSystem),                                   !< total # of slip systems for hex 
   LATTICE_hex_Ntwin = 24_pInt, & ! sum(lattice_hex_NtwinSystem)                                    !< total # of twin systems for hex
   LATTICE_hex_NnonSchmid = 0_pInt, &                                                               !< # of non-Schmid contributions for hex
   LATTICE_hex_Ntrans  = 0_pInt, &                                                                  !< total # of transformations for hex
   LATTICE_hex_Ncleavage  = 3_pInt                                                                  !< total # of transformations for hex

 real(pReal), dimension(4+4,LATTICE_hex_Nslip), parameter, private :: &
   LATTICE_hex_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
    ! Basal systems <11.0>{00.1} (independent of c/a-ratio, Bravais notation (4 coordinate base))
      2, -1, -1,  0,     0,  0,  0,  1, &
     -1,  2, -1,  0,     0,  0,  0,  1, &
     -1, -1,  2,  0,     0,  0,  0,  1, &
    ! 1st type prismatic systems <11.0>{10.0}  (independent of c/a-ratio)
      2, -1, -1,  0,     0,  1, -1,  0, &
     -1,  2, -1,  0,    -1,  0,  1,  0, &
     -1, -1,  2,  0,     1, -1,  0,  0, &
    ! 2nd type prismatic systems <10.0>{11.0} -- a slip; plane normals independent of c/a-ratio
      0,  1,  -1, 0,     2, -1, -1,  0, &
     -1,  0,  1,  0,    -1,  2, -1,  0, &
      1, -1,  0,  0,    -1, -1,  2,  0,  &
    ! 1st type 1st order pyramidal systems <11.0>{-11.1} -- plane normals depend on the c/a-ratio
      2, -1, -1,  0,     0,  1, -1,  1, &
     -1,  2, -1,  0,    -1,  0,  1,  1, &
     -1, -1,  2,  0,     1, -1,  0,  1, &
      1,  1, -2,  0,    -1,  1,  0,  1, &
     -2,  1,  1,  0,     0, -1,  1,  1, &
      1, -2,  1,  0,     1,  0, -1,  1, &
    ! pyramidal system: c+a slip <11.3>{-10.1} -- plane normals depend on the c/a-ratio
      2, -1, -1,  3,    -1,  1,  0,  1, &
      1, -2,  1,  3,    -1,  1,  0,  1, &
     -1, -1,  2,  3,     1,  0, -1,  1, &
     -2,  1,  1,  3,     1,  0, -1,  1, &
     -1,  2, -1,  3,     0, -1,  1,  1, &
      1,  1, -2,  3,     0, -1,  1,  1, &
     -2,  1,  1,  3,     1, -1,  0,  1, &
     -1,  2, -1,  3,     1, -1,  0,  1, &
      1,  1, -2,  3,    -1,  0,  1,  1, &
      2, -1, -1,  3,    -1,  0,  1,  1, &
      1, -2,  1,  3,     0,  1, -1,  1, &
     -1, -1,  2,  3,     0,  1, -1,  1, &
    ! pyramidal system: c+a slip <11.3>{-1-1.2} -- as for hexagonal ice (Castelnau et al. 1996, similar to twin system found below) 
      2, -1, -1,  3,    -2,  1,  1,  2, & ! sorted according to similar twin system
     -1,  2, -1,  3,     1, -2,  1,  2, & ! <11.3>{-1-1.2} shear = 2((c/a)^2-2)/(3 c/a)
     -1, -1,  2,  3,     1,  1, -2,  2, &
     -2,  1,  1,  3,     2, -1, -1,  2, &
      1, -2,  1,  3,    -1,  2, -1,  2, &
      1,  1, -2,  3,    -1, -1,  2,  2  &
     ],pReal),[ 4_pInt + 4_pInt,LATTICE_hex_Nslip])                                                 !< slip systems for hex sorted by A. Alankar & P. Eisenlohr

 real(pReal), dimension(4+4,LATTICE_hex_Ntwin), parameter, private :: &
   LATTICE_hex_systemTwin =  reshape(real([&
    ! Compression or Tension =f(twinning shear=f(c/a)) for each metal ! (according to Yoo 1981)
      1, -1,  0,  1,    -1,  1,  0,  2, & ! <-10.1>{10.2} shear = (3-(c/a)^2)/(sqrt(3) c/a)
     -1,  0,  1,  1,     1,  0, -1,  2, &
      0,  1, -1,  1,     0, -1,  1,  2, &
     -1,  1,  0,  1,     1, -1,  0,  2, &
      1,  0, -1,  1,    -1,  0,  1,  2, &
      0, -1,  1,  1,     0,  1, -1,  2, &
!
      2, -1, -1,  6,    -2,  1,  1,  1, & ! <11.6>{-1-1.1} shear = 1/(c/a)
     -1,  2, -1,  6,     1, -2,  1,  1, &
     -1, -1,  2,  6,     1,  1, -2,  1, &
     -2,  1,  1,  6,     2, -1, -1,  1, &
      1, -2,  1,  6,    -1,  2, -1,  1, &
      1,  1, -2,  6,    -1, -1,  2,  1, &
!
     -1,  1,  0, -2,    -1,  1,  0,  1, & !! <10.-2>{10.1} shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
      1,  0, -1, -2,     1,  0, -1,  1, &
      0, -1,  1, -2,     0, -1,  1,  1, &
      1, -1,  0, -2,     1, -1,  0,  1, &
     -1,  0,  1, -2,    -1,  0,  1,  1, &
      0,  1, -1, -2,     0,  1, -1,  1, &
!
      2, -1, -1, -3,     2, -1, -1,  2, & ! <11.-3>{11.2} shear = 2((c/a)^2-2)/(3 c/a)
     -1,  2, -1, -3,    -1,  2, -1,  2, &
     -1, -1,  2, -3,    -1, -1,  2,  2, &
     -2,  1,  1, -3,    -2,  1,  1,  2, &
      1, -2,  1, -3,     1, -2,  1,  2, &
      1,  1, -2, -3,     1,  1, -2,  2  &
     ],pReal),[ 4_pInt + 4_pInt ,LATTICE_hex_Ntwin])                                                !< twin systems for hex, order follows Prof. Tom Bieler's scheme; but numbering in data was restarted from 1 

 integer(pInt), dimension(LATTICE_hex_Ntwin), parameter, private :: &
   LATTICE_hex_shearTwin = reshape(int( [&   ! indicator to formula further below
     1, &  ! <-10.1>{10.2}
     1, &
     1, &
     1, &
     1, &
     1, &
     2, &  ! <11.6>{-1-1.1}
     2, &
     2, &
     2, &
     2, &
     2, &
     3, &  ! <10.-2>{10.1}
     3, &
     3, &
     3, &
     3, &
     3, &
     4, &  ! <11.-3>{11.2}
     4, &
     4, &
     4, &
     4, &
     4  &
     ],pInt),[LATTICE_hex_Ntwin])
   
 integer(pInt), dimension(LATTICE_hex_Nslip,LATTICE_hex_Nslip), parameter, public :: &
   LATTICE_hex_interactionSlipSlip = reshape(int( [&
      1, 2, 2,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, &  ! ---> slip
      2, 1, 2,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, &  ! |
      2, 2, 1,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, &  ! |
    !                                                                                                                   v slip
      6, 6, 6,   4, 5, 5,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, &
      6, 6, 6,   5, 4, 5,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, &
      6, 6, 6,   5, 5, 4,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, &
    !
     12,12,12,  11,11,11,   9,10,10,  15,15,15,15,15,15,  23,23,23,23,23,23,23,23,23,23,23,23,  33,33,33,33,33,33, &
     12,12,12,  11,11,11,  10, 9,10,  15,15,15,15,15,15,  23,23,23,23,23,23,23,23,23,23,23,23,  33,33,33,33,33,33, &
     12,12,12,  11,11,11,  10,10, 9,  15,15,15,15,15,15,  23,23,23,23,23,23,23,23,23,23,23,23,  33,33,33,33,33,33, &
    !
     20,20,20,  19,19,19,  18,18,18,  16,17,17,17,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,16,17,17,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,16,17,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,17,16,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,17,17,16,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,17,17,17,16,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
    !
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  25,26,26,26,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,25,26,26,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,25,26,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,25,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,25,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,25,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,25,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,25,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,25,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,26,25,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,26,26,25,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,26,26,26,25,  35,35,35,35,35,35, &
    !
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  36,37,37,37,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,36,37,37,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,36,37,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,37,36,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,37,37,36,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,37,37,37,36  &  
    !
           ],pInt),[LATTICE_hex_Nslip,LATTICE_hex_Nslip],order=[2,1])                                     !< Slip--slip interaction types for hex (onion peel naming scheme)
  
 integer(pInt), dimension(LATTICE_hex_Nslip,LATTICE_hex_Ntwin), parameter, public :: &
   LATTICE_hex_interactionSlipTwin = reshape(int( [&
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! --> twin
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
    !                                                                                   v
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & ! slip
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & 
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, &
    !
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
    !
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
    !
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
    !
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24  &
    !
     ],pInt),[LATTICE_hex_Nslip,LATTICE_hex_Ntwin],order=[2,1])                                     !< Slip--twin interaction types for hex (isotropic, 24 in total) 

 integer(pInt), dimension(LATTICE_hex_Ntwin,LATTICE_hex_Nslip), parameter, public :: &
   LATTICE_hex_interactionTwinSlip = reshape(int( [&
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! --> slip
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! v
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! twin
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, &
    !
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
    !
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
    !
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24  &
     ],pInt),[LATTICE_hex_Ntwin,LATTICE_hex_Nslip],order=[2,1])                                     !< Twin--twin interaction types for hex (isotropic, 20 in total)

 integer(pInt), dimension(LATTICE_hex_Ntwin,LATTICE_hex_Ntwin), parameter, public :: &
   LATTICE_hex_interactionTwinTwin = reshape(int( [&
      1, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! ---> twin
      2, 1, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! |
      2, 2, 1, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! |
      2, 2, 2, 1, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! v twin
      2, 2, 2, 2, 1, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &
      2, 2, 2, 2, 2, 1,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &
    !
      6, 6, 6, 6, 6, 6,   4, 5, 5, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 4, 5, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 4, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 5, 4, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 5, 5, 4, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 5, 5, 5, 4,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
    !
     12,12,12,12,12,12,  11,11,11,11,11,11,   9,10,10,10,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10, 9,10,10,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10, 9,10,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10, 9,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10,10, 9,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10,10,10, 9,  15,15,15,15,15,15, &
    !
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  16,17,17,17,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,16,17,17,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,16,17,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,16,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,17,16,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,17,17,16  &
     ],pInt),[lattice_hex_Ntwin,lattice_hex_Ntwin],order=[2,1])                                     !< Twin--slip interaction types for hex (isotropic, 16 in total)

 real(pReal), dimension(4+4,LATTICE_hex_Ncleavage), parameter, private :: &
   LATTICE_hex_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      2,-1,-1, 0,     0, 0, 0, 1, &
      0, 0, 0, 1,     2,-1,-1, 0, &
      0, 0, 0, 1,     0, 1,-1, 0  &
     ],pReal),[ 4_pInt + 4_pInt,LATTICE_hex_Ncleavage])



!--------------------------------------------------------------------------------------------------
! bct
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, public :: &
   LATTICE_bct_NslipSystem = int([2, 2, 2, 4, 2, 4, 2, 2, 4, 8, 4, 8, 8 ],pInt)                     !< # of slip systems per family for bct (Sn) Bieler J. Electr Mater 2009

 integer(pInt), dimension(LATTICE_maxNtwinFamily), parameter, public :: &
   LATTICE_bct_NtwinSystem = int([0, 0, 0, 0], pInt)                                                !< total # of twin systems per family for bct-example

 integer(pInt), dimension(LATTICE_maxNtransFamily), parameter, public :: &
   LATTICE_bct_NtransSystem = int([0,0],pInt)                                                       !< total # of transformation systems per family for bct

 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, public :: &
   LATTICE_bct_NcleavageSystem = int([0,0,0],pInt)                                                  !< total # of cleavage systems per family for bct
 
   
 integer(pInt), parameter , private :: &
   LATTICE_bct_Nslip = 52_pInt, & ! sum(lattice_bct_NslipSystem),                                   !< total # of slip systems for bct 
   LATTICE_bct_Ntwin = 0_pInt, &  ! sum(lattice_bcc_NtwinSystem)                                    !< total # of twin systems for bct
   LATTICE_bct_NnonSchmid = 0_pInt, &                                                               !< # of non-Schmid contributions for bct
   LATTICE_bct_Ntrans  = 0_pInt, &                                                                  !< total # of transformations for bct
   LATTICE_bct_Ncleavage  = 0_pInt                                                                  !< total # of transformations for bct
 
 real(pReal), dimension(3+3,LATTICE_bct_Nslip), parameter, private :: &
   LATTICE_bct_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
    ! Slip family 1 {100)<001] (Bravais notation {hkl)<uvw] for bct c/a = 0.5456)
      0, 0, 1,      1, 0, 0, &
      0, 0, 1,      0, 1, 0, &
    ! Slip family 2 {110)<001]
      0, 0, 1,      1, 1, 0, &
      0, 0, 1,     -1, 1, 0, &
    ! slip family 3 {100)<010] 
      0,  1, 0,     1, 0, 0, &
      1,  0, 0,     0, 1, 0, &
    ! Slip family 4 {110)<1-11]/2
      1,-1, 1,      1, 1, 0, &
      1,-1,-1,      1, 1, 0, &
     -1,-1,-1,     -1, 1, 0, &
     -1,-1, 1,     -1, 1, 0, &
    ! Slip family 5 {110)<1-10]
      1, -1, 0,     1, 1, 0, &
      1,  1, 0,     1,-1, 0, &
    ! Slip family 6 {100)<011] 
      0, 1, 1,      1, 0, 0, &
      0,-1, 1,      1, 0, 0, &
     -1, 0, 1,      0, 1, 0, &
      1, 0, 1,      0, 1, 0, &  
    ! Slip family 7 {001)<010] 
      0, 1, 0,      0, 0, 1, &
      1, 0, 0,      0, 0, 1, &
    ! Slip family 8 {001)<110] 
      1, 1, 0,      0, 0, 1, &
     -1, 1, 0,      0, 0, 1, &  
    ! Slip family 9 {011)<01-1] 
      0, 1,-1,      0, 1, 1, &
      0,-1,-1,      0,-1, 1, &
     -1, 0,-1,     -1, 0, 1, &
      1, 0,-1,      1, 0, 1, &  
    ! Slip family 10 {011)<1-11]/2 
      1,-1, 1,      0, 1, 1, &
      1, 1,-1,      0, 1, 1, &
      1, 1, 1,      0, 1,-1, &
     -1, 1, 1,      0, 1,-1, &  
      1,-1,-1,      1, 0, 1, &
     -1,-1, 1,      1, 0, 1, &
      1, 1, 1,      1, 0,-1, &
      1,-1, 1,      1, 0,-1,  &  
    ! Slip family 11 {011)<100] 
      1, 0, 0,      0, 1, 1, &
      1, 0, 0,      0, 1,-1, &
      0, 1, 0,      1, 0, 1, &
      0, 1, 0,      1, 0,-1, &  
    ! Slip family 12 {211)<01-1] 
      0, 1,-1,      2, 1, 1, &
      0,-1,-1,      2,-1, 1, &
      1, 0,-1,      1, 2, 1, &
     -1, 0,-1,     -1, 2, 1, &  
      0, 1,-1,     -2, 1, 1, &
      0,-1,-1,     -2,-1, 1, &
     -1, 0,-1,     -1,-2, 1, &
      1, 0,-1,      1,-2, 1, &  
    ! Slip family 13 {211)<-111]/2 
     -1, 1, 1,      2, 1, 1, &
     -1,-1, 1,      2,-1, 1, &
      1,-1, 1,      1, 2, 1, &
     -1,-1, 1,     -1, 2, 1, &  
      1, 1, 1,     -2, 1, 1, &
      1,-1, 1,     -2,-1, 1, &
     -1, 1, 1,     -1,-2, 1, &
      1, 1, 1,      1,-2, 1  &  
      ],pReal),[ 3_pInt + 3_pInt,LATTICE_bct_Nslip])                                                 !< slip systems for bct sorted by Bieler

 integer(pInt), dimension(LATTICE_bct_Nslip,LATTICE_bct_Nslip), parameter, public :: &
   LATTICE_bct_interactionSlipSlip = reshape(int( [&
     1,  2,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91,  111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, &
     2,  1,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91,  111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, &
    !
     6,  6,   4,  5,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92,  112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, &
     6,  6,   5,  4,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92,  112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, &
    ! 
    12, 12,  11, 11,   9, 10,  15, 15, 15, 15,  23, 23,  33, 33, 33, 33,  45, 45,  59, 59,  75, 75, 75, 75,  93, 93, 93, 93, 93, 93, 93, 93,  113, 113, 113, 113, 135,135,135,135,135,135,135,135, 159,159,159,159,159,159,159,159, &
    12, 12,  11, 11,  10,  9,  15, 15, 15, 15,  23, 23,  33, 33, 33, 33,  45, 45,  59, 59,  75, 75, 75, 75,  93, 93, 93, 93, 93, 93, 93, 93,  113, 113, 113, 113, 135,135,135,135,135,135,135,135, 159,159,159,159,159,159,159,159, &
    !
    20, 20,  19, 19,  18, 18,  16, 17, 17, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    20, 20,  19, 19,  18, 18,  17, 16, 17, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    20, 20,  19, 19,  18, 18,  17, 17, 16, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    20, 20,  19, 19,  18, 18,  17, 17, 17, 16,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    !
    30, 30,  29, 29,  28, 28,  27, 27, 27, 27,  25, 26,  35, 35, 35, 35,  47, 47,  61, 61,  77, 77, 77, 77,  95, 95, 95, 95, 95, 95, 95, 95,  115, 115, 115, 115, 137,137,137,137,137,137,137,137, 161,161,161,161,161,161,161,161, &
    30, 30,  29, 29,  28, 28,  27, 27, 27, 27,  26, 25,  35, 35, 35, 35,  47, 47,  61, 61,  77, 77, 77, 77,  95, 95, 95, 95, 95, 95, 95, 95,  115, 115, 115, 115, 137,137,137,137,137,137,137,137, 161,161,161,161,161,161,161,161, &
    !
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  36, 37, 37, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 36, 37, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 37, 36, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 37, 37, 36,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    !
    56, 56,  55, 55,  54, 54,  53, 53, 53, 53,  52, 52,  51, 51, 51, 51,  49, 50,  63, 63,  79, 79, 79, 79,  97, 97, 97, 97, 97, 97, 97, 97,  117, 117, 117, 117, 139,139,139,139,139,139,139,139, 163,163,163,163,163,163,163,163, &
    56, 56,  55, 55,  54, 54,  53, 53, 53, 53,  52, 52,  51, 51, 51, 51,  50, 49,  63, 63,  79, 79, 79, 79,  97, 97, 97, 97, 97, 97, 97, 97,  117, 117, 117, 117, 139,139,139,139,139,139,139,139, 163,163,163,163,163,163,163,163, &
    !
    72, 72,  71, 71,  70, 70,  69, 69, 69, 69,  68, 68,  67, 67, 67, 67,  66, 66,  64, 65,  80, 80, 80, 80,  98, 98, 98, 98, 98, 98, 98, 98,  118, 118, 118, 118, 140,140,140,140,140,140,140,140, 164,164,164,164,164,164,164,164, &
    72, 72,  71, 71,  70, 70,  69, 69, 69, 69,  68, 68,  67, 67, 67, 67,  66, 66,  65, 64,  80, 80, 80, 80,  98, 98, 98, 98, 98, 98, 98, 98,  118, 118, 118, 118, 140,140,140,140,140,140,140,140, 164,164,164,164,164,164,164,164, &
    !
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  81, 82, 82, 82,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 81, 82, 82,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 82, 81, 82,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 82, 82, 81,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    !
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 100,101,101,101,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,100,101,101,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,100,101,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,100,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,100,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,100,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,101,100,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,101,101,100,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
    !
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 122, 122, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 121, 122, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 122, 121, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 122, 122, 121, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
    !
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 144,145,145,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,144,145,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,144,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,144,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,144,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,145,144,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,145,145,144,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,145,145,145,144, 168,168,168,168,168,168,168,168, &
    !
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,169,170,170,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,169,170,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,170,169,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,170,170,169,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,169,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,169,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,170,169  &

 ],pInt),[lattice_bct_Nslip,lattice_bct_Nslip],order=[2,1])                                  

!--------------------------------------------------------------------------------------------------
! isotropic
 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, public :: &
   LATTICE_iso_NcleavageSystem = int([3,0,0],pInt)                                                    !< total # of cleavage systems per family for isotropic
   
 integer(pInt), parameter, private  :: &
   LATTICE_iso_Ncleavage  = 3_pInt                                                                  !< total # of cleavage systems for bcc

 real(pReal), dimension(3+3,LATTICE_iso_Ncleavage), parameter, private :: &
   LATTICE_iso_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_iso_Ncleavage])

!--------------------------------------------------------------------------------------------------
! orthorhombic
 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, public :: &
   LATTICE_ortho_NcleavageSystem = int([1,1,1],pInt)                                                    !< total # of cleavage systems per family for orthotropic
   
 integer(pInt), parameter, private  :: &
   LATTICE_ortho_Ncleavage  = 3_pInt                                                                  !< total # of cleavage systems for bcc

 real(pReal), dimension(3+3,LATTICE_ortho_Ncleavage), parameter, private :: &
   LATTICE_ortho_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_ortho_Ncleavage])

 real(pReal),                              dimension(:,:,:),   allocatable, public, protected :: &
   lattice_C66,   lattice_trans_C66
 real(pReal),                              dimension(:,:,:,:,:),   allocatable, public, protected :: &
   lattice_C3333, lattice_trans_C3333
 real(pReal),                              dimension(:),   allocatable, public, protected :: &
   lattice_mu, &
   lattice_nu, &
   lattice_trans_mu, &
   lattice_trans_nu
 real(pReal),                              dimension(:,:,:),   allocatable, public, protected :: &
   lattice_thermalConductivity33, &
   lattice_thermalExpansion33, &
   lattice_damageDiffusion33, &
   lattice_vacancyfluxDiffusion33, &
   lattice_vacancyfluxMobility33, &
   lattice_porosityDiffusion33, &
   lattice_hydrogenfluxDiffusion33, &
   lattice_hydrogenfluxMobility33
 real(pReal),                              dimension(:),       allocatable, public, protected :: &
   lattice_damageMobility, &
   lattice_porosityMobility, &
   lattice_massDensity, &
   lattice_specificHeat, &
   lattice_vacancyFormationEnergy, &
   lattice_vacancySurfaceEnergy, &
   lattice_vacancyVol, &
   lattice_hydrogenFormationEnergy, &
   lattice_hydrogenSurfaceEnergy, &
   lattice_hydrogenVol, &
   lattice_referenceTemperature, &
   lattice_equilibriumVacancyConcentration, &
   lattice_equilibriumHydrogenConcentration
 enum, bind(c)
   enumerator :: LATTICE_undefined_ID, &
                 LATTICE_iso_ID, &
                 LATTICE_fcc_ID, &
                 LATTICE_bcc_ID, &
                 LATTICE_hex_ID, &
                 LATTICE_bct_ID, &
                 LATTICE_ort_ID
 end enum
 integer(kind(LATTICE_undefined_ID)),        dimension(:),       allocatable, public, protected :: &
   lattice_structure, trans_lattice_structure



integer(pInt), dimension(2), parameter, private :: &
   lattice_NsymOperations = [24_pInt,12_pInt]       

real(pReal), dimension(4,36), parameter, private :: &
  lattice_symOperations = reshape([&
     1.0_pReal,                 0.0_pReal,                 0.0_pReal,                 0.0_pReal, &                      ! cubic symmetry operations
     0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal,  0.7071067811865476_pReal, &       !     2-fold symmetry
     0.0_pReal,                 0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal, &
     0.0_pReal,                 0.7071067811865476_pReal,  0.7071067811865476_pReal,  0.0_pReal, &
     0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal, -0.7071067811865476_pReal, &
     0.0_pReal,                -0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal, &
     0.0_pReal,                 0.7071067811865476_pReal, -0.7071067811865476_pReal,  0.0_pReal, &
     0.5_pReal,                 0.5_pReal,                 0.5_pReal,                 0.5_pReal, &                      !     3-fold symmetry
    -0.5_pReal,                 0.5_pReal,                 0.5_pReal,                 0.5_pReal, &
     0.5_pReal,                -0.5_pReal,                 0.5_pReal,                 0.5_pReal, &
    -0.5_pReal,                -0.5_pReal,                 0.5_pReal,                 0.5_pReal, &
     0.5_pReal,                 0.5_pReal,                -0.5_pReal,                 0.5_pReal, &
    -0.5_pReal,                 0.5_pReal,                -0.5_pReal,                 0.5_pReal, &
     0.5_pReal,                 0.5_pReal,                 0.5_pReal,                -0.5_pReal, &
    -0.5_pReal,                 0.5_pReal,                 0.5_pReal,                -0.5_pReal, &
     0.7071067811865476_pReal,  0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal, &                      !     4-fold symmetry
     0.0_pReal,                 1.0_pReal,                 0.0_pReal,                 0.0_pReal, &
    -0.7071067811865476_pReal,  0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal, &
     0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal,  0.0_pReal, &
     0.0_pReal,                 0.0_pReal,                 1.0_pReal,                 0.0_pReal, &
    -0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal,  0.0_pReal, &
     0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal, &
     0.0_pReal,                 0.0_pReal,                 0.0_pReal,                 1.0_pReal, &
    -0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal, &
!
     1.0_pReal,                 0.0_pReal,                 0.0_pReal,                 0.0_pReal, &                      ! hexagonal symmetry operations
     0.0_pReal,                 1.0_pReal,                 0.0_pReal,                 0.0_pReal, &                      !     2-fold symmetry
     0.0_pReal,                 0.0_pReal,                 1.0_pReal,                 0.0_pReal, &
     0.0_pReal,                 0.5_pReal,                 0.866025403784439_pReal,   0.0_pReal, &
     0.0_pReal,                -0.5_pReal,                 0.866025403784439_pReal,   0.0_pReal, &
     0.0_pReal,                 0.866025403784439_pReal,   0.5_pReal,                 0.0_pReal, &
     0.0_pReal,                -0.866025403784439_pReal,   0.5_pReal,                 0.0_pReal, &
     0.866025403784439_pReal,   0.0_pReal,                 0.0_pReal,                 0.5_pReal, &                      !     6-fold symmetry
    -0.866025403784439_pReal,   0.0_pReal,                 0.0_pReal,                 0.5_pReal, &
     0.5_pReal,                 0.0_pReal,                 0.0_pReal,                 0.866025403784439_pReal, &
    -0.5_pReal,                 0.0_pReal,                 0.0_pReal,                 0.866025403784439_pReal, &
     0.0_pReal,                 0.0_pReal,                 0.0_pReal,                 1.0_pReal &
     ],[4,36])  !< Symmetry operations as quaternions 24 for cubic, 12 for hexagonal = 36

 ! use this later on to substitute the matrix above
 !   if self.lattice == 'cubic':
 !     symQuats =  [
 !                   [ 1.0,0.0,0.0,0.0 ],
 !                   [ 0.0,1.0,0.0,0.0 ],
 !                   [ 0.0,0.0,1.0,0.0 ],
 !                   [ 0.0,0.0,0.0,1.0 ],
 !                   [ 0.0, 0.0, 0.5*math.sqrt(2), 0.5*math.sqrt(2) ],
 !                   [ 0.0, 0.0, 0.5*math.sqrt(2),-0.5*math.sqrt(2) ],
 !                   [ 0.0, 0.5*math.sqrt(2), 0.0, 0.5*math.sqrt(2) ],
 !                   [ 0.0, 0.5*math.sqrt(2), 0.0,-0.5*math.sqrt(2) ],
 !                   [ 0.0, 0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0 ],
 !                   [ 0.0,-0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0 ],
 !                   [ 0.5, 0.5, 0.5, 0.5 ],
 !                   [-0.5, 0.5, 0.5, 0.5 ],
 !                   [-0.5, 0.5, 0.5,-0.5 ],
 !                   [-0.5, 0.5,-0.5, 0.5 ],
 !                   [-0.5,-0.5, 0.5, 0.5 ],
 !                   [-0.5,-0.5, 0.5,-0.5 ],
 !                   [-0.5,-0.5,-0.5, 0.5 ],
 !                   [-0.5, 0.5,-0.5,-0.5 ],
 !                   [-0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
 !                   [ 0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
 !                   [-0.5*math.sqrt(2), 0.0, 0.5*math.sqrt(2), 0.0 ],
 !                   [-0.5*math.sqrt(2), 0.0,-0.5*math.sqrt(2), 0.0 ],
 !                   [-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0, 0.0 ],
 !                   [-0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0, 0.0 ],
 !                 ]
 !   elif self.lattice == 'hexagonal':
 !     symQuats =  [
 !                   [ 1.0,0.0,0.0,0.0 ],
 !                   [ 0.0,1.0,0.0,0.0 ],
 !                   [ 0.0,0.0,1.0,0.0 ],
 !                   [ 0.0,0.0,0.0,1.0 ],
 !                   [-0.5*math.sqrt(3), 0.0, 0.0, 0.5 ],
 !                   [-0.5*math.sqrt(3), 0.0, 0.0,-0.5 ],
 !                   [ 0.0, 0.5*math.sqrt(3), 0.5, 0.0 ],
 !                   [ 0.0,-0.5*math.sqrt(3), 0.5, 0.0 ],
 !                   [ 0.0, 0.5,-0.5*math.sqrt(3), 0.0 ],
 !                   [ 0.0,-0.5,-0.5*math.sqrt(3), 0.0 ],
 !                   [ 0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
 !                   [-0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
 !                 ]
 !   elif self.lattice == 'tetragonal':
 !     symQuats =  [
 !                   [ 1.0,0.0,0.0,0.0 ],
 !                   [ 0.0,1.0,0.0,0.0 ],
 !                   [ 0.0,0.0,1.0,0.0 ],
 !                   [ 0.0,0.0,0.0,1.0 ],
 !                   [ 0.0, 0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0 ],
 !                   [ 0.0,-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0 ],
 !                   [ 0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
 !                   [-0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
 !                 ]
 !   elif self.lattice == 'orthorhombic':
 !     symQuats =  [
 !                   [ 1.0,0.0,0.0,0.0 ],
 !                   [ 0.0,1.0,0.0,0.0 ],
 !                   [ 0.0,0.0,1.0,0.0 ],
 !                   [ 0.0,0.0,0.0,1.0 ],
 !                 ]
 !   else:
 !     symQuats =  [
 !                   [ 1.0,0.0,0.0,0.0 ],
 !                 ]

 public :: &
  lattice_init, &
  lattice_qDisorientation, &
  LATTICE_fcc_ID, &
  LATTICE_bcc_ID, &
  LATTICE_bct_ID, &
  LATTICE_hex_ID

contains

!--------------------------------------------------------------------------------------------------
!> @brief Module initialization
!--------------------------------------------------------------------------------------------------
subroutine lattice_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_open_file,&
   IO_open_jobFile_stat, &
   IO_countSections, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF, &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue
 use material, only: &
   material_configfile, &
   material_localFileExt, &
   material_partPhase
 use debug, only: &
   debug_level, &
   debug_lattice, &
   debug_levelBasic
 use numerics, only: &
   worldrank
   
 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: Nphases
 character(len=65536) :: &
   tag  = '', &
   line = ''
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: section = 0_pInt,i
 real(pReal),  dimension(:), allocatable :: &
   CoverA, &                                                                                        !!!!!!< c/a ratio for low symmetry type lattice
   CoverA_trans, &                                                                                  !< c/a ratio for transformed hex type lattice      
   a_fcc, &                                                                                         !< lattice parameter a for fcc austenite
   a_bcc                                                                                            !< lattice paramater a for bcc martensite

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)') ' <<<+-  lattice init  -+>>>'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

!--------------------------------------------------------------------------------------------------
! consistency checks

 if (LATTICE_maxNslip /= maxval([LATTICE_fcc_Nslip,LATTICE_bcc_Nslip,LATTICE_hex_Nslip,LATTICE_bct_Nslip])) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_maxNslip')
 if (LATTICE_maxNtwin /= maxval([LATTICE_fcc_Ntwin,LATTICE_bcc_Ntwin,LATTICE_hex_Ntwin])) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_maxNtwin')
 if (LATTICE_maxNtrans /= maxval([LATTICE_fcc_Ntrans,LATTICE_bcc_Ntrans,LATTICE_hex_Ntrans])) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_maxNtrans')
 if (LATTICE_maxNnonSchmid /= maxval([lattice_fcc_NnonSchmid,lattice_bcc_NnonSchmid,&
   lattice_hex_NnonSchmid])) call IO_error(0_pInt,ext_msg = 'LATTICE_maxNnonSchmid')

 if (LATTICE_fcc_Nslip /= sum(lattice_fcc_NslipSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_fcc_Nslip')
 if (LATTICE_bcc_Nslip /= sum(lattice_bcc_NslipSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bcc_Nslip')
 if (LATTICE_hex_Nslip /= sum(lattice_hex_NslipSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_hex_Nslip')
 if (LATTICE_bct_Nslip /= sum(lattice_bct_NslipSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bct_Nslip')

 if (LATTICE_fcc_Ntwin /= sum(lattice_fcc_NtwinSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_fcc_Ntwin')
 if (LATTICE_bcc_Ntwin /= sum(lattice_bcc_NtwinSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bcc_Ntwin')
 if (LATTICE_hex_Ntwin /= sum(lattice_hex_NtwinSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_hex_Ntwin')
 if (LATTICE_bct_Ntwin /= sum(lattice_bct_NtwinSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bct_Ntwin')

 if (LATTICE_fcc_Ntrans /= sum(lattice_fcc_NtransSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_fcc_Ntrans')
 if (LATTICE_bcc_Ntrans /= sum(lattice_bcc_NtransSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bcc_Ntrans')
 if (LATTICE_hex_Ntrans /= sum(lattice_hex_NtransSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_hex_Ntrans')
 if (LATTICE_bct_Ntrans /= sum(lattice_bct_NtransSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bct_Ntrans')

 if (LATTICE_fcc_Ncleavage /= sum(lattice_fcc_NcleavageSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_fcc_Ncleavage')
 if (LATTICE_bcc_Ncleavage /= sum(lattice_bcc_NcleavageSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bcc_Ncleavage')
 if (LATTICE_hex_Ncleavage /= sum(lattice_hex_NcleavageSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_hex_Ncleavage')
 if (LATTICE_bct_Ncleavage /= sum(lattice_bct_NcleavageSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_bct_Ncleavage')
 if (LATTICE_iso_Ncleavage /= sum(lattice_iso_NcleavageSystem)) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_iso_Ncleavage')

 if (LATTICE_maxNinteraction /= max(&
   maxval(lattice_fcc_interactionSlipSlip), &
   maxval(lattice_bcc_interactionSlipSlip), &
   maxval(lattice_hex_interactionSlipSlip), &
   maxval(lattice_bct_interactionSlipSlip), &
   !
   maxval(lattice_fcc_interactionSlipTwin), &
   maxval(lattice_bcc_interactionSlipTwin), &
   maxval(lattice_hex_interactionSlipTwin), &
!  maxval(lattice_bct_interactionSlipTwin), &
   !
   maxval(lattice_fcc_interactionTwinSlip), &
   maxval(lattice_bcc_interactionTwinSlip), &
   maxval(lattice_hex_interactionTwinSlip), &
!  maxval(lattice_bct_interactionTwinSlip), &
   !
   maxval(lattice_fcc_interactionTwinTwin), &
   maxval(lattice_bcc_interactionTwinTwin), &
   maxval(lattice_hex_interactionTwinTwin))) &
!  maxval(lattice_bct_interactionTwinTwin))) &
   call IO_error(0_pInt,ext_msg = 'LATTICE_maxNinteraction')

!--------------------------------------------------------------------------------------------------
! read from material configuration file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 Nphases = IO_countSections(FILEUNIT,material_partPhase)

 if(Nphases<1_pInt) &
   call IO_error(160_pInt,Nphases, ext_msg='No phases found')

 if (iand(debug_level(debug_lattice),debug_levelBasic) /= 0_pInt) then
   write(6,'(a16,1x,i5)')   ' # phases:',Nphases
 endif
 
 allocate(lattice_structure(Nphases),source = LATTICE_undefined_ID)
 allocate(trans_lattice_structure(Nphases),source = LATTICE_undefined_ID)
 allocate(lattice_C66(6,6,Nphases),  source=0.0_pReal)
 allocate(lattice_C3333(3,3,3,3,Nphases),  source=0.0_pReal)
 allocate(lattice_trans_C66(6,6,Nphases),  source=0.0_pReal)
 allocate(lattice_trans_C3333(3,3,3,3,Nphases),  source=0.0_pReal)
 allocate(lattice_thermalConductivity33  (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_thermalExpansion33     (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_damageDiffusion33      (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_vacancyfluxDiffusion33 (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_vacancyfluxMobility33  (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_PorosityDiffusion33    (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_hydrogenfluxDiffusion33(3,3,Nphases), source=0.0_pReal)
 allocate(lattice_hydrogenfluxMobility33 (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_damageMobility         (    Nphases), source=0.0_pReal)
 allocate(lattice_PorosityMobility       (    Nphases), source=0.0_pReal)
 allocate(lattice_massDensity            (    Nphases), source=0.0_pReal)
 allocate(lattice_specificHeat           (    Nphases), source=0.0_pReal)
 allocate(lattice_vacancyFormationEnergy (    Nphases), source=0.0_pReal)
 allocate(lattice_vacancySurfaceEnergy   (    Nphases), source=0.0_pReal)
 allocate(lattice_vacancyVol             (    Nphases), source=0.0_pReal)
 allocate(lattice_hydrogenFormationEnergy(    Nphases), source=0.0_pReal)
 allocate(lattice_hydrogenSurfaceEnergy  (    Nphases), source=0.0_pReal)
 allocate(lattice_hydrogenVol            (    Nphases), source=0.0_pReal)
 allocate(lattice_referenceTemperature   (    Nphases), source=300.0_pReal)
 allocate(lattice_equilibriumVacancyConcentration(Nphases), source=0.0_pReal)
 allocate(lattice_equilibriumHydrogenConcentration(Nphases),source=0.0_pReal)

 allocate(lattice_mu(Nphases),       source=0.0_pReal)
 allocate(lattice_nu(Nphases),       source=0.0_pReal)
 allocate(lattice_trans_mu(Nphases), source=0.0_pReal)
 allocate(lattice_trans_nu(Nphases), source=0.0_pReal)

 allocate(lattice_NnonSchmid(Nphases), source=0_pInt)
 allocate(lattice_Sslip(3,3,1+2*lattice_maxNnonSchmid,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_Sslip_v(6,1+2*lattice_maxNnonSchmid,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_Scleavage(3,3,3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_Scleavage_v(6,3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_sd(3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_st(3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_sn(3,lattice_maxNslip,Nphases),source=0.0_pReal)

 allocate(lattice_Qtwin(3,3,lattice_maxNtwin,Nphases),source=0.0_pReal)
 allocate(lattice_Stwin(3,3,lattice_maxNtwin,Nphases),source=0.0_pReal)
 allocate(lattice_Stwin_v(6,lattice_maxNtwin,Nphases),source=0.0_pReal)
 allocate(lattice_td(3,lattice_maxNtwin,Nphases),source=0.0_pReal)
 allocate(lattice_tt(3,lattice_maxNtwin,Nphases),source=0.0_pReal)
 allocate(lattice_tn(3,lattice_maxNtwin,Nphases),source=0.0_pReal)

 allocate(lattice_shearTwin(lattice_maxNtwin,Nphases),source=0.0_pReal)
 allocate(lattice_shearTrans(lattice_maxNtrans,Nphases),source=0.0_pReal)

 allocate(lattice_Qtrans(3,3,lattice_maxNtrans,Nphases),source=0.0_pReal)
 allocate(lattice_Strans(3,3,lattice_maxNtrans,Nphases),source=0.0_pReal)
 allocate(lattice_Strans_v(6,lattice_maxNtrans,Nphases),source=0.0_pReal)
 allocate(lattice_projectionTrans(lattice_maxNtrans,lattice_maxNtrans,Nphases),source=0.0_pReal)

 allocate(lattice_NslipSystem(lattice_maxNslipFamily,Nphases),source=0_pInt)
 allocate(lattice_NtwinSystem(lattice_maxNtwinFamily,Nphases),source=0_pInt)
 allocate(lattice_NtransSystem(lattice_maxNtransFamily,Nphases),source=0_pInt)
 allocate(lattice_NcleavageSystem(lattice_maxNcleavageFamily,Nphases),source=0_pInt)

 allocate(lattice_interactionSlipSlip(lattice_maxNslip,lattice_maxNslip,Nphases),source=0_pInt)     ! other:me
 allocate(lattice_interactionSlipTwin(lattice_maxNslip,lattice_maxNtwin,Nphases),source=0_pInt)     ! other:me
 allocate(lattice_interactionTwinSlip(lattice_maxNtwin,lattice_maxNslip,Nphases),source=0_pInt)     ! other:me
 allocate(lattice_interactionTwinTwin(lattice_maxNtwin,lattice_maxNtwin,Nphases),source=0_pInt)     ! other:me
 allocate(lattice_interactionSlipTrans(lattice_maxNslip,lattice_maxNtrans,Nphases),source=0_pInt)   ! other:me
 allocate(lattice_interactionTransSlip(lattice_maxNtrans,lattice_maxNslip,Nphases),source=0_pInt)   ! other:me
 allocate(lattice_interactionTransTrans(lattice_maxNtrans,lattice_maxNtrans,Nphases),source=0_pInt) ! other:me

 allocate(CoverA(Nphases),source=0.0_pReal)
 allocate(CoverA_trans(Nphases),source=0.0_pReal)
 allocate(a_fcc(Nphases),source=0.0_pReal)
 allocate(a_bcc(Nphases),source=0.0_pReal)

 rewind(fileUnit)
 line        = ''                                                                                   ! to have it initialized
 section     = 0_pInt                                                                               !  - " -
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <Phase>
   line = IO_read(fileUnit)
 enddo

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
   endif
   if (section > 0_pInt) then
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('lattice_structure')
         select case(trim(IO_lc(IO_stringValue(line,chunkPos,2_pInt))))
           case('iso','isotropic')
             lattice_structure(section) = LATTICE_iso_ID
           case('fcc')
             lattice_structure(section) = LATTICE_fcc_ID
           case('bcc')
             lattice_structure(section) = LATTICE_bcc_ID
           case('hex','hexagonal')
             lattice_structure(section) = LATTICE_hex_ID
           case('bct')
             lattice_structure(section) = LATTICE_bct_ID
           case('ort','orthorhombic')
             lattice_structure(section) = LATTICE_ort_ID
           case default
             call IO_error(130_pInt,ext_msg=trim(IO_lc(IO_stringValue(line,chunkPos,2_pInt))))
         end select
     case('trans_lattice_structure')
       select case(trim(IO_lc(IO_stringValue(line,chunkPos,2_pInt))))
         case('bcc')
           trans_lattice_structure(section) = LATTICE_bcc_ID
         case('hex','hexagonal','hcp')
           trans_lattice_structure(section) = LATTICE_hex_ID
       end select
     case ('c11')
       lattice_C66(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c12')
       lattice_C66(1,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c13')
       lattice_C66(1,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c22')
       lattice_C66(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c23')
       lattice_C66(2,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c33')
       lattice_C66(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c44')
       lattice_C66(4,4,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c55')
       lattice_C66(5,5,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c66')
       lattice_C66(6,6,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c11_trans')
       lattice_trans_C66(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c12_trans')
       lattice_trans_C66(1,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c13_trans')
       lattice_trans_C66(1,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c22_trans')
       lattice_trans_C66(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c23_trans')
       lattice_trans_C66(2,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c33_trans')
       lattice_trans_C66(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c44_trans')
       lattice_trans_C66(4,4,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c55_trans')
       lattice_trans_C66(5,5,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c66_trans')
       lattice_trans_C66(6,6,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('covera_ratio','c/a_ratio','c/a')
       CoverA(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('c/a_trans','c/a_martensite','c/a_mart')
       CoverA_trans(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('a_fcc')
       a_fcc(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('a_bcc')
       a_bcc(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('thermal_conductivity11')
       lattice_thermalConductivity33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('thermal_conductivity22')
       lattice_thermalConductivity33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('thermal_conductivity33')
       lattice_thermalConductivity33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('thermal_expansion11')
       lattice_thermalExpansion33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('thermal_expansion22')
       lattice_thermalExpansion33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('thermal_expansion33')
       lattice_thermalExpansion33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('specific_heat')
       lattice_specificHeat(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyformationenergy')
       lattice_vacancyFormationEnergy(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancysurfaceenergy')
       lattice_vacancySurfaceEnergy(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyvolume')
       lattice_vacancyVol(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenformationenergy')
       lattice_hydrogenFormationEnergy(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogensurfaceenergy')
       lattice_hydrogenSurfaceEnergy(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenvolume')
       lattice_hydrogenVol(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('mass_density')
       lattice_massDensity(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('reference_temperature')
       lattice_referenceTemperature(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('damage_diffusion11')
       lattice_DamageDiffusion33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('damage_diffusion22')
       lattice_DamageDiffusion33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('damage_diffusion33')
       lattice_DamageDiffusion33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('damage_mobility')
       lattice_DamageMobility(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyflux_diffusion11')
       lattice_vacancyfluxDiffusion33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyflux_diffusion22')
       lattice_vacancyfluxDiffusion33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyflux_diffusion33')
       lattice_vacancyfluxDiffusion33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyflux_mobility11')
       lattice_vacancyfluxMobility33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyflux_mobility22')
       lattice_vacancyfluxMobility33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancyflux_mobility33')
       lattice_vacancyfluxMobility33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('porosity_diffusion11')
       lattice_PorosityDiffusion33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('porosity_diffusion22')
       lattice_PorosityDiffusion33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('porosity_diffusion33')
       lattice_PorosityDiffusion33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('porosity_mobility')
       lattice_PorosityMobility(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenflux_diffusion11')
       lattice_hydrogenfluxDiffusion33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenflux_diffusion22')
       lattice_hydrogenfluxDiffusion33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenflux_diffusion33')
       lattice_hydrogenfluxDiffusion33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenflux_mobility11')
       lattice_hydrogenfluxMobility33(1,1,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenflux_mobility22')
       lattice_hydrogenfluxMobility33(2,2,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogenflux_mobility33')
       lattice_hydrogenfluxMobility33(3,3,section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('vacancy_eqcv')
       lattice_equilibriumVacancyConcentration(section) = IO_floatValue(line,chunkPos,2_pInt)
     case ('hydrogen_eqch')
       lattice_equilibriumHydrogenConcentration(section) = IO_floatValue(line,chunkPos,2_pInt)
     end select
   endif
 enddo
 
 do i = 1_pInt,Nphases
   if ((CoverA(i) < 1.0_pReal .or. CoverA(i) > 2.0_pReal) &
       .and. lattice_structure(i) == LATTICE_hex_ID) call IO_error(131_pInt,el=i)                        ! checking physical significance of c/a
   if ((CoverA(i) > 2.0_pReal) &
       .and. lattice_structure(i) == LATTICE_bct_ID) call IO_error(131_pInt,el=i)                        ! checking physical significance of c/a
   call lattice_initializeStructure(i, CoverA(i), CoverA_trans(i), a_fcc(i), a_bcc(i))
 enddo

 deallocate(CoverA,CoverA_trans,a_fcc,a_bcc)

end subroutine lattice_init


!--------------------------------------------------------------------------------------------------
!> @brief   Calculation of Schmid matrices, etc.
!--------------------------------------------------------------------------------------------------
subroutine lattice_initializeStructure(myPhase,CoverA,CoverA_trans,a_fcc,a_bcc)
 use prec, only: &
  tol_math_check
 use math, only: &
   math_crossproduct, &
   math_tensorproduct33, &
   math_mul33x33, &
   math_mul33x3, &
   math_transpose33, &
   math_trace33, &
   math_symmetric33, &
   math_Mandel33to6, &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
   math_axisAngleToR, &
   INRAD, &
   MATH_I3
 use IO, only: &
   IO_error, &
   IO_warning
 
 implicit none
 integer(pInt), intent(in) :: myPhase
 real(pReal), intent(in) :: &
   CoverA, &
   CoverA_trans, &
   a_fcc, &
   a_bcc

 real(pReal), dimension(3) :: &
   sdU, snU, &
   np,  nn
 real(pReal), dimension(3,3) :: &
   sstr, sdtr, sttr
 real(pReal), dimension(3,lattice_maxNslip) :: &
   sd,  sn
 real(pReal), dimension(3,3,2,lattice_maxNnonSchmid,lattice_maxNslip) :: &
   sns
 real(pReal), dimension(3,lattice_maxNtwin) :: &
   td, tn
 real(pReal), dimension(lattice_maxNtwin) :: &
   ts
 real(pReal), dimension(lattice_maxNtrans) :: &
   trs
 real(pReal), dimension(3,lattice_maxNtrans) :: &
   xtr, ytr, ztr
 real(pReal), dimension(3,3,lattice_maxNtrans) :: &
   Rtr, Utr, Btr, Qtr, Str
 real(pReal), dimension(3,lattice_maxNcleavage) :: &
   cd,  cn, ct
 integer(pInt) :: &
  i,j, &
  myNslip = 0_pInt, myNtwin = 0_pInt, myNtrans = 0_pInt, myNcleavage = 0_pInt
 real(pReal) :: c11bar, c12bar, c13bar, c14bar, c33bar, c44bar, A, B

 lattice_C66(1:6,1:6,myPhase) = lattice_symmetrizeC66(lattice_structure(myPhase),&
                                                      lattice_C66(1:6,1:6,myPhase))
 
 lattice_mu(myPhase) = 0.2_pReal *(  lattice_C66(1,1,myPhase) &
                                   - lattice_C66(1,2,myPhase) &
                                   + 3.0_pReal*lattice_C66(4,4,myPhase))                             ! (C11iso-C12iso)/2 with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
 lattice_nu(myPhase) = (  lattice_C66(1,1,myPhase) &
                        + 4.0_pReal*lattice_C66(1,2,myPhase) &
                        - 2.0_pReal*lattice_C66(4,4,myPhase)) &
                                                             /(  4.0_pReal*lattice_C66(1,1,myPhase) &
                                                               + 6.0_pReal*lattice_C66(1,2,myPhase) &
                                                               + 2.0_pReal*lattice_C66(4,4,myPhase))! C12iso/(C11iso+C12iso) with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
 lattice_C3333(1:3,1:3,1:3,1:3,myPhase) = math_Voigt66to3333(lattice_C66(1:6,1:6,myPhase))          ! Literature data is Voigt
 lattice_C66(1:6,1:6,myPhase) = math_Mandel3333to66(lattice_C3333(1:3,1:3,1:3,1:3,myPhase))         ! DAMASK uses Mandel
 do i = 1_pInt, 6_pInt
   if (abs(lattice_C66(i,i,myPhase))<tol_math_check) &
     call IO_error(135_pInt,el=i,ip=myPhase,ext_msg='matrix diagonal "el"ement of phase "ip"')
 enddo

! Elasticity matrices for transformed phase
 select case(lattice_structure(myPhase))
   case (LATTICE_fcc_ID)
     select case(trans_lattice_structure(myPhase))
       case (LATTICE_bcc_ID)
         lattice_trans_C66(1:6,1:6,myPhase) = lattice_C66(1:6,1:6,myPhase)
         lattice_trans_mu(myPhase) = lattice_mu(myPhase)
         lattice_trans_nu(myPhase) = lattice_nu(myPhase)
         lattice_trans_C3333(1:3,1:3,1:3,1:3,myPhase) = lattice_C3333(1:3,1:3,1:3,1:3,myPhase)
         lattice_trans_C66(1:6,1:6,myPhase) = math_Mandel3333to66(lattice_trans_C3333(1:3,1:3,1:3,1:3,myPhase))
         do i = 1_pInt, 6_pInt
           if (abs(lattice_trans_C66(i,i,myPhase))<tol_math_check) &
           call IO_error(135_pInt,el=i,ip=myPhase,ext_msg='matrix diagonal "el"ement of phase "ip"')
         enddo
       case (LATTICE_hex_ID)
         c11bar = (lattice_C66(1,1,myPhase) + lattice_C66(1,2,myPhase) + 2.0_pReal*lattice_C66(4,4,myPhase))/2.0_pReal
         c12bar = (lattice_C66(1,1,myPhase) + 5.0_pReal*lattice_C66(1,2,myPhase) - 2.0_pReal*lattice_C66(4,4,myPhase))/6.0_pReal
         c33bar = (lattice_C66(1,1,myPhase) + 2.0_pReal*lattice_C66(1,2,myPhase) + 4.0_pReal*lattice_C66(4,4,myPhase))/3.0_pReal
         c13bar = (lattice_C66(1,1,myPhase) + 2.0_pReal*lattice_C66(1,2,myPhase) - 2.0_pReal*lattice_C66(4,4,myPhase))/3.0_pReal
         c44bar = (lattice_C66(1,1,myPhase) - lattice_C66(1,2,myPhase) + lattice_C66(4,4,myPhase))/3.0_pReal
         c14bar = (lattice_C66(1,1,myPhase) - lattice_C66(1,2,myPhase) - 2.0_pReal*lattice_C66(4,4,myPhase)) &
                                                             /(3.0_pReal*sqrt(2.0_pReal))
         A = c14bar**(2.0_pReal)/c44bar
         B = c14bar**(2.0_pReal)/(0.5_pReal*(c11bar - c12bar))
         lattice_trans_C66(1,1,myPhase) = c11bar - A
         lattice_trans_C66(1,2,myPhase) = c12bar + A
         lattice_trans_C66(1,3,myPhase) = c13bar
         lattice_trans_C66(3,3,myPhase) = c33bar
         lattice_trans_C66(4,4,myPhase) = c44bar - B
         
         lattice_trans_C66(1:6,1:6,myPhase) = lattice_symmetrizeC66(trans_lattice_structure(myPhase),&
                                                                    lattice_trans_C66(1:6,1:6,myPhase))
         lattice_trans_mu(myPhase) = 0.2_pReal *(  lattice_trans_C66(1,1,myPhase) &
                                                 - lattice_trans_C66(1,2,myPhase) &
                                                 + 3.0_pReal*lattice_trans_C66(4,4,myPhase))
         lattice_trans_nu(myPhase) = (  lattice_trans_C66(1,1,myPhase) &
                                      + 4.0_pReal*lattice_trans_C66(1,2,myPhase) &
                                      - 2.0_pReal*lattice_trans_C66(4,4,myPhase)) &
                                                             /(  4.0_pReal*lattice_trans_C66(1,1,myPhase) &
                                                               + 6.0_pReal*lattice_trans_C66(1,2,myPhase) &
                                                               + 2.0_pReal*lattice_trans_C66(4,4,myPhase))
         lattice_trans_C3333(1:3,1:3,1:3,1:3,myPhase) = math_Voigt66to3333(lattice_trans_C66(1:6,1:6,myPhase))
         lattice_trans_C66(1:6,1:6,myPhase) = math_Mandel3333to66(lattice_trans_C3333(1:3,1:3,1:3,1:3,myPhase))
         do i = 1_pInt, 6_pInt
           if (abs(lattice_trans_C66(i,i,myPhase))<tol_math_check) &
           call IO_error(135_pInt,el=i,ip=myPhase,ext_msg='matrix diagonal "el"ement of phase "ip"')
         enddo
     end select
 end select

 lattice_thermalConductivity33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                       lattice_thermalConductivity33(1:3,1:3,myPhase))
 lattice_thermalExpansion33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                    lattice_thermalExpansion33(1:3,1:3,myPhase))
 lattice_DamageDiffusion33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                 lattice_DamageDiffusion33(1:3,1:3,myPhase))
 lattice_vacancyfluxDiffusion33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                 lattice_vacancyfluxDiffusion33(1:3,1:3,myPhase))
 lattice_vacancyfluxMobility33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                 lattice_vacancyfluxMobility33(1:3,1:3,myPhase))
 lattice_PorosityDiffusion33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                 lattice_PorosityDiffusion33(1:3,1:3,myPhase))
 lattice_hydrogenfluxDiffusion33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                 lattice_hydrogenfluxDiffusion33(1:3,1:3,myPhase))
 lattice_hydrogenfluxMobility33(1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                 lattice_hydrogenfluxMobility33(1:3,1:3,myPhase))
 
 select case(lattice_structure(myPhase))
!--------------------------------------------------------------------------------------------------
! fcc
   case (LATTICE_fcc_ID)
     myNslip  = lattice_fcc_Nslip
     myNtwin  = lattice_fcc_Ntwin
     myNtrans = lattice_fcc_Ntrans
     myNcleavage = lattice_fcc_Ncleavage
     do i = 1_pInt,myNslip                                                                           ! assign slip system vectors
       sd(1:3,i) = lattice_fcc_systemSlip(1:3,i)
       sn(1:3,i) = lattice_fcc_systemSlip(4:6,i)
     enddo  
     do i = 1_pInt,myNtwin                                                                           ! assign twin system vectors and shears
       td(1:3,i) = lattice_fcc_systemTwin(1:3,i)
       tn(1:3,i) = lattice_fcc_systemTwin(4:6,i)
       ts(i)     = lattice_fcc_shearTwin(i)
     enddo
     do i = 1_pInt, myNcleavage                                                                      ! assign cleavage system vectors
       cd(1:3,i) = lattice_fcc_systemCleavage(1:3,i)/norm2(lattice_fcc_systemCleavage(1:3,i))
       cn(1:3,i) = lattice_fcc_systemCleavage(4:6,i)/norm2(lattice_fcc_systemCleavage(4:6,i))
       ct(1:3,i) = math_crossproduct(cd(1:3,i),cn(1:3,i))
     enddo  

     ! Phase transformation
     select case(trans_lattice_structure(myPhase))
       case (LATTICE_bcc_ID)                                                                         ! fcc to bcc transformation                         
         do i = 1_pInt,myNtrans
           Rtr(1:3,1:3,i) = math_axisAngleToR(lattice_fccTobcc_systemTrans(1:3,i), &                 ! Pitsch rotation
                              lattice_fccTobcc_systemTrans(4,i)*INRAD)
           Btr(1:3,1:3,i) = math_axisAngleToR(lattice_fccTobcc_bainRot(1:3,i), &                     ! Rotation of fcc to Bain coordinate system
                              lattice_fccTobcc_bainRot(4,i)*INRAD)
           xtr(1:3,i) = real(LATTICE_fccTobcc_bainVariant(1:3,i),pReal)
           ytr(1:3,i) = real(LATTICE_fccTobcc_bainVariant(4:6,i),pReal)
           ztr(1:3,i) = real(LATTICE_fccTobcc_bainVariant(7:9,i),pReal)
           Utr(1:3,1:3,i) = 0.0_pReal                                                                ! Bain deformation
           if ((a_fcc > 0.0_pReal) .and. (a_bcc > 0.0_pReal)) then
             Utr(1:3,1:3,i) = (a_bcc/a_fcc)*math_tensorproduct33(xtr(1:3,i), xtr(1:3,i)) + &
              sqrt(2.0_pReal)*(a_bcc/a_fcc)*math_tensorproduct33(ytr(1:3,i), ytr(1:3,i)) + &
              sqrt(2.0_pReal)*(a_bcc/a_fcc)*math_tensorproduct33(ztr(1:3,i), ztr(1:3,i))
           endif
           Qtr(1:3,1:3,i) = math_mul33x33(Rtr(1:3,1:3,i), Btr(1:3,1:3,i))
           Str(1:3,1:3,i) = math_mul33x33(Rtr(1:3,1:3,i), Utr(1:3,1:3,i)) - MATH_I3
         enddo
       case (LATTICE_hex_ID)
         sstr(1:3,1:3) = MATH_I3
         sstr(1,3)     = sqrt(2.0_pReal)/4.0_pReal
         sdtr(1:3,1:3) = MATH_I3
         if (CoverA_trans > 1.0_pReal .and. CoverA_trans < 2.0_pReal) then
           sdtr(3,3) = CoverA_trans/sqrt(8.0_pReal/3.0_pReal)
         endif
         sttr = math_mul33x33(sdtr, sstr)
         do i = 1_pInt,myNtrans
           xtr(1:3,i) = lattice_fccTohex_systemTrans(1:3,i)/norm2(lattice_fccTohex_systemTrans(1:3,i))
           ztr(1:3,i) = lattice_fccTohex_systemTrans(4:6,i)/norm2(lattice_fccTohex_systemTrans(4:6,i))
           ytr(1:3,i) = -math_crossproduct(xtr(1:3,i), ztr(1:3,i))
           Rtr(1:3,1,i) = xtr(1:3,i)
           Rtr(1:3,2,i) = ytr(1:3,i)
           Rtr(1:3,3,i) = ztr(1:3,i)
           Qtr(1:3,1:3,i) = Rtr(1:3,1:3,i)
           Str(1:3,1:3,i) = math_mul33x33(Rtr(1:3,1:3,i), math_mul33x33(sttr, math_transpose33(Rtr(1:3,1:3,i))))
           Str(1:3,1:3,i) = Str(1:3,1:3,i) - MATH_I3
           trs(i)     = lattice_fccTohex_shearTrans(i)
         enddo
       case default
         Qtr = 0.0_pReal
         Str = 0.0_pReal
     end select

     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)          = lattice_fcc_NslipSystem
     lattice_NtwinSystem(1:lattice_maxNtwinFamily,myPhase)          = lattice_fcc_NtwinSystem
     lattice_NtransSystem(1:lattice_maxNtransFamily,myPhase)        = lattice_fcc_NtransSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_fcc_NcleavageSystem
     lattice_NnonSchmid(myPhase)                                    = lattice_fcc_NnonSchmid
     lattice_interactionSlipSlip(1:myNslip,1:myNslip,myPhase)       = lattice_fcc_interactionSlipSlip
     lattice_interactionSlipTwin(1:myNslip,1:myNtwin,myPhase)       = lattice_fcc_interactionSlipTwin
     lattice_interactionTwinSlip(1:myNtwin,1:myNslip,myPhase)       = lattice_fcc_interactionTwinSlip
     lattice_interactionTwinTwin(1:myNtwin,1:myNtwin,myPhase)       = lattice_fcc_interactionTwinTwin
     lattice_interactionSlipTrans(1:myNslip,1:myNtrans,myPhase)     = lattice_fccTohex_interactionSlipTrans
     lattice_interactionTransSlip(1:myNtrans,1:myNslip,myPhase)     = lattice_fccTohex_interactionTransSlip
     lattice_interactionTransTrans(1:myNtrans,1:myNtrans,myPhase)   = lattice_fccTohex_interactionTransTrans
     lattice_projectionTrans(1:myNtrans,1:myNtrans,myPhase)         = LATTICE_fccTobcc_projectionTrans*&
       LATTICE_fccTobcc_projectionTransFactor

!--------------------------------------------------------------------------------------------------
! bcc
   case (LATTICE_bcc_ID)
     myNslip = lattice_bcc_Nslip
     myNtwin = lattice_bcc_Ntwin
     myNtrans = lattice_bcc_Ntrans
     myNcleavage = lattice_bcc_Ncleavage
     do i = 1_pInt,myNslip                                                                          ! assign slip system vectors
       sd(1:3,i) = lattice_bcc_systemSlip(1:3,i)
       sn(1:3,i) = lattice_bcc_systemSlip(4:6,i)
       sdU = sd(1:3,i) / norm2(sd(1:3,i))
       snU = sn(1:3,i) / norm2(sn(1:3,i))
       ! "np" and "nn" according to Gröger_etal2008, Acta Materialia 56 (2008) 5412–5425, table 1 (corresponds to their "n1" for positive and negative slip direction respectively)
       np = math_mul33x3(math_axisAngleToR(sdU,60.0_pReal*INRAD), snU)
       nn = math_mul33x3(math_axisAngleToR(-sdU,60.0_pReal*INRAD), snU)
         ! Schmid matrices with non-Schmid contributions according to Koester_etal2012, Acta Materialia 60 (2012) 3894–3901, eq. (17) ("n1" is replaced by either "np" or "nn" according to either positive or negative slip direction)
       sns(1:3,1:3,1,1,i) = math_tensorproduct33(sdU, np)
       sns(1:3,1:3,2,1,i) = math_tensorproduct33(-sdU, nn)
       sns(1:3,1:3,1,2,i) = math_tensorproduct33(math_crossproduct(snU, sdU), snU)
       sns(1:3,1:3,2,2,i) = math_tensorproduct33(math_crossproduct(snU, -sdU), snU)
       sns(1:3,1:3,1,3,i) = math_tensorproduct33(math_crossproduct(np, sdU), np)
       sns(1:3,1:3,2,3,i) = math_tensorproduct33(math_crossproduct(nn, -sdU), nn)
       sns(1:3,1:3,1,4,i) = math_tensorproduct33(snU, snU)
       sns(1:3,1:3,2,4,i) = math_tensorproduct33(snU, snU)
       sns(1:3,1:3,1,5,i) = math_tensorproduct33(math_crossproduct(snU, sdU), math_crossproduct(snU, sdU))
       sns(1:3,1:3,2,5,i) = math_tensorproduct33(math_crossproduct(snU, -sdU), math_crossproduct(snU, -sdU))
       sns(1:3,1:3,1,6,i) = math_tensorproduct33(sdU, sdU)
       sns(1:3,1:3,2,6,i) = math_tensorproduct33(-sdU, -sdU)
     enddo
     do i = 1_pInt,myNtwin                                                                          ! assign twin system vectors and shears
       td(1:3,i) = lattice_bcc_systemTwin(1:3,i)
       tn(1:3,i) = lattice_bcc_systemTwin(4:6,i)
       ts(i)     = lattice_bcc_shearTwin(i)
     enddo
     do i = 1_pInt, myNcleavage                                                                      ! assign cleavage system vectors
       cd(1:3,i) = lattice_bcc_systemCleavage(1:3,i)/norm2(lattice_bcc_systemCleavage(1:3,i))
       cn(1:3,i) = lattice_bcc_systemCleavage(4:6,i)/norm2(lattice_bcc_systemCleavage(4:6,i))
       ct(1:3,i) = math_crossproduct(cd(1:3,i),cn(1:3,i))
     enddo  
     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)          = lattice_bcc_NslipSystem
     lattice_NtwinSystem(1:lattice_maxNtwinFamily,myPhase)          = lattice_bcc_NtwinSystem
     lattice_NtransSystem(1:lattice_maxNtransFamily,myPhase)        = lattice_bcc_NtransSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_bcc_NcleavageSystem
     lattice_NnonSchmid(myPhase)                                    = lattice_bcc_NnonSchmid
     lattice_interactionSlipSlip(1:myNslip,1:myNslip,myPhase)       = lattice_bcc_interactionSlipSlip
     lattice_interactionSlipTwin(1:myNslip,1:myNtwin,myPhase)       = lattice_bcc_interactionSlipTwin
     lattice_interactionTwinSlip(1:myNtwin,1:myNslip,myPhase)       = lattice_bcc_interactionTwinSlip
     lattice_interactionTwinTwin(1:myNtwin,1:myNtwin,myPhase)       = lattice_bcc_interactionTwinTwin
     
!--------------------------------------------------------------------------------------------------
! hex (including conversion from miller-bravais (a1=a2=a3=c) to miller (a, b, c) indices)
   case (LATTICE_hex_ID)
     myNslip = lattice_hex_Nslip
     myNtwin = lattice_hex_Ntwin
     myNtrans = lattice_hex_Ntrans
     myNcleavage = lattice_hex_Ncleavage
     do i = 1_pInt,myNslip                                                                          ! assign slip system vectors                                                         
       sd(1,i) =  lattice_hex_systemSlip(1,i)*1.5_pReal                                             ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]
       sd(2,i) = (lattice_hex_systemSlip(1,i)+2.0_pReal*lattice_hex_systemSlip(2,i))*&
                                                                      0.5_pReal*sqrt(3.0_pReal)
       sd(3,i) =  lattice_hex_systemSlip(4,i)*CoverA
       sn(1,i) =  lattice_hex_systemSlip(5,i)                                                       ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))
       sn(2,i) = (lattice_hex_systemSlip(5,i)+2.0_pReal*lattice_hex_systemSlip(6,i))/sqrt(3.0_pReal)
       sn(3,i) =  lattice_hex_systemSlip(8,i)/CoverA
     enddo  
     do i = 1_pInt,myNtwin                                                                          ! assign twin system vectors and shears
       td(1,i) =  lattice_hex_systemTwin(1,i)*1.5_pReal
       td(2,i) = (lattice_hex_systemTwin(1,i)+2.0_pReal*lattice_hex_systemTwin(2,i))*&
                                                                     0.5_pReal*sqrt(3.0_pReal)
       td(3,i) =  lattice_hex_systemTwin(4,i)*CoverA
       tn(1,i) =  lattice_hex_systemTwin(5,i)
       tn(2,i) = (lattice_hex_systemTwin(5,i)+2.0_pReal*lattice_hex_systemTwin(6,i))/sqrt(3.0_pReal)
       tn(3,i) =  lattice_hex_systemTwin(8,i)/CoverA
       select case(lattice_hex_shearTwin(i))                                                        ! from Christian & Mahajan 1995 p.29
         case (1_pInt)                                                                              ! <-10.1>{10.2}
           ts(i) = (3.0_pReal-CoverA*CoverA)/sqrt(3.0_pReal)/CoverA
         case (2_pInt)                                                                              ! <11.6>{-1-1.1}
           ts(i) = 1.0_pReal/CoverA
         case (3_pInt)                                                                              ! <10.-2>{10.1}
           ts(i) = (4.0_pReal*CoverA*CoverA-9.0_pReal)/4.0_pReal/sqrt(3.0_pReal)/CoverA
         case (4_pInt)                                                                              ! <11.-3>{11.2}
           ts(i) = 2.0_pReal*(CoverA*CoverA-2.0_pReal)/3.0_pReal/CoverA
       end select
     enddo
     do i = 1_pInt, myNcleavage                                                                     ! cleavage system vectors                                                         
       cd(1,i) =  lattice_hex_systemCleavage(1,i)*1.5_pReal                                         ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]
       cd(2,i) = (lattice_hex_systemCleavage(1,i)+2.0_pReal*lattice_hex_systemCleavage(2,i))*&
                                                                      0.5_pReal*sqrt(3.0_pReal)
       cd(3,i) =  lattice_hex_systemCleavage(4,i)*CoverA
       cd(1:3,1) = cd(1:3,i)/norm2(cd(1:3,i))
       cn(1,i) =  lattice_hex_systemCleavage(5,i)                                                   ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))
       cn(2,i) = (lattice_hex_systemCleavage(5,i)+2.0_pReal*lattice_hex_systemCleavage(6,i))/sqrt(3.0_pReal)
       cn(3,i) =  lattice_hex_systemCleavage(8,i)/CoverA
       cn(1:3,1) = cn(1:3,i)/norm2(cn(1:3,i))
       ct(1:3,i) = math_crossproduct(cd(1:3,i),cn(1:3,i))
     enddo  
     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)          = lattice_hex_NslipSystem
     lattice_NtwinSystem(1:lattice_maxNtwinFamily,myPhase)          = lattice_hex_NtwinSystem
     lattice_NtransSystem(1:lattice_maxNtransFamily,myPhase)        = lattice_hex_NtransSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_hex_NcleavageSystem
     lattice_NnonSchmid(myPhase)                                    = lattice_hex_NnonSchmid
     lattice_interactionSlipSlip(1:myNslip,1:myNslip,myPhase)       = lattice_hex_interactionSlipSlip
     lattice_interactionSlipTwin(1:myNslip,1:myNtwin,myPhase)       = lattice_hex_interactionSlipTwin
     lattice_interactionTwinSlip(1:myNtwin,1:myNslip,myPhase)       = lattice_hex_interactionTwinSlip
     lattice_interactionTwinTwin(1:myNtwin,1:myNtwin,myPhase)       = lattice_hex_interactionTwinTwin

!--------------------------------------------------------------------------------------------------
! bct
   case (LATTICE_bct_ID)
     myNslip = lattice_bct_Nslip
     myNtwin = lattice_bct_Ntwin
     myNcleavage = lattice_bct_Ncleavage
     do i = 1_pInt,myNslip                                                                          ! assign slip system vectors
       sd(1:2,i) = lattice_bct_systemSlip(1:2,i)
       sd(3,i) = lattice_bct_systemSlip(3,i)*CoverA
       sn(1:2,i) = lattice_bct_systemSlip(4:5,i)
       sn(3,i) =  lattice_bct_systemSlip(6,i)/CoverA
       sdU = sd(1:3,i) / norm2(sd(1:3,i))
       snU = sn(1:3,i) / norm2(sn(1:3,i))
     enddo
     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)          = lattice_bct_NslipSystem
     lattice_NtwinSystem(1:lattice_maxNtwinFamily,myPhase)          = lattice_bct_NtwinSystem
     lattice_NtransSystem(1:lattice_maxNtransFamily,myPhase)        = lattice_bct_NtransSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_bct_NcleavageSystem
     lattice_NnonSchmid(myPhase)                                    = lattice_bct_NnonSchmid
     lattice_interactionSlipSlip(1:myNslip,1:myNslip,myPhase)       = lattice_bct_interactionSlipSlip

!--------------------------------------------------------------------------------------------------
! orthorhombic (no crystal plasticity)
   case (LATTICE_ort_ID)
     myNslip       = 0_pInt
     myNtwin       = 0_pInt
     myNtrans      = 0_pInt
     myNcleavage   = lattice_ortho_Ncleavage
     do i = 1_pInt, myNcleavage                                                                      ! assign cleavage system vectors
       cd(1:3,i) = lattice_iso_systemCleavage(1:3,i)/norm2(LATTICE_ortho_systemCleavage(1:3,i))
       cn(1:3,i) = lattice_iso_systemCleavage(4:6,i)/norm2(LATTICE_ortho_systemCleavage(4:6,i))
       ct(1:3,i) = math_crossproduct(cd(1:3,i),cn(1:3,i))
     enddo  
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_iso_NcleavageSystem

!--------------------------------------------------------------------------------------------------
! isotropic (no crystal plasticity)
   case (LATTICE_iso_ID)
     myNslip       = 0_pInt
     myNtwin       = 0_pInt
     myNtrans      = 0_pInt
     myNcleavage   = lattice_iso_Ncleavage
     do i = 1_pInt, myNcleavage                                                                      ! assign cleavage system vectors
       cd(1:3,i) = lattice_iso_systemCleavage(1:3,i)/norm2(lattice_iso_systemCleavage(1:3,i))
       cn(1:3,i) = lattice_iso_systemCleavage(4:6,i)/norm2(lattice_iso_systemCleavage(4:6,i))
       ct(1:3,i) = math_crossproduct(cd(1:3,i),cn(1:3,i))
     enddo  
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_iso_NcleavageSystem

!--------------------------------------------------------------------------------------------------
! something went wrong
   case default
     call IO_error(130_pInt,ext_msg='lattice_initializeStructure')
 end select


 do i = 1_pInt,myNslip                                                                              ! store slip system vectors and Schmid matrix for my structure
   lattice_sd(1:3,i,myPhase) = sd(1:3,i)/norm2(sd(1:3,i))                                           ! make unit vector
   lattice_sn(1:3,i,myPhase) = sn(1:3,i)/norm2(sn(1:3,i))                                           ! make unit vector
   lattice_st(1:3,i,myPhase) = math_crossproduct(lattice_sd(1:3,i,myPhase), &
                                                  lattice_sn(1:3,i,myPhase))
   lattice_Sslip(1:3,1:3,1,i,myPhase) = math_tensorproduct33(lattice_sd(1:3,i,myPhase), &
                                                           lattice_sn(1:3,i,myPhase))               ! calculate Schmid matrix d \otimes n
   do j = 1_pInt,lattice_NnonSchmid(myPhase)
     lattice_Sslip(1:3,1:3,2*j  ,i,myPhase) = sns(1:3,1:3,1,j,i)
     lattice_Sslip(1:3,1:3,2*j+1,i,myPhase) = sns(1:3,1:3,2,j,i)
   enddo 
   do j = 1_pInt,1_pInt+2_pInt*lattice_NnonSchmid(myPhase)
     lattice_Sslip_v(1:6,j,i,myPhase) = &
       math_Mandel33to6(math_symmetric33(lattice_Sslip(1:3,1:3,j,i,myPhase)))
   enddo
   if (abs(math_trace33(lattice_Sslip(1:3,1:3,1,i,myPhase))) > tol_math_check) &
     call IO_error(0_pInt,myPhase,i,0_pInt,ext_msg = 'dilatational slip Schmid matrix')
 enddo
 do i = 1_pInt,myNtwin                                                                              ! store twin system vectors and Schmid plus rotation matrix for my structure
   lattice_td(1:3,i,myPhase) = td(1:3,i)/norm2(td(1:3,i))                                           ! make unit vector
   lattice_tn(1:3,i,myPhase) = tn(1:3,i)/norm2(tn(1:3,i))                                           ! make unit vector
   lattice_tt(1:3,i,myPhase) = math_crossproduct(lattice_td(1:3,i,myPhase), &
                                                  lattice_tn(1:3,i,myPhase))
   lattice_Stwin(1:3,1:3,i,myPhase) = math_tensorproduct33(lattice_td(1:3,i,myPhase), &
                                                         lattice_tn(1:3,i,myPhase))
   lattice_Stwin_v(1:6,i,myPhase)   = math_Mandel33to6(math_symmetric33(lattice_Stwin(1:3,1:3,i,myPhase)))
   lattice_Qtwin(1:3,1:3,i,myPhase) = math_axisAngleToR(tn(1:3,i),180.0_pReal*INRAD)
   lattice_shearTwin(i,myPhase)     = ts(i)
   if (abs(math_trace33(lattice_Stwin(1:3,1:3,i,myPhase))) > tol_math_check) &
     call IO_error(301_pInt,myPhase,ext_msg = 'dilatational twin Schmid matrix')
 enddo
 do i = 1_pInt,myNtrans
   lattice_Qtrans(1:3,1:3,i,myPhase) = Qtr(1:3,1:3,i)
   lattice_Strans(1:3,1:3,i,myPhase) = Str(1:3,1:3,i)
   lattice_Strans_v(1:6,i,myPhase)   = math_Mandel33to6(math_symmetric33(lattice_Strans(1:3,1:3,i,myPhase)))
   lattice_shearTrans(i,myPhase)     = trs(i)
 enddo
 do i = 1_pInt,myNcleavage                                                                         ! store slip system vectors and Schmid matrix for my structure
   lattice_Scleavage(1:3,1:3,1,i,myPhase) = math_tensorproduct33(cd(1:3,i),cn(1:3,i))
   lattice_Scleavage(1:3,1:3,2,i,myPhase) = math_tensorproduct33(ct(1:3,i),cn(1:3,i))
   lattice_Scleavage(1:3,1:3,3,i,myPhase) = math_tensorproduct33(cn(1:3,i),cn(1:3,i))
   do j = 1_pInt,3_pInt
     lattice_Scleavage_v(1:6,j,i,myPhase) = &
       math_Mandel33to6(math_symmetric33(lattice_Scleavage(1:3,1:3,j,i,myPhase)))
   enddo
 enddo
      
end subroutine lattice_initializeStructure


!--------------------------------------------------------------------------------------------------
!> @brief Symmetrizes stiffness matrix according to lattice type
!--------------------------------------------------------------------------------------------------
pure function lattice_symmetrizeC66(struct,C66)

 implicit none
 integer(kind(LATTICE_undefined_ID)), intent(in) :: struct
 real(pReal), dimension(6,6), intent(in) :: C66
 real(pReal), dimension(6,6) :: lattice_symmetrizeC66
 integer(pInt) :: j,k

 lattice_symmetrizeC66 = 0.0_pReal
 
 select case(struct)
   case (LATTICE_iso_ID)
     forall(k=1_pInt:3_pInt)
       forall(j=1_pInt:3_pInt) lattice_symmetrizeC66(k,j) = C66(1,2)
       lattice_symmetrizeC66(k,k) = C66(1,1)
       lattice_symmetrizeC66(k+3,k+3) = 0.5_pReal*(C66(1,1)-C66(1,2))
     end forall
   case (LATTICE_fcc_ID,LATTICE_bcc_ID)
     forall(k=1_pInt:3_pInt)
       forall(j=1_pInt:3_pInt) lattice_symmetrizeC66(k,j) =   C66(1,2)
       lattice_symmetrizeC66(k,k) =     C66(1,1)
       lattice_symmetrizeC66(k+3_pInt,k+3_pInt) = C66(4,4)
     end forall    
   case (LATTICE_hex_ID)
     lattice_symmetrizeC66(1,1) = C66(1,1)
     lattice_symmetrizeC66(2,2) = C66(1,1)
     lattice_symmetrizeC66(3,3) = C66(3,3)
     lattice_symmetrizeC66(1,2) = C66(1,2)
     lattice_symmetrizeC66(2,1) = C66(1,2)
     lattice_symmetrizeC66(1,3) = C66(1,3)
     lattice_symmetrizeC66(3,1) = C66(1,3)
     lattice_symmetrizeC66(2,3) = C66(1,3)
     lattice_symmetrizeC66(3,2) = C66(1,3)
     lattice_symmetrizeC66(4,4) = C66(4,4)
     lattice_symmetrizeC66(5,5) = C66(4,4)
     lattice_symmetrizeC66(6,6) = 0.5_pReal*(C66(1,1)-C66(1,2))
   case (LATTICE_ort_ID)
     lattice_symmetrizeC66(1,1) = C66(1,1)
     lattice_symmetrizeC66(2,2) = C66(2,2)
     lattice_symmetrizeC66(3,3) = C66(3,3)
     lattice_symmetrizeC66(1,2) = C66(1,2)
     lattice_symmetrizeC66(2,1) = C66(1,2)
     lattice_symmetrizeC66(1,3) = C66(1,3)
     lattice_symmetrizeC66(3,1) = C66(1,3)
     lattice_symmetrizeC66(2,3) = C66(2,3)
     lattice_symmetrizeC66(3,2) = C66(2,3)
     lattice_symmetrizeC66(4,4) = C66(4,4)
     lattice_symmetrizeC66(5,5) = C66(5,5)
     lattice_symmetrizeC66(6,6) = C66(6,6)
   case (LATTICE_bct_ID)
     lattice_symmetrizeC66(1,1) = C66(1,1)
     lattice_symmetrizeC66(2,2) = C66(1,1)
     lattice_symmetrizeC66(3,3) = C66(3,3)
     lattice_symmetrizeC66(1,2) = C66(1,2)
     lattice_symmetrizeC66(2,1) = C66(1,2)
     lattice_symmetrizeC66(1,3) = C66(1,3)
     lattice_symmetrizeC66(3,1) = C66(1,3)
     lattice_symmetrizeC66(2,3) = C66(1,3)
     lattice_symmetrizeC66(3,2) = C66(1,3)
     lattice_symmetrizeC66(4,4) = C66(4,4)
     lattice_symmetrizeC66(5,5) = C66(4,4)
     lattice_symmetrizeC66(6,6) = C66(6,6)        !J. A. Rayne and B. S. Chandrasekhar Phys. Rev. 120, 1658 Erratum Phys. Rev. 122, 1962
   case default
     lattice_symmetrizeC66 = C66
  end select
  
 end function lattice_symmetrizeC66

!--------------------------------------------------------------------------------------------------
!> @brief Symmetrizes 2nd order tensor according to lattice type
!--------------------------------------------------------------------------------------------------
pure function lattice_symmetrize33(struct,T33)

 implicit none
 integer(kind(LATTICE_undefined_ID)), intent(in) :: struct
 real(pReal), dimension(3,3), intent(in) :: T33
 real(pReal), dimension(3,3) :: lattice_symmetrize33
 integer(pInt) :: k

 lattice_symmetrize33 = 0.0_pReal
 
 select case(struct)
   case (LATTICE_iso_ID,LATTICE_fcc_ID,LATTICE_bcc_ID)
     forall(k=1_pInt:3_pInt) lattice_symmetrize33(k,k) = T33(1,1)
   case (LATTICE_hex_ID)
     lattice_symmetrize33(1,1) = T33(1,1)
     lattice_symmetrize33(2,2) = T33(1,1)
     lattice_symmetrize33(3,3) = T33(3,3)
   case (LATTICE_ort_ID,lattice_bct_ID)
     lattice_symmetrize33(1,1) = T33(1,1)
     lattice_symmetrize33(2,2) = T33(2,2)
     lattice_symmetrize33(3,3) = T33(3,3)
   case default
     lattice_symmetrize33 = T33
  end select
  
 end function lattice_symmetrize33


!--------------------------------------------------------------------------------------------------
!> @brief figures whether unit quat falls into stereographic standard triangle
!--------------------------------------------------------------------------------------------------
logical pure function lattice_qInSST(Q, struct)
 use prec, only: &
   prec_isNaN
 use math, only: &
   math_qToRodrig

 implicit none
 real(pReal), dimension(4), intent(in) ::      Q                                                    ! orientation
 integer(kind(LATTICE_undefined_ID)), intent(in) :: struct                                          ! lattice structure
 real(pReal), dimension(3) ::                  Rodrig                                               ! Rodrigues vector of Q

 Rodrig = math_qToRodrig(Q)
 if (any(prec_isNaN(Rodrig))) then
   lattice_qInSST = .false.
 else
   select case (struct)
     case (LATTICE_bcc_ID,LATTICE_fcc_ID)
       lattice_qInSST = Rodrig(1) > Rodrig(2) .and. &
                        Rodrig(2) > Rodrig(3) .and. &
                        Rodrig(3) > 0.0_pReal
     case (LATTICE_hex_ID)
       lattice_qInSST = Rodrig(1) > sqrt(3.0_pReal)*Rodrig(2) .and. &
                        Rodrig(2) > 0.0_pReal .and. &
                        Rodrig(3) > 0.0_pReal
     case default
       lattice_qInSST = .true.
   end select
 endif

end function lattice_qInSST


!--------------------------------------------------------------------------------------------------
!> @brief calculates the disorientation for 2 unit quaternions
!--------------------------------------------------------------------------------------------------
pure function lattice_qDisorientation(Q1, Q2, struct)
 use prec, only: &
  tol_math_check
 use math, only: &
   math_qMul, &
   math_qConj

 implicit none
 real(pReal), dimension(4) ::                  lattice_qDisorientation
 real(pReal), dimension(4), intent(in) :: &
   Q1, &                                                                                            ! 1st orientation
                                               Q2                                                   ! 2nd orientation
 integer(kind(LATTICE_undefined_ID)), optional, intent(in) :: &                                     ! if given, symmetries between the two orientation will be considered
   struct

 real(pReal), dimension(4) ::                  dQ,dQsymA,mis
 integer(pInt)    ::                           i,j,k,s,symmetry
 integer(kind(LATTICE_undefined_ID)) :: myStruct

!--------------------------------------------------------------------------------------------------
! check if a structure with known symmetries is given
 if (present(struct)) then
   myStruct = struct
   select case (struct)
     case(LATTICE_fcc_ID,LATTICE_bcc_ID)
       symmetry = 1_pInt
     case(LATTICE_hex_ID)
       symmetry = 2_pInt
     case default
       symmetry = 0_pInt
   end select
 else
   symmetry = 0_pInt
   myStruct = LATTICE_undefined_ID
 endif


!--------------------------------------------------------------------------------------------------
! calculate misorientation, for cubic and hexagonal structure find symmetries
 dQ = math_qMul(math_qConj(Q1),Q2)
 lattice_qDisorientation = dQ

 select case(symmetry) 

    case (1_pInt,2_pInt)
     s = sum(lattice_NsymOperations(1:symmetry-1_pInt))
      do i = 1_pInt,2_pInt
        dQ = math_qConj(dQ)                                                                         ! switch order of "from -- to"
       do j = 1_pInt,lattice_NsymOperations(symmetry)                                               ! run through first crystal's symmetries
          dQsymA = math_qMul(lattice_symOperations(1:4,s+j),dQ)                                     ! apply sym
         do k = 1_pInt,lattice_NsymOperations(symmetry)                                             ! run through 2nd crystal's symmetries
            mis = math_qMul(dQsymA,lattice_symOperations(1:4,s+k))                                  ! apply sym
            if (mis(1) < 0.0_pReal) &                                                               ! want positive angle
              mis = -mis
            if (mis(1)-lattice_qDisorientation(1) > -tol_math_check &
             .and. lattice_qInSST(mis,LATTICE_undefined_ID)) lattice_qDisorientation = mis          ! found better one
      enddo; enddo; enddo
   case (0_pInt)
     if (lattice_qDisorientation(1) < 0.0_pReal) lattice_qDisorientation = -lattice_qDisorientation ! keep omega within 0 to 180 deg
  end select

end function lattice_qDisorientation

end module lattice
