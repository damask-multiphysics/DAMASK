!* $Id$
!************************************
!*         Module: LATTICE          *
!************************************
!* contains:                        *
!* - Lattice structure definition   *
!* - Slip system definition         *
!* - Schmid matrices calculation    *
!************************************

MODULE lattice

!*** Include other modules ***
use prec, only: pReal,pInt
implicit none

!************************************
!*      Lattice structures          *
!************************************

integer(pInt) lattice_Nhexagonal, &                               ! # of hexagonal lattice structure (from tag CoverA_ratio)
              lattice_Nstructure                                  ! # of lattice structures (1: fcc,2: bcc,3+: hexagonal)
integer(pInt), parameter :: lattice_maxNslipFamily = 4            ! max # of slip system families over lattice structures
integer(pInt), parameter :: lattice_maxNtwinFamily = 4            ! max # of twin system families over lattice structures
integer(pInt), parameter :: lattice_maxNslip = 48                 ! max # of slip systems over lattice structures
integer(pInt), parameter :: lattice_maxNtwin = 24                 ! max # of twin systems over lattice structures
integer(pInt), parameter :: lattice_maxNinteraction = 20          ! max # of interaction types (in hardening matrix part)

integer(pInt), pointer, dimension(:,:) :: interactionSlipSlip, &
                                          interactionSlipTwin, &
                                          interactionTwinSlip, &
                                          interactionTwinTwin

! Schmid matrices, normal, shear direction and d x n of slip systems
real(pReal), allocatable, dimension(:,:,:,:) :: lattice_Sslip
real(pReal), allocatable, dimension(:,:,:)   :: lattice_Sslip_v
real(pReal), allocatable, dimension(:,:,:)   :: lattice_sn, &
                                                lattice_sd, &
                                                lattice_st

! rotation and Schmid matrices, normal, shear direction and d x n of twin systems
real(pReal), allocatable, dimension(:,:,:,:) :: lattice_Qtwin
real(pReal), allocatable, dimension(:,:,:,:) :: lattice_Stwin
real(pReal), allocatable, dimension(:,:,:)   :: lattice_Stwin_v
real(pReal), allocatable, dimension(:,:,:)   :: lattice_tn, &
                                                lattice_td, &
                                                lattice_tt

! characteristic twin shear
real(pReal), allocatable, dimension(:,:)     :: lattice_shearTwin

! number of slip and twin systems in each family
integer(pInt), allocatable, dimension(:,:)   :: lattice_NslipSystem, &
                                                lattice_NtwinSystem

! interaction type of slip and twin systems among each other
integer(pInt), allocatable, dimension(:,:,:) :: lattice_interactionSlipSlip, &
                                                lattice_interactionSlipTwin, &
                                                lattice_interactionTwinSlip, &
                                                lattice_interactionTwinTwin


!============================== fcc (1) =================================

 integer(pInt), parameter, dimension(lattice_maxNslipFamily) :: lattice_fcc_NslipSystem = (/12, 0, 0, 0/)
 integer(pInt), parameter, dimension(lattice_maxNtwinFamily) :: lattice_fcc_NtwinSystem = (/12, 0, 0, 0/)
 integer(pInt), parameter :: lattice_fcc_Nslip = 12                                       ! sum(lattice_fcc_NslipSystem)
 integer(pInt), parameter :: lattice_fcc_Ntwin = 12                                       ! sum(lattice_fcc_NtwinSystem)
 integer(pInt) ::            lattice_fcc_Nstructure = 0_pInt

 real(pReal), dimension(3+3,lattice_fcc_Nslip), parameter :: lattice_fcc_systemSlip = &
 reshape((/&
! Slip system <110>{111}  Sorted according to Eisenlohr & Hantcherli
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
   /),(/3+3,lattice_fcc_Nslip/))

 real(pReal), dimension(3+3,lattice_fcc_Ntwin), parameter :: lattice_fcc_systemTwin = &
 reshape((/&
! Twin system <112>{111}  Sorted according to Eisenlohr & Hantcherli
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
   /),(/3+3,lattice_fcc_Ntwin/))

 real(pReal), dimension(lattice_fcc_Ntwin), parameter :: lattice_fcc_shearTwin = &
 reshape((/&
! Twin system <112>{111}  Sorted according to Eisenlohr & Hantcherli
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812, &
 0.7071067812  &
   /),(/lattice_fcc_Ntwin/))

 integer(pInt), target, dimension(lattice_fcc_Nslip,lattice_fcc_Nslip) :: lattice_fcc_interactionSlipSlip = &
 reshape((/&
! Interaction types
! 1 --- self interaction
! 2 --- coplanar interaction
! 3 --- collinear interaction
! 4 --- Hirth locks
! 5 --- glissile junctions
! 6 --- Lomer locks
 1,2,2,4,6,5,3,5,5,4,5,6, &
 2,1,2,6,4,5,5,4,6,5,3,5, &
 2,2,1,5,5,3,5,6,4,6,5,4, &
 4,6,5,1,2,2,4,5,6,3,5,5, &
 6,4,5,2,1,2,5,3,5,5,4,6, &
 5,5,3,2,2,1,6,5,4,5,6,4, &
 3,5,5,4,5,6,1,2,2,4,6,5, &
 5,4,6,5,3,5,2,1,2,6,4,5, &
 5,6,4,6,5,4,2,2,1,5,5,3, &
 4,5,6,3,5,5,4,6,5,1,2,2, &
 5,3,5,5,4,6,6,4,5,2,1,2, &
 6,5,4,5,6,4,5,5,3,2,2,1  &
   /),(/lattice_fcc_Nslip,lattice_fcc_Nslip/))

 integer(pInt), target, dimension(lattice_fcc_Ntwin,lattice_fcc_Nslip) :: lattice_fcc_interactionSlipTwin = &
 reshape((/&
 1,1,1,2,2,1,1,2,2,2,1,2, &
 1,1,1,2,2,1,1,2,2,2,1,2, &
 1,1,1,2,2,1,1,2,2,2,1,2, &
 2,2,1,1,1,1,2,1,2,1,2,2, &
 2,2,1,1,1,1,2,1,2,1,2,2, &
 2,2,1,1,1,1,2,1,2,1,2,2, &
 1,2,2,2,1,2,1,1,1,2,2,1, &
 1,2,2,2,1,2,1,1,1,2,2,1, &
 1,2,2,2,1,2,1,1,1,2,2,1, &
 2,1,2,1,2,2,2,2,1,1,1,1, &
 2,1,2,1,2,2,2,2,1,1,1,1, &
 2,1,2,1,2,2,2,2,1,1,1,1  &
   /),(/lattice_fcc_Ntwin,lattice_fcc_Nslip/))

 integer(pInt), target, dimension(lattice_fcc_Nslip,lattice_fcc_Ntwin) :: lattice_fcc_interactionTwinSlip = 0

 integer(pInt), target, dimension(lattice_fcc_Ntwin,lattice_fcc_Ntwin) :: lattice_fcc_interactionTwinTwin = &
 reshape((/&
 1,1,1,2,2,2,2,2,2,2,2,2, &
 1,1,1,2,2,2,2,2,2,2,2,2, &
 1,1,1,2,2,2,2,2,2,2,2,2, &
 2,2,2,1,1,1,2,2,2,2,2,2, &
 2,2,2,1,1,1,2,2,2,2,2,2, &
 2,2,2,1,1,1,2,2,2,2,2,2, &
 2,2,2,2,2,2,1,1,1,2,2,2, &
 2,2,2,2,2,2,1,1,1,2,2,2, &
 2,2,2,2,2,2,1,1,1,2,2,2, &
 2,2,2,2,2,2,2,2,2,1,1,1, &
 2,2,2,2,2,2,2,2,2,1,1,1, &
 2,2,2,2,2,2,2,2,2,1,1,1  &
   /),(/lattice_fcc_Ntwin,lattice_fcc_Ntwin/))
   

!============================== bcc (2) =================================

 integer(pInt), parameter, dimension(lattice_maxNslipFamily) :: lattice_bcc_NslipSystem = (/12,12,24, 0/)
 integer(pInt), parameter, dimension(lattice_maxNtwinFamily) :: lattice_bcc_NtwinSystem = (/12, 0, 0, 0/)
 integer(pInt), parameter :: lattice_bcc_Nslip = 48                                       ! sum(lattice_bcc_NslipSystem)
 integer(pInt), parameter :: lattice_bcc_Ntwin = 12                                       ! sum(lattice_bcc_NtwinSystem)
 integer(pInt) ::            lattice_bcc_Nstructure = 0_pInt

 real(pReal), dimension(3+3,lattice_bcc_Nslip), parameter :: lattice_bcc_systemSlip = &
 reshape((/&
! Slip system <111>{110}  meaningful sorting?
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
! Slip system <111>{112} meaningful sorting ?
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
  1, 1, 1,     1, 1,-2, &
! Slip system <111>{123} meaningful sorting ?
  1, 1,-1,     1, 2, 3, &
  1,-1, 1,    -1, 2, 3, &
 -1, 1, 1,     1,-2, 3, &
  1, 1, 1,     1, 2,-3, &
  1,-1, 1,     1, 3, 2, &
  1, 1,-1,    -1, 3, 2, &
  1, 1, 1,     1,-3, 2, &
 -1, 1, 1,     1, 3,-2, &
  1, 1,-1,     2, 1, 3, &
  1,-1, 1,    -2, 1, 3, &
 -1, 1, 1,     2,-1, 3, &
  1, 1, 1,     2, 1,-3, &
  1,-1, 1,     2, 3, 1, &
  1, 1,-1,    -2, 3, 1, &
  1, 1, 1,     2,-3, 1, &
 -1, 1, 1,     2, 3,-1, &
 -1, 1, 1,     3, 1, 2, &
  1, 1, 1,    -3, 1, 2, &
  1, 1,-1,     3,-1, 2, &
  1,-1, 1,     3, 1,-2, &
 -1, 1, 1,     3, 2, 1, &
  1, 1, 1,    -3, 2, 1, &
  1, 1,-1,     3,-2, 1, &
  1,-1, 1,     3, 2,-1  &
   /),(/3+3,lattice_bcc_Nslip/))

! twin system <111>{112}
! MISSING: not implemented yet -- now dummy copy from fcc !!
 real(pReal), dimension(3+3,lattice_bcc_Ntwin), parameter :: lattice_bcc_systemTwin = &
 reshape((/&
! Twin system <112>{111}  Sorted according to Eisenlohr & Hantcherli
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
   /),(/3+3,lattice_bcc_Ntwin/))

 real(pReal), dimension(lattice_bcc_Ntwin), parameter :: lattice_bcc_shearTwin = &
 reshape((/&
! Twin system {111}<112>  just a dummy
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123, &
 0.123  &
   /),(/lattice_bcc_Ntwin/))

!*** slip--slip interactions for BCC structures (2) ***
 integer(pInt), target, dimension(lattice_bcc_Nslip,lattice_bcc_Nslip) :: lattice_bcc_interactionSlipSlip = &
 reshape((/&
 1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1  &
   /),(/lattice_bcc_Nslip,lattice_bcc_Nslip/))

!*** slip--twin interactions for BCC structures (2) ***
! MISSING: not implemented yet
 integer(pInt), target, dimension(lattice_bcc_Ntwin,lattice_bcc_Nslip) :: lattice_bcc_interactionSlipTwin = &
 reshape((/&
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0  &
   /),(/lattice_bcc_Ntwin,lattice_bcc_Nslip/))


!*** twin--slip interactions for BCC structures (2) ***
! MISSING: not implemented yet
 integer(pInt), target, dimension(lattice_bcc_Nslip,lattice_bcc_Ntwin) :: lattice_bcc_interactionTwinSlip = &
 reshape((/&
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  &
   /),(/lattice_bcc_Nslip,lattice_bcc_Ntwin/))

!*** twin-twin interactions for BCC structures (2) ***
! MISSING: not implemented yet
 integer(pInt), target, dimension(lattice_bcc_Ntwin,lattice_bcc_Ntwin) :: lattice_bcc_interactionTwinTwin = &
 reshape((/&
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0, &
 0,0,0,0,0,0,0,0,0,0,0,0  &
   /),(/lattice_bcc_Ntwin,lattice_bcc_Ntwin/))
 

!============================== hex (3+) =================================

 integer(pInt), parameter, dimension(lattice_maxNslipFamily) :: lattice_hex_NslipSystem = (/ 3, 3, 6,12/)
 integer(pInt), parameter, dimension(lattice_maxNtwinFamily) :: lattice_hex_NtwinSystem = (/ 6, 6, 6, 6/)
 integer(pInt), parameter :: lattice_hex_Nslip = 24                                       ! sum(lattice_hex_NslipSystem)
 integer(pInt), parameter :: lattice_hex_Ntwin = 24                                       ! sum(lattice_hex_NtwinSystem)
 integer(pInt) ::            lattice_hex_Nstructure = 0_pInt

 !* sorted by YJ.Ro and Philip
 real(pReal), dimension(4+4,lattice_hex_Nslip), parameter :: lattice_hex_systemSlip = &
 reshape((/&
! Basal systems <1120>{0001} (independent of c/a-ratio, Bravais notation (4 coordinate base))
  2, -1, -1,  0,     0,  0,  0,  1, &
  1,  1, -2,  0,     0,  0,  0,  1, &
 -1,  2, -1,  0,     0,  0,  0,  1, &
! 1st type prismatic systems <1120>{1010}  (independent of c/a-ratio)
  2, -1, -1,  0,     0, -1,  1,  0, &
  1,  1,  2,  0,     1, -1,  0,  0, &
 -1,  2, -1,  0,     1,  0, -1,  0, &
! 1st type 1st order pyramidal systems <1120>{1011}
  2, -1, -1,  0,     0, -1,  1,  1, &
  1,  1, -2,  0,     1, -1,  0,  1, &
 -1,  2, -1,  0,     1,  0, -1,  1, &
 -2,  1,  1,  0,     0,  1, -1,  1, &
 -1, -1,  2,  0,    -1,  1,  0,  1, &
  1, -2,  1,  0,    -1,  0,  1,  1, &
! pyramidal system: c+a slip <2113>{1011} -- plane normals depend on the c/a-ratio
  2, -1, -1, -3,     1, -1,  0,  1, &
  1,  1, -2, -3,     1,  0, -1,  1, &
 -1,  2, -1, -3,     1, -1,  0,  1, &
 -2,  1,  1, -3,     1,  0, -1,  1, &
 -1, -1,  2, -3,     0,  1, -1,  1, &
  1, -2,  1, -3,     0, -1,  1,  1, &
 -2,  1,  1, -3,    -1,  0,  1,  1, &
 -1, -1,  2, -3,     0, -1,  1,  1, &
  1, -2,  1, -3,     0,  1, -1,  1, &
  2, -1, -1, -3,    -1,  1,  0,  1, &
  1,  1, -2, -3,    -1,  0,  1,  1, &
 -1,  2, -1, -3,    -1,  1,  0,  1  &
   /),(/4+4,lattice_hex_Nslip/))

 real(pReal), dimension(4+4,lattice_hex_Ntwin), parameter :: lattice_hex_systemTwin = &
 reshape((/&
  0,  1, -1,  1,     0, -1,  1,  2, & ! <-10.1>{10.2} shear = (3-(c/a)^2)/(sqrt(3) c/a)
 -1,  1,  0,  1,     1, -1,  0,  2, &
 -1,  0,  1,  1,     1,  0, -1,  2, & !!
  0, -1,  1,  1,     0,  1, -1,  2, &
  1, -1,  0,  1,    -1,  1,  0,  2, &
  1,  0, -1,  1,    -1,  0,  1,  2, &
  2, -1, -1, -3,     2, -1, -1,  2, & ! <11.-3>{11.2} shear = 2((c/a)^2-2)/(3 c/a)
  1,  1, -2, -3,     1,  1, -2,  2, & !!
 -1,  2, -1, -3,    -1,  2, -1,  2, &
 -2,  1,  1, -3,    -2,  1,  1,  2, &
 -1, -1,  2, -3,    -1, -1,  2,  2, &
  1, -2,  1, -3,     1, -2,  1,  2, &
 -2,  1,  1,  6,     2, -1, -1,  1, & ! <-1-1.6>{11.1} shear = 1/(c/a)
 -1, -1,  2,  6,     1,  1, -2,  1, & !!
  1, -2,  1,  6,    -1,  2, -1,  1, &
  2, -1, -1,  6,    -2,  1,  1,  1, &
  1,  1, -2,  6,    -1, -1,  2,  1, &
 -1,  2, -1,  6,     1, -2,  1,  1, &
  1,  0, -1, -2,     1,  0, -1,  1, & !! <10.-2>{10.1} shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
 -1,  0,  1, -2,    -1,  0,  1,  1, &
  0,  1, -1, -2,     0,  1, -1,  1, &
  0, -1,  1, -2,     0, -1,  1,  1, &
  1, -1,  0, -2,     1, -1,  0,  1, &
 -1,  1,  0, -2,    -1,  1,  0,  1  &
   /),(/4+4,lattice_hex_Ntwin/))                 !* Sort? Numbering of twin system follows Prof. Tom Bieler's scheme (to be consistent with his work); but numbering in data was restarted from 1 &

 integer(pInt), dimension(lattice_hex_Ntwin), parameter :: lattice_hex_shearTwin = &    ! indicator to formula further below
 reshape((/&
 1, &  ! {10.2}<-10.1>
 1, &
 1, &
 1, &
 1, &
 1, &
 2, &  ! {11.2}<11.-3>
 2, &
 2, &
 2, &
 2, &
 2, &
 3, &  ! {11.1}<-1-1.6>
 3, &
 3, &
 3, &
 3, &
 3, &
 4, &  ! {10.1}<10.-2>
 4, &
 4, &
 4, &
 4, &
 4  &
   /),(/lattice_hex_Ntwin/))

!* four different interaction type matrix
 !* 1. slip-slip interaction - 20 types
 !* 2. slip-twin interaction - 16 types
 !* 3. twin-twin interaction - 20 types
 !* 4. twin-slip interaction - 16 types
   
 integer(pInt), target, dimension(lattice_hex_Nslip,lattice_hex_Nslip) :: lattice_hex_interactionSlipSlip = &
 reshape((/&
  1, 5, 5,   9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14,14,14,14,14,14,14, &
  5, 1, 5,   9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14,14,14,14,14,14,14, &
  5, 5, 1,   9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14,14,14,14,14,14,14, &
!
 15,15,15,   2, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13,13,13,13,13,13,13, &
 15,15,15,   6, 2, 6,  10,10,10,10,10,10,  13,13,13,13,13,13,13,13,13,13,13,13, &
 15,15,15,   6, 6, 2,  10,10,10,10,10,10,  13,13,13,13,13,13,13,13,13,13,13,13, &
!
 18,18,18,  16,16,16,   3, 7, 7, 7, 7, 7,  11,11,11,11,11,11,11,11,11,11,11,11, &
 18,18,18,  16,16,16,   7, 3, 7, 7, 7, 7,  11,11,11,11,11,11,11,11,11,11,11,11, &
 18,18,18,  16,16,16,   7, 7, 3, 7, 7, 7,  11,11,11,11,11,11,11,11,11,11,11,11, &
 18,18,18,  16,16,16,   7, 7, 7, 3, 7, 7,  11,11,11,11,11,11,11,11,11,11,11,11, &
 18,18,18,  16,16,16,   7, 7, 7, 7, 3, 7,  11,11,11,11,11,11,11,11,11,11,11,11, &
 18,18,18,  16,16,16,   7, 7, 7, 7, 7, 3,  11,11,11,11,11,11,11,11,11,11,11,11, &
!
 20,20,20,  19,19,19,  17,17,17,17,17,17,   4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, &
 20,20,20,  19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4  &
   /),(/lattice_hex_Nslip,lattice_hex_Nslip/))
  
!* isotropic interaction at the moment
 integer(pInt), target, dimension(lattice_hex_Ntwin,lattice_hex_Nslip) :: lattice_hex_interactionSlipTwin = &
 reshape((/&
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
 13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
 13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
 13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
 13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
 13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
 13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16  &
   /),(/lattice_hex_Ntwin,lattice_hex_Nslip/))  

 !* isotropic interaction at the moment
 integer(pInt), target, dimension(lattice_hex_Nslip,lattice_hex_Ntwin) :: lattice_hex_interactionTwinSlip = &
 reshape((/&
  1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13, & ! --> slip
  1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13, & ! |
  1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13, & ! |
  1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13, & ! v
  1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13, & ! twin
  1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13, &
!
  2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14, &
  2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14, &
  2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14, &
  2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14, &
  2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14, &
  2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14, &
!
  3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15, &
  3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15, &
  3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15, &
  3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15, &
  3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15, &
  3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15, &
!
  4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16, &
  4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16, &
  4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16, &
  4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16, &
  4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16, &
  4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16  &
   /),(/lattice_hex_Nslip,lattice_hex_Ntwin/))  


 integer(pInt), target, dimension(lattice_hex_Ntwin,lattice_hex_Ntwin) :: lattice_hex_interactionTwinTwin = &
 reshape((/&
  1, 5, 5, 5, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
  5, 1, 5, 5, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
  5, 5, 1, 5, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
  5, 5, 5, 1, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
  5, 5, 5, 5, 1, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
  5, 5, 5, 5, 5, 1,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
!
 15,15,15,15,15,15,   2, 6, 6, 6, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
 15,15,15,15,15,15,   6, 2, 6, 6, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
 15,15,15,15,15,15,   6, 6, 2, 6, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
 15,15,15,15,15,15,   6, 6, 6, 2, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
 15,15,15,15,15,15,   6, 6, 6, 6, 2, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
 15,15,15,15,15,15,   6, 6, 6, 6, 6, 2,  10,10,10,10,10,10,  13,13,13,13,13,13, &
!
 18,18,18,18,18,18,  16,16,16,16,16,16,   3, 7, 7, 7, 7, 7,  11,11,11,11,11,11, &
 18,18,18,18,18,18,  16,16,16,16,16,16,   7, 3, 7, 7, 7, 7,  11,11,11,11,11,11, &
 18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 3, 7, 7, 7,  11,11,11,11,11,11, &
 18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 7, 3, 7, 7,  11,11,11,11,11,11, &
 18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 7, 7, 3, 7,  11,11,11,11,11,11, &
 18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 7, 7, 7, 3,  11,11,11,11,11,11, &
!
 20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   4, 8, 8, 8, 8, 8, &
 20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 4, 8, 8, 8, 8, &
 20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 4, 8, 8, 8, &
 20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 8, 4, 8, 8, &
 20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 4, 8, &
 20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 4  &
   /),(/lattice_hex_Ntwin,lattice_hex_Ntwin/))  


CONTAINS
!****************************************
!* - lattice_init
!* - lattice_initializeStructure
!****************************************


subroutine lattice_init()
!**************************************
!*      Module initialization         *
!**************************************
 use IO, only: IO_open_file,IO_countSections,IO_countTagInPart,IO_error
 use material, only: material_configfile,material_partPhase
 implicit none
 
 integer(pInt), parameter :: fileunit = 200
 integer(pInt) i,Nsections

 write(6,*)
 write(6,*) '<<<+-  lattice init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)

 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error(100) ! cannot open config file
 Nsections = IO_countSections(fileunit,material_partPhase)
 lattice_Nstructure = 2_pInt + sum(IO_countTagInPart(fileunit,material_partPhase,'covera_ratio',Nsections)) ! fcc + bcc + all hex
! lattice_Nstructure = Nsections + 2_pInt                                                ! most conservative assumption
 close(fileunit)

 write(6,'(a16,x,i5)') '# sections:',Nsections
 write(6,'(a16,x,i5)') '# structures:',lattice_Nstructure
 write(6,*)

 allocate(lattice_Sslip(3,3,lattice_maxNslip,lattice_Nstructure)); lattice_Sslip   = 0.0_pReal
 allocate(lattice_Sslip_v(6,lattice_maxNslip,lattice_Nstructure)); lattice_Sslip_v = 0.0_pReal
 allocate(lattice_sd(3,lattice_maxNslip,lattice_Nstructure));      lattice_sd      = 0.0_pReal
 allocate(lattice_st(3,lattice_maxNslip,lattice_Nstructure));      lattice_st      = 0.0_pReal
 allocate(lattice_sn(3,lattice_maxNslip,lattice_Nstructure));      lattice_sn      = 0.0_pReal

 allocate(lattice_Qtwin(3,3,lattice_maxNtwin,lattice_Nstructure)); lattice_Qtwin   = 0.0_pReal
 allocate(lattice_Stwin(3,3,lattice_maxNtwin,lattice_Nstructure)); lattice_Stwin   = 0.0_pReal
 allocate(lattice_Stwin_v(6,lattice_maxNtwin,lattice_Nstructure)); lattice_Stwin_v = 0.0_pReal
 allocate(lattice_td(3,lattice_maxNtwin,lattice_Nstructure));      lattice_td      = 0.0_pReal
 allocate(lattice_tt(3,lattice_maxNtwin,lattice_Nstructure));      lattice_tt      = 0.0_pReal
 allocate(lattice_tn(3,lattice_maxNtwin,lattice_Nstructure));      lattice_tn      = 0.0_pReal

 allocate(lattice_shearTwin(lattice_maxNtwin,lattice_Nstructure)); lattice_shearTwin = 0.0_pReal

 allocate(lattice_NslipSystem(lattice_maxNslipFamily,lattice_Nstructure)); lattice_NslipSystem = 0.0_pReal
 allocate(lattice_NtwinSystem(lattice_maxNtwinFamily,lattice_Nstructure)); lattice_NtwinSystem = 0.0_pReal

 allocate(lattice_interactionSlipSlip(lattice_maxNslip,lattice_maxNslip,lattice_Nstructure)); lattice_interactionSlipSlip = 0_pInt  ! other:me
 allocate(lattice_interactionSlipTwin(lattice_maxNtwin,lattice_maxNslip,lattice_Nstructure)); lattice_interactionSlipTwin = 0_pInt  ! other:me
 allocate(lattice_interactionTwinSlip(lattice_maxNslip,lattice_maxNtwin,lattice_Nstructure)); lattice_interactionTwinSlip = 0_pInt  ! other:me
 allocate(lattice_interactionTwinTwin(lattice_maxNtwin,lattice_maxNtwin,lattice_Nstructure)); lattice_interactionTwinTwin = 0_pInt  ! other:me

end subroutine


function lattice_initializeStructure(struct,CoverA)
!**************************************
!*   Calculation of Schmid            *
!*   matrices, etc.                   *
!**************************************
 use prec, only: pReal,pInt
 use math
 use IO, only: IO_error
 implicit none

 character(len=*) struct
 real(pReal) CoverA
 real(pReal), dimension(3,lattice_maxNslip) :: sd = 0.0_pReal, &
                                               sn = 0.0_pReal, &
                                               st = 0.0_pReal
 real(pReal), dimension(3,lattice_maxNtwin) :: td = 0.0_pReal, &
                                               tn = 0.0_pReal, &
                                               tt = 0.0_pReal
 real(pReal), dimension(lattice_maxNtwin) ::   ts = 0.0_pReal
 real(pReal), dimension(3) ::                  hex_d = 0.0_pReal, &
                                               hex_n = 0.0_pReal
 integer(pInt), dimension(lattice_maxNslipFamily) :: myNslipSystem = 0_pInt
 integer(pInt), dimension(lattice_maxNtwinFamily) :: myNtwinSystem = 0_pInt
 integer(pInt) :: i,myNslip,myNtwin,myStructure = 0_pInt
 logical :: processMe

 integer(pInt) lattice_initializeStructure

 processMe = .false.

 select case(struct(1:3))                          ! check first three chars of structure name
   case ('fcc')
     myStructure = 1_pInt
     myNslipSystem = lattice_fcc_NslipSystem       ! size of slip system families
     myNtwinSystem = lattice_fcc_NtwinSystem       ! size of twin system families
     myNslip = lattice_fcc_Nslip                   ! overall number of slip systems
     myNtwin = lattice_fcc_Ntwin                   ! overall number of twin systems
     lattice_fcc_Nstructure = lattice_fcc_Nstructure + 1_pInt    ! count fcc instances
     if (lattice_fcc_Nstructure == 1_pInt) then    ! me is first fcc structure
       processMe = .true.
       do i = 1,myNslip                            ! calculate slip system vectors
         sd(:,i) = lattice_fcc_systemSlip(1:3,i)/dsqrt(math_mul3x3(lattice_fcc_systemSlip(1:3,i),lattice_fcc_systemSlip(1:3,i)))
         sn(:,i) = lattice_fcc_systemSlip(4:6,i)/dsqrt(math_mul3x3(lattice_fcc_systemSlip(4:6,i),lattice_fcc_systemSlip(4:6,i)))
         st(:,i) = math_vectorproduct(sd(:,i),sn(:,i))
       enddo
       do i = 1,myNtwin                            ! calculate twin system vectors and (assign) shears
         td(:,i) = lattice_fcc_systemTwin(1:3,i)/dsqrt(math_mul3x3(lattice_fcc_systemTwin(1:3,i),lattice_fcc_systemTwin(1:3,i)))
         tn(:,i) = lattice_fcc_systemTwin(4:6,i)/dsqrt(math_mul3x3(lattice_fcc_systemTwin(4:6,i),lattice_fcc_systemTwin(4:6,i)))
         tt(:,i) = math_vectorproduct(td(:,i),tn(:,i))
         ts(i)   = lattice_fcc_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_fcc_interactionSlipSlip
       interactionSlipTwin => lattice_fcc_interactionSlipTwin
       interactionTwinSlip => lattice_fcc_interactionTwinSlip
       interactionTwinTwin => lattice_fcc_interactionTwinTwin
     endif
     
   case ('bcc')
     myStructure = 2_pInt
     myNslipSystem = lattice_bcc_NslipSystem       ! size of slip system families
     myNtwinSystem = lattice_bcc_NtwinSystem       ! size of twin system families
     myNslip = lattice_bcc_Nslip                   ! overall number of slip systems
     myNtwin = lattice_bcc_Ntwin                   ! overall number of twin systems
     lattice_bcc_Nstructure = lattice_bcc_Nstructure + 1_pInt    ! count bcc instances
     if (lattice_bcc_Nstructure == 1_pInt) then    ! me is first bcc structure
       processMe = .true.
       do i = 1,myNslip                            ! calculate slip system vectors
         sd(:,i) = lattice_bcc_systemSlip(1:3,i)/dsqrt(math_mul3x3(lattice_bcc_systemSlip(1:3,i),lattice_bcc_systemSlip(1:3,i)))
         sn(:,i) = lattice_bcc_systemSlip(4:6,i)/dsqrt(math_mul3x3(lattice_bcc_systemSlip(4:6,i),lattice_bcc_systemSlip(4:6,i)))
         st(:,i) = math_vectorproduct(sd(:,i),sn(:,i))
       enddo
       do i = 1,myNtwin                            ! calculate twin system vectors and (assign) shears
         td(:,i) = lattice_bcc_systemTwin(1:3,i)/dsqrt(math_mul3x3(lattice_bcc_systemTwin(1:3,i),lattice_bcc_systemTwin(1:3,i)))
         tn(:,i) = lattice_bcc_systemTwin(4:6,i)/dsqrt(math_mul3x3(lattice_bcc_systemTwin(4:6,i),lattice_bcc_systemTwin(4:6,i)))
         tt(:,i) = math_vectorproduct(td(:,i),tn(:,i))
         ts(i)   = lattice_bcc_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_bcc_interactionSlipSlip
       interactionSlipTwin => lattice_bcc_interactionSlipTwin
       interactionTwinSlip => lattice_bcc_interactionTwinSlip
       interactionTwinTwin => lattice_bcc_interactionTwinTwin
     endif
     
   case ('hex')
     if (CoverA >= 1.0_pReal) then                 ! checking physical significance of c/a
       lattice_hex_Nstructure = lattice_hex_Nstructure + 1_pInt  ! count instances of hex structures
       myStructure = 2_pInt + lattice_hex_Nstructure             ! 3,4,5,.. for hex
       myNslipSystem = lattice_hex_NslipSystem     ! size of slip system families
       myNtwinSystem = lattice_hex_NtwinSystem     ! size of twin system families
       myNslip = lattice_hex_Nslip                 ! overall number of slip systems
       myNtwin = lattice_hex_Ntwin                 ! overall number of twin systems
       processMe = .true.
! converting from 4 axes coordinate system (a1=a2=a3=c) to ortho-hexgonal system (a, b, c)
       do i = 1,myNslip
         hex_d(1) =  lattice_hex_systemSlip(1,i)*1.5_pReal ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]
         hex_d(2) = (lattice_hex_systemSlip(1,i)+2.0_pReal*lattice_hex_systemSlip(2,i))*(0.5_pReal*dsqrt(3.0_pReal))
         hex_d(3) =  lattice_hex_systemSlip(4,i)*CoverA
         hex_n(1) =  lattice_hex_systemSlip(5,i)           ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))
         hex_n(2) = (lattice_hex_systemSlip(5,i)+2.0_pReal*lattice_hex_systemSlip(6,i))/dsqrt(3.0_pReal)
         hex_n(3) =  lattice_hex_systemSlip(8,i)/CoverA

         sd(:,i) = hex_d/dsqrt(math_mul3x3(hex_d,hex_d))
         sn(:,i) = hex_n/dsqrt(math_mul3x3(hex_n,hex_n))
         st(:,i) = math_vectorproduct(sd(:,i),sn(:,i))
       enddo
       do i = 1,myNtwin
         hex_d(1) =  lattice_hex_systemTwin(1,i)*1.5_pReal
         hex_d(2) = (lattice_hex_systemTwin(1,i)+2.0_pReal*lattice_hex_systemTwin(2,i))*(0.5_pReal*dsqrt(3.0_pReal))
         hex_d(3) =  lattice_hex_systemTwin(4,i)*CoverA
         hex_n(1) =  lattice_hex_systemTwin(5,i)
         hex_n(2) = (lattice_hex_systemTwin(5,i)+2.0_pReal*lattice_hex_systemTwin(6,i))/dsqrt(3.0_pReal)
         hex_n(3) =  lattice_hex_systemTwin(8,i)/CoverA

         td(:,i) = hex_d/dsqrt(math_mul3x3(hex_d,hex_d))
         tn(:,i) = hex_n/dsqrt(math_mul3x3(hex_n,hex_n))
         tt(:,i) = math_vectorproduct(td(:,i),tn(:,i))

         select case(lattice_hex_shearTwin(i))                                     ! from Christian & Mahajan 1995 p.29
           case (1)                                                                ! {10.2}<-10.1>
                    ts(i) = (3.0_pReal-CoverA*CoverA)/dsqrt(3.0_pReal)/CoverA
           case (2)                                                                ! {11.2}<11.-3>
                    ts(i) = 2.0_pReal*(CoverA*CoverA-2.0_pReal)/3.0_pReal/CoverA
           case (3)                                                                ! {11.1}<-1-1.6>
                    ts(i) = 1.0_pReal/CoverA
           case (4)                                                                ! {10.1}<10.-2>
                    ts(i) = (4.0_pReal*CoverA*CoverA-9.0_pReal)/4.0_pReal/dsqrt(3.0_pReal)/CoverA
         end select

       enddo
       interactionSlipSlip => lattice_hex_interactionSlipSlip
       interactionSlipTwin => lattice_hex_interactionSlipTwin
       interactionTwinSlip => lattice_hex_interactionTwinSlip
       interactionTwinTwin => lattice_hex_interactionTwinTwin
     endif
 end select

 if (processMe) then
   if  (myStructure > lattice_Nstructure) &
     call IO_error(666,0,0,0,'structure index too large')        ! check for memory leakage
   do i = 1,myNslip                                              ! store slip system vectors and Schmid matrix for my structure
     lattice_sd(:,i,myStructure) = sd(:,i)
     lattice_st(:,i,myStructure) = st(:,i)
     lattice_sn(:,i,myStructure) = sn(:,i)
     lattice_Sslip(:,:,i,myStructure) = math_tensorproduct(sd(:,i),sn(:,i))
     lattice_Sslip_v(:,i,myStructure) = math_Mandel33to6(math_symmetric3x3(lattice_Sslip(:,:,i,myStructure)))
   enddo
   do i = 1,myNtwin                                              ! store twin system vectors and Schmid plus rotation matrix for my structure
     lattice_td(:,i,myStructure) = td(:,i)
     lattice_tt(:,i,myStructure) = tt(:,i)
     lattice_tn(:,i,myStructure) = tn(:,i)
     lattice_Stwin(:,:,i,myStructure) = math_tensorproduct(td(:,i),tn(:,i))
     lattice_Stwin_v(:,i,myStructure) = math_Mandel33to6(math_symmetric3x3(lattice_Stwin(:,:,i,myStructure)))
     lattice_Qtwin(:,:,i,myStructure) = math_RodrigToR(tn(:,i),180.0_pReal*inRad)
     lattice_shearTwin(i,myStructure) = ts(i)
   enddo
   lattice_NslipSystem(1:lattice_maxNslipFamily,myStructure) = myNslipSystem                                ! number of slip systems in each family
   lattice_NtwinSystem(1:lattice_maxNtwinFamily,myStructure) = myNtwinSystem                                ! number of twin systems in each family
   lattice_interactionSlipSlip(1:myNslip,1:myNslip,myStructure) = interactionSlipSlip(1:myNslip,1:myNslip)
   lattice_interactionSlipTwin(1:myNtwin,1:myNslip,myStructure) = interactionSlipTwin(1:myNtwin,1:myNslip)
   lattice_interactionTwinSlip(1:myNslip,1:myNtwin,myStructure) = interactionTwinSlip(1:myNslip,1:myNtwin)
   lattice_interactionTwinTwin(1:myNtwin,1:myNtwin,myStructure) = interactionTwinTwin(1:myNtwin,1:myNtwin)
 endif

 lattice_initializeStructure = myStructure        ! report my structure index back

end function


END MODULE
