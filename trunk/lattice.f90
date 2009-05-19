
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
integer(pInt), parameter :: lattice_maxNslip = 48                 ! max # of slip systems over lattice structures
integer(pInt), parameter :: lattice_maxNtwin = 24                 ! max # of twin systems over lattice structures

integer(pInt), pointer, dimension(:,:) :: interactionSlipSlip, &
                                          interactionSlipTwin, &
                                          interactionTwinTwin, &
										  interactionTwinSlip

! Schmid matrices, normal, shear direction and nxd of slip systems
real(pReal), allocatable, dimension(:,:,:,:) :: lattice_Sslip
real(pReal), allocatable, dimension(:,:,:)   :: lattice_Sslip_v
real(pReal), allocatable, dimension(:,:,:)   :: lattice_sn
real(pReal), allocatable, dimension(:,:,:)   :: lattice_sd
real(pReal), allocatable, dimension(:,:,:)   :: lattice_st

! Rotation and Schmid matrices, normal, shear direction and nxd of twin systems
real(pReal), allocatable, dimension(:,:,:,:) :: lattice_Qtwin
real(pReal), allocatable, dimension(:,:,:,:) :: lattice_Stwin
real(pReal), allocatable, dimension(:,:,:)   :: lattice_Stwin_v
real(pReal), allocatable, dimension(:,:,:)   :: lattice_tn
real(pReal), allocatable, dimension(:,:,:)   :: lattice_td
real(pReal), allocatable, dimension(:,:,:)   :: lattice_tt
real(pReal), allocatable, dimension(:,:)     :: lattice_shearTwin

integer(pInt), allocatable, dimension(:,:,:) :: lattice_interactionSlipSlip
integer(pInt), allocatable, dimension(:,:,:) :: lattice_interactionSlipTwin
integer(pInt), allocatable, dimension(:,:,:) :: lattice_interactionTwinTwin
integer(pInt), allocatable, dimension(:,:,:) :: lattice_interactionTwinSlip



!============================== fcc (1) =================================

 integer(pInt), parameter :: lattice_fcc_Nslip = 12_pInt
 integer(pInt), parameter :: lattice_fcc_Ntwin = 12_pInt
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

 integer(pInt), target, dimension(lattice_fcc_Nslip,lattice_fcc_Ntwin) :: lattice_fcc_interactionSlipTwin = &
 reshape((/&
 0,0,0,1,1,0,0,1,1,1,0,1, &
 0,0,0,1,1,0,0,1,1,1,0,1, &
 0,0,0,1,1,0,0,1,1,1,0,1, &
 1,1,0,0,0,0,1,0,1,0,1,1, &
 1,1,0,0,0,0,1,0,1,0,1,1, &
 1,1,0,0,0,0,1,0,1,0,1,1, &
 0,1,1,1,0,1,0,0,0,1,1,0, &
 0,1,1,1,0,1,0,0,0,1,1,0, &
 0,1,1,1,0,1,0,0,0,1,1,0, &
 1,0,1,0,1,1,1,1,0,0,0,0, &
 1,0,1,0,1,1,1,1,0,0,0,0, &
 1,0,1,0,1,1,1,1,0,0,0,0  &
   /),(/lattice_fcc_Nslip,lattice_fcc_Ntwin/))

 integer(pInt), target, dimension(lattice_fcc_Ntwin,lattice_fcc_Ntwin) :: lattice_fcc_interactionTwinTwin = &
 reshape((/&
 0,0,0,1,1,1,1,1,1,1,1,1, &
 0,0,0,1,1,1,1,1,1,1,1,1, &
 0,0,0,1,1,1,1,1,1,1,1,1, &
 1,1,1,0,0,0,1,1,1,1,1,1, &
 1,1,1,0,0,0,1,1,1,1,1,1, &
 1,1,1,0,0,0,1,1,1,1,1,1, &
 1,1,1,1,1,1,0,0,0,1,1,1, &
 1,1,1,1,1,1,0,0,0,1,1,1, &
 1,1,1,1,1,1,0,0,0,1,1,1, &
 1,1,1,1,1,1,1,1,1,0,0,0, &
 1,1,1,1,1,1,1,1,1,0,0,0, &
 1,1,1,1,1,1,1,1,1,0,0,0  &
   /),(/lattice_fcc_Ntwin,lattice_fcc_Ntwin/))


!============================== bcc (2) =================================

 integer(pInt), parameter :: lattice_bcc_Nslip = 48_pInt
 integer(pInt), parameter :: lattice_bcc_Ntwin = 12_pInt
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

!*** Slip-Slip interactions for BCC structures (2) ***
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

!*** Slip-twin interactions for BCC structures (2) ***
! MISSING: not implemented yet
 integer(pInt), target, dimension(lattice_bcc_Nslip,lattice_bcc_Ntwin) :: lattice_bcc_interactionSlipTwin = &
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

!*** Twin-twin interactions for BCC structures (2) ***
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

 integer(pInt), parameter :: lattice_hex_Nslip = 24_pInt
 integer(pInt), parameter :: lattice_hex_Ntwin = 24_pInt
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
  0,  1, -1,  1,     0, -1,  1,  2, & ! <1011>{1012} Twin: shear 0.169 -1.26 compression
 -1,  1,  0,  1,     1, -1,  0,  2, &
 -1,  0,  1,  1,     1,  0, -1,  2, &
  0, -1,  1,  1,     0,  1, -1,  2, &
  1, -1,  0,  1,    -1,  1,  0,  2, &
  1,  0, -1,  1,    -1,  0,  1,  2, &
  2, -1, -1, -3,     2, -1, -1,  2, & ! <211-2>{2112} Twin: shear 0.224 1.19 tension
  1,  1, -2, -3,     1,  1, -2,  2, &
 -1,  2, -1, -3,    -1,  2, -1,  2, &
 -2,  1,  1, -3,    -2,  1,  1,  2, &
 -1, -1,  2, -3,    -1, -1,  2,  2, &
  1, -2,  1, -3,     1, -2,  1,  2, &
 -2,  1,  1,  6,     2, -1, -1,  1, & ! <211-6>{2111} Twin: shear 0.628 -0.39 compression
 -1, -1,  2,  6,     1,  1, -2,  1, &
  1, -2,  1,  6,    -1,  2, -1,  1, &
  2, -1, -1,  6,    -2,  1,  1,  1, &
  1,  1, -2,  6,    -1, -1,  2,  1, &
 -1,  2, -1,  6,     1, -2,  1,  1, &
  1,  0, -1, -2,     1,  0, -1,  1, & ! <101-2>{1011} Twin: shear 0.103 1.09 tension
 -1,  0,  1, -2,    -1,  0,  1,  1, &
  0,  1, -1, -2,     0,  1, -1,  1, &
  0, -1,  1, -2,     0, -1,  1,  1, &
  1, -1,  0, -2,     1, -1,  0,  1, &
 -1,  1,  0, -2,    -1,  1,  0,  1  &
   /),(/4+4,lattice_hex_Ntwin/))                 !* Sort? Numbering of twin system follows Prof. Tom Bieler's scheme (to be consistent with his work); but numbering in data was restarted from 1 &

 real(pReal), dimension(lattice_hex_Ntwin), parameter :: lattice_hex_shearTwin = &
 reshape((/&
 0.169, &  ! <1011>{1012} Twin: shear 0.169 -1.26 compression
 0.169, &
 0.169, &
 0.169, &
 0.169, &
 0.169, &
 0.224, &  ! <211-2>{2112} Twin: shear 0.224 1.19 tension
 0.224, &
 0.224, &
 0.224, &
 0.224, &
 0.224, &
 0.628, &  ! <211-6>{2111} Twin: shear 0.628 -0.39 compression
 0.628, &
 0.628, &
 0.628, &
 0.628, &
 0.628, &
 0.103, &  ! <101-2>{1011} Twin: shear 0.103 1.09 tension
 0.103, &
 0.103, &
 0.103, &
 0.103, &
 0.103  &
   /),(/lattice_hex_Ntwin/))

!* four different interaction type matrix
 !* 1. slip-slip interaction
    !- wihtin the same slip familiy, self or latent hardening; 1 (self), 2 ~ 5 (latent for each slip; bas, pri, pyr_a, pyr_c+a in sequence)
	!- interaction between different slip families; 6 ~ 11
 !* 2. slip-twin interaction: current scheme -- indirect effect of twin volume
	!- all 1 
 !* 3. twin-twin interaction
    !- wihtin the same slip familiy, self or latent hardening; 1 (self), 2 ~ 5
	!- interaction between different slip families; 6 ~ 11
 !* 4. twin-slip interaction
	!- T1 interacting slip; 1 ~ 4, C1 interacting slip; 5 ~ 8
	!- T2 interacting slip; 9 ~ 12, C1 interacting slip; 13 ~ 16
   
 integer(pInt), target, dimension(lattice_hex_Nslip,lattice_hex_Nslip) :: lattice_hex_interactionSlipSlip = &
 reshape((/&
 1,2,2,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8, &
 2,1,2,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8, &
 2,2,1,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8, &
 6,6,6,1,3,3,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10, &
 6,6,6,3,1,3,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10, &
 6,6,6,3,3,1,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10, &
 7,7,7,9,9,9,1,4,4,4,4,4,11,11,11,11,11,11,11,11,11,11,11,11, &
 7,7,7,9,9,9,4,1,4,4,4,4,11,11,11,11,11,11,11,11,11,11,11,11, &
 7,7,7,9,9,9,4,4,1,4,4,4,11,11,11,11,11,11,11,11,11,11,11,11, &
 7,7,7,9,9,9,4,4,4,1,4,4,11,11,11,11,11,11,11,11,11,11,11,11, &
 7,7,7,9,9,9,4,4,4,4,1,4,11,11,11,11,11,11,11,11,11,11,11,11, &
 7,7,7,9,9,9,4,4,4,4,4,1,11,11,11,11,11,11,11,11,11,11,11,11, &
 8,8,8,10,10,10,11,11,11,11,11,11,1,5,1,5,5,5,5,5,5,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,1,5,1,5,5,5,5,5,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,1,5,1,5,5,5,5,5,5,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,1,5,1,5,5,5,5,5,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,1,5,5,5,1,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,5,1,5,1,5,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,5,5,1,5,5,5,1,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,5,1,5,1,5,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,1,5,5,5,1,5,5,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,5,5,5,5,5,1,5,1, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,5,5,1,5,5,5,1,5, &
 8,8,8,10,10,10,11,11,11,11,11,11,5,5,5,5,5,5,5,5,5,1,5,1  &
   /),(/lattice_hex_Nslip,lattice_hex_Nslip/))
  
!* isotropic interaction at the moment
 integer(pInt), target, dimension(lattice_hex_Nslip,lattice_hex_Ntwin) :: lattice_hex_interactionSlipTwin = &
 reshape((/&
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1  &
   /),(/lattice_hex_Nslip,lattice_hex_Ntwin/))  


 integer(pInt), target, dimension(lattice_hex_Ntwin,lattice_hex_Ntwin) :: lattice_hex_interactionTwinTwin = &
 reshape((/&
 1,2,2,2,2,2,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8, &
 2,1,2,2,2,2,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8, &
 2,2,1,2,2,2,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8, &
 2,2,2,1,2,2,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8, &
 2,2,2,2,1,2,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8, &
 2,2,2,2,2,1,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8, &
 6,6,6,6,6,6,1,3,3,3,3,3,9,9,9,9,9,9,10,10,10,10,10,10, &
 6,6,6,6,6,6,3,1,3,3,3,3,9,9,9,9,9,9,10,10,10,10,10,10, &
 6,6,6,6,6,6,3,3,1,3,3,3,9,9,9,9,9,9,10,10,10,10,10,10, &
 6,6,6,6,6,6,3,3,3,1,3,3,9,9,9,9,9,9,10,10,10,10,10,10, &
 6,6,6,6,6,6,3,3,3,3,1,3,9,9,9,9,9,9,10,10,10,10,10,10, &
 6,6,6,6,6,6,3,3,3,3,3,1,9,9,9,9,9,9,10,10,10,10,10,10, &
 7,7,7,7,7,7,9,9,9,9,9,9,1,4,4,4,4,4,11,11,11,11,11,11, &
 7,7,7,7,7,7,9,9,9,9,9,9,4,1,4,4,4,4,11,11,11,11,11,11, &
 7,7,7,7,7,7,9,9,9,9,9,9,4,4,1,4,4,4,11,11,11,11,11,11, &
 7,7,7,7,7,7,9,9,9,9,9,9,4,4,4,1,4,4,11,11,11,11,11,11, &
 7,7,7,7,7,7,9,9,9,9,9,9,4,4,4,4,1,4,11,11,11,11,11,11, &
 7,7,7,7,7,7,9,9,9,9,9,9,4,4,4,4,4,1,11,11,11,11,11,11, &
 8,8,8,8,8,8,10,10,10,10,10,10,11,11,11,11,11,11,1,5,5,5,5,5, &
 8,8,8,8,8,8,10,10,10,10,10,10,11,11,11,11,11,11,5,1,5,5,5,5, &
 8,8,8,8,8,8,10,10,10,10,10,10,11,11,11,11,11,11,5,5,1,5,5,5, &
 8,8,8,8,8,8,10,10,10,10,10,10,11,11,11,11,11,11,5,5,5,1,5,5, &
 8,8,8,8,8,8,10,10,10,10,10,10,11,11,11,11,11,11,5,5,5,5,1,5, &
 8,8,8,8,8,8,10,10,10,10,10,10,11,11,11,11,11,11,5,5,5,5,5,1  &
   /),(/lattice_hex_Ntwin,lattice_hex_Ntwin/))  

 !* isotropic interaction at the moment
 integer(pInt), target, dimension(lattice_hex_Ntwin,lattice_hex_Nslip) :: lattice_hex_interactionTwinSlip = &
 reshape((/&
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1  &
   /),(/lattice_hex_Ntwin,lattice_hex_Nslip/))  


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

 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file
 Nsections = IO_countSections(fileunit,material_partPhase)
 lattice_Nstructure = 2_pInt + sum(IO_countTagInPart(fileunit,material_partPhase,'covera_ratio',Nsections)) ! fcc + bcc + all hex
 close(fileunit)

 allocate(lattice_Sslip(3,3,lattice_maxNslip,lattice_Nstructure)); lattice_Sslip = 0.0_pReal
 allocate(lattice_Sslip_v(6,lattice_maxNslip,lattice_Nstructure)); lattice_Sslip_v = 0.0_pReal
 allocate(lattice_sd(3,lattice_maxNslip,lattice_Nstructure)); lattice_sd = 0.0_pReal
 allocate(lattice_st(3,lattice_maxNslip,lattice_Nstructure)); lattice_st = 0.0_pReal
 allocate(lattice_sn(3,lattice_maxNslip,lattice_Nstructure)); lattice_sn = 0.0_pReal

 allocate(lattice_Qtwin(3,3,lattice_maxNtwin,lattice_Nstructure)); lattice_Qtwin = 0.0_pReal
 allocate(lattice_Stwin(3,3,lattice_maxNtwin,lattice_Nstructure)); lattice_Stwin = 0.0_pReal
 allocate(lattice_Stwin_v(6,lattice_maxNtwin,lattice_Nstructure)); lattice_Stwin_v = 0.0_pReal
 allocate(lattice_td(3,lattice_maxNtwin,lattice_Nstructure)); lattice_td = 0.0_pReal
 allocate(lattice_tt(3,lattice_maxNtwin,lattice_Nstructure)); lattice_tt = 0.0_pReal
 allocate(lattice_tn(3,lattice_maxNtwin,lattice_Nstructure)); lattice_tn = 0.0_pReal
 allocate(lattice_shearTwin(lattice_maxNtwin,lattice_Nstructure)); lattice_shearTwin = 0.0_pReal

 allocate(lattice_interactionSlipSlip(lattice_maxNslip,lattice_maxNslip,lattice_Nstructure)); lattice_interactionSlipSlip = 0_pInt
 allocate(lattice_interactionSlipTwin(lattice_maxNslip,lattice_maxNtwin,lattice_Nstructure)); lattice_interactionSlipTwin = 0_pInt
 allocate(lattice_interactionTwinTwin(lattice_maxNtwin,lattice_maxNtwin,lattice_Nstructure)); lattice_interactionTwinTwin = 0_pInt
 allocate(lattice_interactionTwinTwin(lattice_maxNtwin,lattice_maxNtwin,lattice_Nstructure)); lattice_interactionTwinSlip = 0_pInt

end subroutine


function lattice_initializeStructure(struct,CoverA)
!**************************************
!*   Calculation of Schmid            *
!*   matrices, etc.                   *
!**************************************
 use prec, only: pReal,pInt
 use math
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
 integer(pInt) :: i,myNslip,myNtwin,myStructure = 0_pInt
 logical :: processMe = .false.

 integer(pInt) lattice_initializeStructure

 select case(struct(1:3))                          ! check first three chars of structure name
   case ('fcc')
     myStructure = 1_pInt
     myNslip = lattice_fcc_Nslip
     myNtwin = lattice_fcc_Ntwin
     lattice_fcc_Nstructure = lattice_fcc_Nstructure + 1_pInt
     if (lattice_fcc_Nstructure == 1_pInt) then
       processMe = .true.
       do i = 1,myNslip
         sn(:,i) = lattice_fcc_systemSlip(1:3,i)/dsqrt(math_mul3x3(lattice_fcc_systemSlip(1:3,i),lattice_fcc_systemSlip(1:3,i)))
         sd(:,i) = lattice_fcc_systemSlip(4:6,i)/dsqrt(math_mul3x3(lattice_fcc_systemSlip(4:6,i),lattice_fcc_systemSlip(4:6,i)))
         st(:,i) = math_vectorproduct(sn(:,i),sd(:,i))
       enddo
       do i = 1,myNtwin
         tn(:,i) = lattice_fcc_systemTwin(1:3,i)/dsqrt(math_mul3x3(lattice_fcc_systemTwin(1:3,i),lattice_fcc_systemTwin(1:3,i)))
         td(:,i) = lattice_fcc_systemTwin(4:6,i)/dsqrt(math_mul3x3(lattice_fcc_systemTwin(4:6,i),lattice_fcc_systemTwin(4:6,i)))
         tt(:,i) = math_vectorproduct(tn(:,i),td(:,i))
         ts(i)   = lattice_fcc_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_fcc_interactionSlipSlip
       interactionSlipTwin => lattice_fcc_interactionSlipTwin
       interactionTwinTwin => lattice_fcc_interactionTwinTwin
     endif
     
   case ('bcc')
     myStructure = 2_pInt
     myNslip = lattice_bcc_Nslip
     myNtwin = lattice_bcc_Ntwin
     lattice_bcc_Nstructure = lattice_bcc_Nstructure + 1_pInt
     if (lattice_bcc_Nstructure == 1_pInt) then
       processMe = .true.
       do i = 1,myNslip
         sn(:,i) = lattice_bcc_systemSlip(1:3,i)/dsqrt(math_mul3x3(lattice_bcc_systemSlip(1:3,i),lattice_bcc_systemSlip(1:3,i)))
         sd(:,i) = lattice_bcc_systemSlip(4:6,i)/dsqrt(math_mul3x3(lattice_bcc_systemSlip(4:6,i),lattice_bcc_systemSlip(4:6,i)))
         st(:,i) = math_vectorproduct(sn(:,i),sd(:,i))
       enddo
       do i = 1,myNtwin
         tn(:,i) = lattice_bcc_systemTwin(1:3,i)/dsqrt(math_mul3x3(lattice_bcc_systemTwin(1:3,i),lattice_bcc_systemTwin(1:3,i)))
         td(:,i) = lattice_bcc_systemTwin(4:6,i)/dsqrt(math_mul3x3(lattice_bcc_systemTwin(4:6,i),lattice_bcc_systemTwin(4:6,i)))
         tt(:,i) = math_vectorproduct(tn(:,i),td(:,i))
         ts(i)   = lattice_bcc_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_bcc_interactionSlipSlip
       interactionSlipTwin => lattice_bcc_interactionSlipTwin
       interactionTwinTwin => lattice_bcc_interactionTwinTwin
     endif
     
   case ('hex')
     if (CoverA > 0.0_pReal) then
       lattice_hex_Nstructure = lattice_hex_Nstructure + 1_pInt
       myStructure = 2_pInt + lattice_hex_Nstructure
       myNslip = lattice_hex_Nslip
       myNtwin = lattice_hex_Ntwin
       processMe = .true.
! converting from 4 axes coordinate system (a1=a2=a3=c) to ortho-hexgonal system (a, b, c)
       do i = 1,myNslip
         hex_n(1) =  lattice_hex_systemSlip(1,i)           ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))
         hex_n(2) = (lattice_hex_systemSlip(1,i)+2.0_pReal*lattice_hex_systemSlip(2,i))/dsqrt(3.0_pReal)
         hex_n(3) =  lattice_hex_systemSlip(4,i)/CoverA
         hex_d(1) =  lattice_hex_systemSlip(5,i)*1.5_pReal ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]
         hex_d(2) = (lattice_hex_systemSlip(5,i)+2.0_pReal*lattice_hex_systemSlip(6,i))*(0.5_pReal*dsqrt(3.0_pReal))
         hex_d(3) =  lattice_hex_systemSlip(8,i)*CoverA

         sn(:,i) = hex_n/dsqrt(math_mul3x3(hex_n,hex_n))
         sd(:,i) = hex_d/dsqrt(math_mul3x3(hex_d,hex_d))
         st(:,i) = math_vectorproduct(sn(:,i),sd(:,i))
       enddo
       do i = 1,myNtwin
         hex_n(1) =  lattice_hex_systemTwin(1,i)
         hex_n(2) = (lattice_hex_systemTwin(1,i)+2.0_pReal*lattice_hex_systemTwin(2,i))/dsqrt(3.0_pReal)
         hex_n(3) =  lattice_hex_systemTwin(4,i)/CoverA
         hex_d(1) =  lattice_hex_systemTwin(5,i)*1.5_pReal
         hex_d(2) = (lattice_hex_systemTwin(5,i)+2.0_pReal*lattice_hex_systemTwin(6,i))*(0.5_pReal*dsqrt(3.0_pReal))
         hex_d(3) =  lattice_hex_systemTwin(8,i)*CoverA

         tn(:,i) = hex_n/dsqrt(math_mul3x3(hex_n,hex_n))
         td(:,i) = hex_d/dsqrt(math_mul3x3(hex_d,hex_d))
         tt(:,i) = math_vectorproduct(tn(:,i),td(:,i))
         ts(i)   = lattice_hex_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_hex_interactionSlipSlip
       interactionSlipTwin => lattice_hex_interactionSlipTwin
       interactionTwinTwin => lattice_hex_interactionTwinTwin
       interactionTwinTwin => lattice_hex_interactionTwinSlip
     endif
   end select

 if (processMe) then
   do i = 1,myNslip
     lattice_sd(:,i,myStructure) = sd(:,i)
     lattice_st(:,i,myStructure) = st(:,i)
     lattice_sn(:,i,myStructure) = sn(:,i)
     lattice_Sslip(:,:,i,myStructure) = math_tensorproduct(sd(:,i),sn(:,i))
     lattice_Sslip_v(:,i,myStructure) = math_Mandel33to6(math_symmetric3x3(lattice_Sslip(:,:,i,myStructure)))
   enddo
   do i = 1,myNtwin
     lattice_td(:,i,myStructure) = td(:,i)
     lattice_tt(:,i,myStructure) = tt(:,i)
     lattice_tn(:,i,myStructure) = tn(:,i)
     lattice_Stwin(:,:,i,myStructure) = math_tensorproduct(td(:,i),tn(:,i))
     lattice_Stwin_v(:,i,myStructure) = math_Mandel33to6(math_symmetric3x3(lattice_Stwin(:,:,i,myStructure)))
     lattice_Qtwin(:,:,i,myStructure) = math_RodrigToR(tn(:,i),180.0_pReal*inRad)
     lattice_shearTwin(i,myStructure) = ts(i)
   enddo
   lattice_interactionSlipSlip(1:myNslip,1:myNslip,myStructure) = interactionSlipSlip(1:myNslip,1:myNslip)
   lattice_interactionSlipTwin(1:myNslip,1:myNtwin,myStructure) = interactionSlipTwin(1:myNslip,1:myNtwin)
   lattice_interactionTwinTwin(1:myNtwin,1:myNtwin,myStructure) = interactionTwinTwin(1:myNtwin,1:myNtwin)
   lattice_interactionTwinSlip(1:myNtwin,1:myNslip,myStructure) = interactionTwinSlip(1:myNtwin,1:myNslip)
 endif

 lattice_initializeStructure = myStructure

end function


END MODULE
