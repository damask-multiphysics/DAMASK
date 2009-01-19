
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
!* Number of lattice structures (1-FCC,2-BCC,3-HCP)
integer(pInt), parameter :: lattice_MaxLatticeStructure = 3
!* Total number of slip systems per lattice structure
!* (has to be changed according the definition of slip systems)
integer(pInt), dimension(lattice_MaxLatticeStructure), parameter :: lattice_MaxNslipOfStructure = &
reshape((/12,48,24/),(/lattice_MaxLatticeStructure/))
!* Total number of twin systems per lattice structure
!* (has to be changed according the definition of twin systems)
integer(pInt), dimension(lattice_MaxLatticeStructure), parameter :: lattice_MaxNtwinOfStructure = &
reshape((/12,0,24/),(/lattice_MaxLatticeStructure/))
!* Maximum number of slip systems over lattice structures
integer(pInt), parameter :: lattice_MaxMaxNslipOfStructure = 48
!* Maximum number of twin systems over lattice structures,  changed form 12 to 24 (yj.ro)
integer(pInt), parameter :: lattice_MaxMaxNtwinOfStructure = 24
!* Slip direction, slip normales and Schmid matrices
real(pReal), dimension(3,3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_Sslip
real(pReal), dimension(6,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_Sslip_v
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_sn
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_sd
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_st

!* HCP - slip direction, slip normal (4 indices): Prof. Tom Bieler, Leeyun, YJRO 
real(pReal), dimension(4,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: Hlattice_sn
real(pReal), dimension(4,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: Hlattice_sd
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: H_lattice_sn
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: H_lattice_sd

!* twin direction, twin normales, Schmid matrices and transformation matrices
real(pReal), dimension(3,3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_Stwin
real(pReal), dimension(6,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_Stwin_v
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_tn
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_td
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_tt
real(pReal), dimension(3,3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_Qtwin

!* HCP - twin direction, twin normales for 4 indices: Prof. Tom Bieler, Leeyun, YJR 
real(pReal), dimension(4,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: Hlattice_tn
real(pReal), dimension(4,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: Hlattice_td
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: H_lattice_tn
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: H_lattice_td


real(pReal), dimension(lattice_MaxLatticeStructure), parameter :: lattice_TwinShear = &
reshape((/0.7071067812,0.7071067812,0.7071067812/),(/lattice_MaxLatticeStructure/)) ! Depends surely on c/a ratio for HCP


!* Slip_slip interaction matrices
integer(pInt), dimension(lattice_MaxMaxNslipOfStructure,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: &
lattice_SlipIntType
!* Slip_twin interaction matrices
integer(pInt), dimension(lattice_MaxMaxNslipOfStructure,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: &
lattice_SlipTwinIntType
!* Twin-twin interaction matrices
integer(pInt), dimension(lattice_MaxMaxNtwinOfStructure,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: &
lattice_TwinIntType

!*** Slip systems for FCC structures (1) ***
!* System {111}<110>  Sort according Eisenlohr&Hantcherli
data lattice_sd(:, 1,1)/ 0, 1,-1/ ; data lattice_sn(:, 1,1)/ 1, 1, 1/
data lattice_sd(:, 2,1)/-1, 0, 1/ ; data lattice_sn(:, 2,1)/ 1, 1, 1/
data lattice_sd(:, 3,1)/ 1,-1, 0/ ; data lattice_sn(:, 3,1)/ 1, 1, 1/
data lattice_sd(:, 4,1)/ 0,-1,-1/ ; data lattice_sn(:, 4,1)/-1,-1, 1/
data lattice_sd(:, 5,1)/ 1, 0, 1/ ; data lattice_sn(:, 5,1)/-1,-1, 1/
data lattice_sd(:, 6,1)/-1, 1, 0/ ; data lattice_sn(:, 6,1)/-1,-1, 1/
data lattice_sd(:, 7,1)/ 0,-1, 1/ ; data lattice_sn(:, 7,1)/ 1,-1,-1/
data lattice_sd(:, 8,1)/-1, 0,-1/ ; data lattice_sn(:, 8,1)/ 1,-1,-1/
data lattice_sd(:, 9,1)/ 1, 1, 0/ ; data lattice_sn(:, 9,1)/ 1,-1,-1/
data lattice_sd(:,10,1)/ 0, 1, 1/ ; data lattice_sn(:,10,1)/-1, 1,-1/
data lattice_sd(:,11,1)/ 1, 0,-1/ ; data lattice_sn(:,11,1)/-1, 1,-1/
data lattice_sd(:,12,1)/-1,-1, 0/ ; data lattice_sn(:,12,1)/-1, 1,-1/

!*** Twin systems for FCC structures (1) ***
!* System {111}<112>  Sort according Eisenlohr&Hantcherli
data lattice_td(:, 1,1)/-2, 1, 1/ ; data lattice_tn(:, 1,1)/ 1, 1, 1/
data lattice_td(:, 2,1)/ 1,-2, 1/ ; data lattice_tn(:, 2,1)/ 1, 1, 1/
data lattice_td(:, 3,1)/ 1, 1,-2/ ; data lattice_tn(:, 3,1)/ 1, 1, 1/
data lattice_td(:, 4,1)/ 2,-1, 1/ ; data lattice_tn(:, 4,1)/-1,-1, 1/
data lattice_td(:, 5,1)/-1, 2, 1/ ; data lattice_tn(:, 5,1)/-1,-1, 1/
data lattice_td(:, 6,1)/-1,-1,-2/ ; data lattice_tn(:, 6,1)/-1,-1, 1/
data lattice_td(:, 7,1)/-2,-1,-1/ ; data lattice_tn(:, 7,1)/ 1,-1,-1/
data lattice_td(:, 8,1)/ 1, 2,-1/ ; data lattice_tn(:, 8,1)/ 1,-1,-1/
data lattice_td(:, 9,1)/ 1,-1, 2/ ; data lattice_tn(:, 9,1)/ 1,-1,-1/
data lattice_td(:,10,1)/ 2, 1,-1/ ; data lattice_tn(:,10,1)/-1, 1,-1/
data lattice_td(:,11,1)/-1,-2,-1/ ; data lattice_tn(:,11,1)/-1, 1,-1/
data lattice_td(:,12,1)/-1, 1, 2/ ; data lattice_tn(:,12,1)/-1, 1,-1/

!*** Slip-Slip interactions for FCC structures (1) ***
data lattice_SlipIntType( 1,1:lattice_MaxNslipOfStructure(1),1)/1,2,2,4,6,5,3,5,5,4,5,6/
data lattice_SlipIntType( 2,1:lattice_MaxNslipOfStructure(1),1)/2,1,2,6,4,5,5,4,6,5,3,5/
data lattice_SlipIntType( 3,1:lattice_MaxNslipOfStructure(1),1)/2,2,1,5,5,3,5,6,4,6,5,4/
data lattice_SlipIntType( 4,1:lattice_MaxNslipOfStructure(1),1)/4,6,5,1,2,2,4,5,6,3,5,5/
data lattice_SlipIntType( 5,1:lattice_MaxNslipOfStructure(1),1)/6,4,5,2,1,2,5,3,5,5,4,6/
data lattice_SlipIntType( 6,1:lattice_MaxNslipOfStructure(1),1)/5,5,3,2,2,1,6,5,4,5,6,4/
data lattice_SlipIntType( 7,1:lattice_MaxNslipOfStructure(1),1)/3,5,5,4,5,6,1,2,2,4,6,5/
data lattice_SlipIntType( 8,1:lattice_MaxNslipOfStructure(1),1)/5,4,6,5,3,5,2,1,2,6,4,5/
data lattice_SlipIntType( 9,1:lattice_MaxNslipOfStructure(1),1)/5,6,4,6,5,4,2,2,1,5,5,3/
data lattice_SlipIntType(10,1:lattice_MaxNslipOfStructure(1),1)/4,5,6,3,5,5,4,6,5,1,2,2/
data lattice_SlipIntType(11,1:lattice_MaxNslipOfStructure(1),1)/5,3,5,5,4,6,6,4,5,2,1,2/
data lattice_SlipIntType(12,1:lattice_MaxNslipOfStructure(1),1)/6,5,4,5,6,4,5,5,3,2,2,1/

!*** Slip-Twin interactions for FCC structures (1) ***
data lattice_SlipTwinIntType( 1,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,0,0,0,1,1,1/
data lattice_SlipTwinIntType( 2,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,1,1,1,0,0,0/
data lattice_SlipTwinIntType( 3,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,0,0,0,1,1,1,1,1,1/
data lattice_SlipTwinIntType( 4,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,1,1,1,0,0,0/
data lattice_SlipTwinIntType( 5,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,0,0,0,1,1,1/
data lattice_SlipTwinIntType( 6,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,0,0,0,1,1,1,1,1,1/
data lattice_SlipTwinIntType( 7,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,0,0,0,1,1,1/
data lattice_SlipTwinIntType( 8,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,0,0,0,1,1,1/
data lattice_SlipTwinIntType( 9,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,0,0,0,0,0,0/
data lattice_SlipTwinIntType(10,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,1,1,1,0,0,0/
data lattice_SlipTwinIntType(11,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,1,1,1,0,0,0/
data lattice_SlipTwinIntType(12,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,0,0,0,0,0,0/

!*** Twin-Twin interactions for FCC structures (1) ***
data lattice_TwinIntType( 1,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,1,1,1,1,1,1/
data lattice_TwinIntType( 2,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,1,1,1,1,1,1/
data lattice_TwinIntType( 3,1:lattice_MaxNtwinOfStructure(1),1)/0,0,0,1,1,1,1,1,1,1,1,1/
data lattice_TwinIntType( 4,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,1,1,1,1,1,1/
data lattice_TwinIntType( 5,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,1,1,1,1,1,1/
data lattice_TwinIntType( 6,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,0,0,0,1,1,1,1,1,1/
data lattice_TwinIntType( 7,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,0,0,0,1,1,1/
data lattice_TwinIntType( 8,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,0,0,0,1,1,1/
data lattice_TwinIntType( 9,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,0,0,0,1,1,1/
data lattice_TwinIntType(10,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,1,1,1,0,0,0/
data lattice_TwinIntType(11,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,1,1,1,0,0,0/
data lattice_TwinIntType(12,1:lattice_MaxNtwinOfStructure(1),1)/1,1,1,1,1,1,1,1,1,0,0,0/

!*** Slip systems for BCC structures (2) ***
!* System {110}<111>
!* Sort?
data lattice_sd(:, 1,2)/ 1,-1, 1/ ; data lattice_sn(:, 1,2)/ 0, 1, 1/
data lattice_sd(:, 2,2)/-1,-1, 1/ ; data lattice_sn(:, 2,2)/ 0, 1, 1/
data lattice_sd(:, 3,2)/ 1, 1, 1/ ; data lattice_sn(:, 3,2)/ 0,-1, 1/
data lattice_sd(:, 4,2)/-1, 1, 1/ ; data lattice_sn(:, 4,2)/ 0,-1, 1/
data lattice_sd(:, 5,2)/-1, 1, 1/ ; data lattice_sn(:, 5,2)/ 1, 0, 1/
data lattice_sd(:, 6,2)/-1,-1, 1/ ; data lattice_sn(:, 6,2)/ 1, 0, 1/
data lattice_sd(:, 7,2)/ 1, 1, 1/ ; data lattice_sn(:, 7,2)/-1, 0, 1/
data lattice_sd(:, 8,2)/ 1,-1, 1/ ; data lattice_sn(:, 8,2)/-1, 0, 1/
data lattice_sd(:, 9,2)/-1, 1, 1/ ; data lattice_sn(:, 9,2)/ 1, 1, 0/
data lattice_sd(:,10,2)/-1, 1,-1/ ; data lattice_sn(:,10,2)/ 1, 1, 0/
data lattice_sd(:,11,2)/ 1, 1, 1/ ; data lattice_sn(:,11,2)/-1, 1, 0/
data lattice_sd(:,12,2)/ 1, 1,-1/ ; data lattice_sn(:,12,2)/-1, 1, 0/
!* System {112}<111>
!* Sort?
data lattice_sd(:,13,2)/-1, 1, 1/ ; data lattice_sn(:,13,2)/ 2, 1, 1/
data lattice_sd(:,14,2)/ 1, 1, 1/ ; data lattice_sn(:,14,2)/-2, 1, 1/
data lattice_sd(:,15,2)/ 1, 1,-1/ ; data lattice_sn(:,15,2)/ 2,-1, 1/
data lattice_sd(:,16,2)/ 1,-1, 1/ ; data lattice_sn(:,16,2)/ 2, 1,-1/
data lattice_sd(:,17,2)/ 1,-1, 1/ ; data lattice_sn(:,17,2)/ 1, 2, 1/
data lattice_sd(:,18,2)/ 1, 1,-1/ ; data lattice_sn(:,18,2)/-1, 2, 1/
data lattice_sd(:,19,2)/ 1, 1, 1/ ; data lattice_sn(:,19,2)/ 1,-2, 1/
data lattice_sd(:,20,2)/-1, 1, 1/ ; data lattice_sn(:,20,2)/ 1, 2,-1/
data lattice_sd(:,21,2)/ 1, 1,-1/ ; data lattice_sn(:,21,2)/ 1, 1, 2/
data lattice_sd(:,22,2)/ 1,-1, 1/ ; data lattice_sn(:,22,2)/-1, 1, 2/
data lattice_sd(:,23,2)/-1, 1, 1/ ; data lattice_sn(:,23,2)/ 1,-1, 2/
data lattice_sd(:,24,2)/ 1, 1, 1/ ; data lattice_sn(:,24,2)/ 1, 1,-2/
!* System {123}<111>
!* Sort?
data lattice_sd(:,25,2)/ 1, 1,-1/ ; data lattice_sn(:,25,2)/ 1, 2, 3/
data lattice_sd(:,26,2)/ 1,-1, 1/ ; data lattice_sn(:,26,2)/-1, 2, 3/
data lattice_sd(:,27,2)/-1, 1, 1/ ; data lattice_sn(:,27,2)/ 1,-2, 3/
data lattice_sd(:,28,2)/ 1, 1, 1/ ; data lattice_sn(:,28,2)/ 1, 2,-3/
data lattice_sd(:,29,2)/ 1,-1, 1/ ; data lattice_sn(:,29,2)/ 1, 3, 2/
data lattice_sd(:,30,2)/ 1, 1,-1/ ; data lattice_sn(:,30,2)/-1, 3, 2/
data lattice_sd(:,31,2)/ 1, 1, 1/ ; data lattice_sn(:,31,2)/ 1,-3, 2/
data lattice_sd(:,32,2)/-1, 1, 1/ ; data lattice_sn(:,32,2)/ 1, 3,-2/
data lattice_sd(:,33,2)/ 1, 1,-1/ ; data lattice_sn(:,33,2)/ 2, 1, 3/
data lattice_sd(:,34,2)/ 1,-1, 1/ ; data lattice_sn(:,34,2)/-2, 1, 3/
data lattice_sd(:,35,2)/-1, 1, 1/ ; data lattice_sn(:,35,2)/ 2,-1, 3/
data lattice_sd(:,36,2)/ 1, 1, 1/ ; data lattice_sn(:,36,2)/ 2, 1,-3/
data lattice_sd(:,37,2)/ 1,-1, 1/ ; data lattice_sn(:,37,2)/ 2, 3, 1/
data lattice_sd(:,38,2)/ 1, 1,-1/ ; data lattice_sn(:,38,2)/-2, 3, 1/
data lattice_sd(:,39,2)/ 1, 1, 1/ ; data lattice_sn(:,39,2)/ 2,-3, 1/
data lattice_sd(:,40,2)/-1, 1, 1/ ; data lattice_sn(:,40,2)/ 2, 3,-1/
data lattice_sd(:,41,2)/-1, 1, 1/ ; data lattice_sn(:,41,2)/ 3, 1, 2/
data lattice_sd(:,42,2)/ 1, 1, 1/ ; data lattice_sn(:,42,2)/-3, 1, 2/
data lattice_sd(:,43,2)/ 1, 1,-1/ ; data lattice_sn(:,43,2)/ 3,-1, 2/
data lattice_sd(:,44,2)/ 1,-1, 1/ ; data lattice_sn(:,44,2)/ 3, 1,-2/
data lattice_sd(:,45,2)/-1, 1, 1/ ; data lattice_sn(:,45,2)/ 3, 2, 1/
data lattice_sd(:,46,2)/ 1, 1, 1/ ; data lattice_sn(:,46,2)/-3, 2, 1/
data lattice_sd(:,47,2)/ 1, 1,-1/ ; data lattice_sn(:,47,2)/ 3,-2, 1/
data lattice_sd(:,48,2)/ 1,-1, 1/ ; data lattice_sn(:,48,2)/ 3, 2,-1/

!*** Twin systems for BCC structures (2) ***
!* System {112}<111>
!* Sort?
!* MISSING: not implemented yet

!*** Slip-Slip interactions for BCC structures (2) ***
data lattice_SlipIntType( 1,:,2)/1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 2,:,2)/2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 3,:,2)/2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 4,:,2)/2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 5,:,2)/2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 6,:,2)/2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 7,:,2)/2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 8,:,2)/2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 9,:,2)/2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(10,:,2)/2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(11,:,2)/2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(12,:,2)/2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(13,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(14,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(15,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(16,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(17,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(18,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(19,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(20,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(21,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(22,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(23,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(24,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(25,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(26,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(27,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(28,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(29,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(30,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(31,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(32,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(33,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(34,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(35,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(36,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(37,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(38,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(39,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(40,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(41,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2/
data lattice_SlipIntType(42,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2/
data lattice_SlipIntType(43,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2/
data lattice_SlipIntType(44,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2/
data lattice_SlipIntType(45,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2/
data lattice_SlipIntType(46,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2/
data lattice_SlipIntType(47,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2/
data lattice_SlipIntType(48,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1/

!*** Slip-twin interactions for BCC structures (2) ***
! MISSING: not implemented yet

!*** Twin-twin interactions for BCC structures (2) ***
! MISSING: not implemented yet

!*** Slip systems for HCP structures (3) ***
!* Basal systems {0001}<1120> (independent of c/a-ratio)
!* 1- [2  -1  -1  0](0 0 0 1)
!* 2- [-1  2  -1  0](0 0 0 1)
!* 3- [-1 -1   2  0](0 0 0 1)
!* Automatical transformation from Bravais (4 axes coorinate system) to Miller (in ortho-hexagonal coordinate system)
!* not done for the moment
!* Sort? Changed order of slip system and sign of Burges vector (Tom Bieler, yj.ro)
data Hlattice_sd(:, 1,3)/ 2, -1, -1,  0/ ; data Hlattice_sn(:, 1,3)/ 0,  0,  0,  1/
data Hlattice_sd(:, 2,3)/-1,  2, -1,  0/ ; data Hlattice_sn(:, 2,3)/ 0,  0,  0,  1/
data Hlattice_sd(:, 3,3)/-1, -1,  2,  0/ ; data Hlattice_sn(:, 3,3)/ 0,  0,  0,  1/
!* 1st type prismatic systems {1010}<1120>  (independent of c/a-ratio)
!* 4- [ 2 -1 -1  0]( 0  1 -1  0)
!* 5- [-1  2 -1  0]( 1  0 -1  0)
!* 6- [-1 -1  2  0](-1  1  0  0)
!* Sort? Changed order of slip system and sign of Burges vector (yj.ro)
data Hlattice_sd(:, 4,3)/ 2, -1, -1,  0/ ; data Hlattice_sn(:, 4,3)/ 0,  1, -1,  0/
data Hlattice_sd(:, 5,3)/-1,  2, -1,  0/ ; data Hlattice_sn(:, 5,3)/ 1,  0, -1,  0/
data Hlattice_sd(:, 6,3)/-1, -1,  2,  0/ ; data Hlattice_sn(:, 6,3)/-1,  1,  0,  0/
!* 1st type 1st order pyramidal systems {1011}<1120>
!* plane normales depend on the c/a-ratio
!*	7- [ 2 -1 -1  0]( 0  1 -1  1)
!*  8- [-1  2 -1  0]( 1  0 -1  1)
!*  9- [-1 -1  2  0](-1  1  0  1)
!* 10- [ 2 -1 -1  0]( 0 -1  1  1)
!* 11- [-1  2 -1  0](-1  0  1  1)
!* 12- [-1 -1  2  0]( 1 -1  0  1)
!* Sort? Changed order of slip system and sign of Burges vector (Tom Bieler, yj.ro)
data Hlattice_sd(:, 7,3)/ 2, -1, -1,  0/ ; data Hlattice_sn(:, 7,3)/ 0,  1, -1,  1/
data Hlattice_sd(:, 8,3)/-1,  2, -1,  0/ ; data Hlattice_sn(:, 8,3)/ 1,  0, -1,  1/
data Hlattice_sd(:, 9,3)/-1, -1,  2,  0/ ; data Hlattice_sn(:, 9,3)/-1,  1,  0,  1/
data Hlattice_sd(:,10,3)/ 2, -1, -1,  0/ ; data Hlattice_sn(:,10,3)/ 0, -1,  1,  1/
data Hlattice_sd(:,11,3)/-1,  2, -1,  0/ ; data Hlattice_sn(:,11,3)/-1,  0,  1,  1/
data Hlattice_sd(:,12,3)/-1, -1,  2,  0/ ; data Hlattice_sn(:,12,3)/ 1, -1,  0,  1/
!* pyramidal system: c+a slip {1011}<2113>
!* plane normales depend on the c/a-ratio
!* added by Tom Bieler, yj.ro
!* 13- [ 2 -1 -1 -3]( 1  0 -1  1)
!* 14- [ 1  1 -2 -3]( 1  0 -1  1)
!* 15- [ 1  1 -2 -3]( 0  1 -1  1)
!* 16- [-1  2 -1 -3]( 0  1 -1  1)
!* 17- [-1  2 -1 -3](-1  1  0  1)
!* 18- [-2  1  1 -3](-1  1  0  1)
!* 19- [-2  1  1 -3](-1  0  1  1)
!* 20- [-1 -1  2 -3](-1  0  1  1)
!* 21- [-1 -1  2 -3]( 0 -1  1  1)
!* 22- [ 1 -2  1 -3]( 0 -1  1  1)
!* 23- [ 1 -2  1 -3]( 1 -1  0  1)
!* 24- [ 2 -1 -1 -3]( 1 -1  0  1)
data Hlattice_sd(:,13,3)/ 2, -1, -1, -3/ ; data Hlattice_sn(:,13,3)/ 1,  0, -1,  1/
data Hlattice_sd(:,14,3)/ 1,  1, -2, -3/ ; data Hlattice_sn(:,14,3)/ 1,  0, -1,  1/
data Hlattice_sd(:,15,3)/ 1,  1, -2, -3/ ; data Hlattice_sn(:,15,3)/ 0,  1, -1,  1/
data Hlattice_sd(:,16,3)/-1,  2, -1, -3/ ; data Hlattice_sn(:,16,3)/ 0,  1, -1,  1/
data Hlattice_sd(:,17,3)/-1,  2, -1, -3/ ; data Hlattice_sn(:,17,3)/-1,  1,  0,  1/
data Hlattice_sd(:,18,3)/-2,  1,  1, -3/ ; data Hlattice_sn(:,18,3)/-1,  1,  0,  1/
data Hlattice_sd(:,19,3)/-2,  1,  1, -3/ ; data Hlattice_sn(:,19,3)/-1,  0,  1,  1/
data Hlattice_sd(:,20,3)/-1, -1,  2, -3/ ; data Hlattice_sn(:,20,3)/-1,  0,  1,  1/
data Hlattice_sd(:,21,3)/-1, -1,  2, -3/ ; data Hlattice_sn(:,21,3)/ 0, -1,  1,  1/
data Hlattice_sd(:,22,3)/ 1, -2,  1, -3/ ; data Hlattice_sn(:,22,3)/ 0, -1,  1,  1/
data Hlattice_sd(:,23,3)/ 1, -2,  1, -3/ ; data Hlattice_sn(:,23,3)/ 1, -1,  0,  1/
data Hlattice_sd(:,24,3)/ 2, -1, -1, -3/ ; data Hlattice_sn(:,24,3)/ 1, -1,  0,  1/

!*** Twin systems for HCP structures (3) ***
!* Sort? Numbering of twin system follows Prof. Tom Bieler's scheme (to be consistent with his work); but numbering in data was restarted from 1 &
!*(to be consistent with this code structure).
!* MISSING: not implemented yet
!* added by Tom Bieler, yj.ro

!* (1012)<1011> Twin: shear 0.169 -1.26 compression
!* 25- [-1  0  1  1]( 1  0 -1  2)
!* 26- [ 0 -1  1  1]( 0  1 -1  2)
!* 27- [ 1 -1  0  1](-1  1  0  2)
!* 28- [ 1  0 -1  1](-1  0  1  2)
!* 29- [ 0  1 -1  1]( 0 -1  1  2)
!* 30- [-1  1  0  1]( 1 -1  0  2)
data Hlattice_td(:, 1,3)/-1,  0,  1,  1/ ; data Hlattice_tn(:, 1,3)/ 1,  0, -1,  2/
data Hlattice_td(:, 2,3)/ 0, -1,  1,  1/ ; data Hlattice_tn(:, 2,3)/ 0,  1, -1,  2/
data Hlattice_td(:, 3,3)/ 1, -1,  0,  1/ ; data Hlattice_tn(:, 3,3)/-1,  1,  0,  2/
data Hlattice_td(:, 4,3)/ 1,  0, -1,  1/ ; data Hlattice_tn(:, 4,3)/-1,  0,  1,  2/
data Hlattice_td(:, 5,3)/ 0,  1, -1,  1/ ; data Hlattice_tn(:, 5,3)/ 0, -1,  1,  2/
data Hlattice_td(:, 6,3)/-1,  1,  0,  1/ ; data Hlattice_tn(:, 6,3)/ 1, -1,  0,  2/


!*(2112)<211-2> Twin: shear 0.224 1.19 tension
!* 31- [ 2 -1 -1 -3]( 2 -1 -1  2)
!* 32- [ 1  1 -2 -3]( 1  1 -2  2)
!* 33- [-1  2 -1 -3](-1  2 -1  2)
!* 34- [-2  1  1 -3](-2  1  1  2)
!* 35- [-1 -1  2 -3](-1 -1  2  2)
!* 36- [ 1 -2  1 -3]( 1 -2  1  2)
data Hlattice_td(:, 7,3)/ 2, -1, -1, -3/ ; data Hlattice_tn(:, 7,3)/ 2, -1, -1,  2/
data Hlattice_td(:, 8,3)/ 1,  1, -2, -3/ ; data Hlattice_tn(:, 8,3)/ 1,  1, -2,  2/
data Hlattice_td(:, 9,3)/-1,  2, -1, -3/ ; data Hlattice_tn(:, 9,3)/-1,  2, -1,  2/
data Hlattice_td(:,10,3)/-2,  1,  1, -3/ ; data Hlattice_tn(:,10,3)/-2,  1,  1,  2/
data Hlattice_td(:,11,3)/-1, -1,  2, -3/ ; data Hlattice_tn(:,11,3)/-1, -1,  2,  2/
data Hlattice_td(:,12,3)/1,  -2,  1, -3/ ; data Hlattice_tn(:,12,3)/ 1, -2,  1,  2/


!* (2111)<211-6> Twin: shear 0.628 -0.39 compressio
!* 37- [-2  1  1  6]( 2 -1 -1  1)
!* 38- [-1 -1  2  6]( 1  1 -2  1)
!* 39- [ 1 -2  1  6](-1  2 -1  1)
!* 40- [ 2 -1 -1  6](-2  1  1  1)
!* 41- [ 1  1 -2  6](-1 -1  2  1)
!* 42- [-1  2 -1  6]( 1 -2  1  1)
data Hlattice_td(:,13,3)/-2,  1,  1,  6/ ; data Hlattice_tn(:,13,3)/ 2, -1, -1,  1/
data Hlattice_td(:,14,3)/-1, -1,  2,  6/ ; data Hlattice_tn(:,14,3)/ 1,  1, -2,  1/
data Hlattice_td(:,15,3)/ 1, -2,  1,  6/ ; data Hlattice_tn(:,15,3)/-1,  2, -1,  1/
data Hlattice_td(:,16,3)/ 2, -1, -1,  6/ ; data Hlattice_tn(:,16,3)/-2,  1,  1,  1/
data Hlattice_td(:,17,3)/ 1,  1, -2,  6/ ; data Hlattice_tn(:,17,3)/-1, -1,  2,  1/
data Hlattice_td(:,18,3)/-1,  2, -1,  6/ ; data Hlattice_tn(:,18,3)/ 1, -2,  1,  1/


!* (1011)<101-2> Twin: shear 0.103 1.09 tension
!* 43- [ 1  0 -1 -2]( 1  0 -1  1)
!* 44- [-1  0  1 -2](-1  0  1  1)
!* 45- [ 0  1 -1 -2]( 0  1 -1  1)
!* 46- [ 0 -1  2 -2]( 0 -1  1  1)
!* 47- [ 1 -1  0 -2]( 1 -1  0  1)
!* 48- [-1  1  0 -2](-1  1  0  1)
data Hlattice_td(:,19,3)/ 1,  0, -1, -2/ ; data Hlattice_tn(:,19,3)/ 1,  0, -1,  1/
data Hlattice_td(:,20,3)/-1,  0,  1, -2/ ; data Hlattice_tn(:,20,3)/-1,  0,  1,  1/
data Hlattice_td(:,21,3)/ 0,  1, -1, -2/ ; data Hlattice_tn(:,21,3)/ 0,  1, -1,  1/
data Hlattice_td(:,22,3)/ 0, -1,  1, -2/ ; data Hlattice_tn(:,22,3)/ 0, -1,  1,  1/
data Hlattice_td(:,23,3)/ 1, -1,  0, -2/ ; data Hlattice_tn(:,23,3)/ 1, -1,  0,  1/
data Hlattice_td(:,24,3)/-1,  1,  0, -2/ ; data Hlattice_tn(:,24,3)/-1,  1,  0,  1/


!*** Slip-Slip interactions for HCP structures (3) ***
data lattice_SlipIntType( 1,1:lattice_MaxNslipOfStructure(3),3)/1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 2,1:lattice_MaxNslipOfStructure(3),3)/1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 3,1:lattice_MaxNslipOfStructure(3),3)/1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 4,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 5,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 6,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 7,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 8,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 9,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(10,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(11,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(12,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(13,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(14,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(15,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(16,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2/
data lattice_SlipIntType(17,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2/
data lattice_SlipIntType(18,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2/
data lattice_SlipIntType(19,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2/
data lattice_SlipIntType(20,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2/
data lattice_SlipIntType(21,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2/
data lattice_SlipIntType(22,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2/
data lattice_SlipIntType(23,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2/
data lattice_SlipIntType(24,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1/

!*** slip-twin interactions for HCP structures (3) ***
data lattice_SlipTwinIntType( 1,1:lattice_MaxNtwinOfStructure(3),3)/1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 2,1:lattice_MaxNtwinOfStructure(3),3)/2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 3,1:lattice_MaxNtwinOfStructure(3),3)/2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 4,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 5,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 6,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 7,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 8,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType( 9,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(10,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(11,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(12,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(13,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(14,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(15,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(16,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(17,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2/
data lattice_SlipTwinIntType(18,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2/
data lattice_SlipTwinIntType(19,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2/
data lattice_SlipTwinIntType(20,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2/
data lattice_SlipTwinIntType(21,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2/
data lattice_SlipTwinIntType(22,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2/
data lattice_SlipTwinIntType(23,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2/
data lattice_SlipTwinIntType(24,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1/

!*** Twin-twin interactions for HCP structures (3) ***
data lattice_TwinIntType( 1,1:lattice_MaxNtwinOfStructure(3),3)/1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 2,1:lattice_MaxNtwinOfStructure(3),3)/2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 3,1:lattice_MaxNtwinOfStructure(3),3)/2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 4,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 5,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 6,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 7,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 8,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType( 9,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(10,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(11,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(12,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(13,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(14,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(15,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(16,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2/
data lattice_TwinIntType(17,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2/
data lattice_TwinIntType(18,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2/
data lattice_TwinIntType(19,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2/
data lattice_TwinIntType(20,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2/
data lattice_TwinIntType(21,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2/
data lattice_TwinIntType(22,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2/
data lattice_TwinIntType(23,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2/
data lattice_TwinIntType(24,1:lattice_MaxNtwinOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1/



CONTAINS
!****************************************
!* - lattice_Init
!* - lattice_SchmidMatrices
!****************************************


subroutine lattice_init()
!**************************************
!*      Module initialization         *
!**************************************
call lattice_SchmidMatrices()
end subroutine


subroutine lattice_SchmidMatrices()
!**************************************
!*   Calculation of Schmid matrices   *
!**************************************
use prec, only: pReal,pInt
use math, only: math_I3,nrmMandel,mapMandel
implicit none

!* Definition of variables
integer(pInt) i,j,k,l
real(pReal) norm_d,norm_t,norm_n
real(pReal) norm_sn, norm_sd, norm_tn, norm_td, ratio

!*** Only HCP crystal: converting from 4 axes coordinate system (a1=a2=a3=c) to ortho-hexgonal system (a, b, c)
!* Plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a)): this has been changed to unit vector afterward.
!* Direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]: this has been changed to unit vector afterward.
!* Equations provided by Prof. Tom Bieler
!* need to input "c/a rati"o in somewhere in mattax.mpie file(I am not sure where to insert this value for now????).

ratio = 1.56

do i = 1,lattice_MaxNslipOfStructure(3)

!* slip system conversion
    H_lattice_sn(1,i,3) = Hlattice_sn(1,i,3)
    H_lattice_sn(2,i,3) = (Hlattice_sn(1,i,3)+ 2.0*Hlattice_sn(2,i,3))/sqrt(3.0)
    H_lattice_sn(3,i,3) = Hlattice_sn(4,i,3)/ratio
    
    norm_sn = dsqrt(H_lattice_sn(1,i,3)**2 + H_lattice_sn(2,i,3)**2 + H_lattice_sn(3,i,3)**2)
    
    lattice_sn(1,i,3) = H_lattice_sn(1,i,3)/norm_sn
    lattice_sn(2,i,3) = H_lattice_sn(2,i,3)/norm_sn
    lattice_sn(3,i,3) = H_lattice_sn(3,i,3)/norm_sn
    
    H_lattice_sd(1,i,3) = 1.5*Hlattice_sd(1,i,3)
    H_lattice_sd(2,i,3) = (Hlattice_sd(1,i,3) + 2.0*Hlattice_sd(2,i,3))*(sqrt(3.0)/2.0)
    H_lattice_sd(3,i,3) = Hlattice_sd(4,i,3)*ratio

    norm_sd = dsqrt(H_lattice_sd(1,i,3)**2 + H_lattice_sd(2,i,3)**2 + H_lattice_sd(3,i,3)**2)
    
    lattice_sd(1,i,3) = H_lattice_sd(1,i,3)/norm_sd
    lattice_sd(2,i,3) = H_lattice_sd(2,i,3)/norm_sd
    lattice_sd(3,i,3) = H_lattice_sd(3,i,3)/norm_sd

!* twin system conversion
    H_lattice_tn(1,i,3) = Hlattice_tn(1,i,3)
    H_lattice_tn(2,i,3) = (Hlattice_tn(1,i,3)+ 2.0*Hlattice_tn(2,i,3))/sqrt(3.0)
    H_lattice_tn(3,i,3) = Hlattice_tn(4,i,3)/ratio
    
    norm_tn = dsqrt(H_lattice_tn(1,i,3)**2 + H_lattice_tn(2,i,3)**2 + H_lattice_tn(3,i,3)**2)
    
    lattice_tn(1,i,3) = H_lattice_tn(1,i,3)/norm_tn
    lattice_tn(2,i,3) = H_lattice_tn(2,i,3)/norm_tn
    lattice_tn(3,i,3) = H_lattice_tn(3,i,3)/norm_tn
    
    H_lattice_td(1,i,3) = 1.5*Hlattice_td(1,i,3)
    H_lattice_td(2,i,3) = (Hlattice_td(1,i,3)+ 2.0*Hlattice_td(2,i,3))*(sqrt(3.0)/2.0)
    H_lattice_td(3,i,3) = Hlattice_td(4,i,3)*ratio
    
    norm_td = dsqrt(H_lattice_td(1,i,3)**2 + H_lattice_td(2,i,3)**2 + H_lattice_td(3,i,3)**2)
    
    lattice_td(1,i,3) = H_lattice_td(1,i,3)/norm_td
    lattice_td(2,i,3) = H_lattice_td(2,i,3)/norm_td
    lattice_td(3,i,3) = H_lattice_td(3,i,3)/norm_td
    
enddo

!* Iteration over the lattice structures
do l=1,lattice_MaxLatticeStructure
!* Iteration over the slip systems
   do k=1,lattice_MaxNslipOfStructure(l)
!* Definition of transverse direction st for the frame (sd,st,sn)
      lattice_st(1,k,l)=lattice_sn(2,k,l)*lattice_sd(3,k,l)-lattice_sn(3,k,l)*lattice_sd(2,k,l)
	  lattice_st(2,k,l)=lattice_sn(3,k,l)*lattice_sd(1,k,l)-lattice_sn(1,k,l)*lattice_sd(3,k,l)
	  lattice_st(3,k,l)=lattice_sn(1,k,l)*lattice_sd(2,k,l)-lattice_sn(2,k,l)*lattice_sd(1,k,l)
	  norm_d=dsqrt(lattice_sd(1,k,l)**2+lattice_sd(2,k,l)**2+lattice_sd(3,k,l)**2)
      norm_t=dsqrt(lattice_st(1,k,l)**2+lattice_st(2,k,l)**2+lattice_st(3,k,l)**2)
      norm_n=dsqrt(lattice_sn(1,k,l)**2+lattice_sn(2,k,l)**2+lattice_sn(3,k,l)**2)
      lattice_sd(:,k,l)=lattice_sd(:,k,l)/norm_d
	  lattice_st(:,k,l)=lattice_st(:,k,l)/norm_t
	  lattice_sn(:,k,l)=lattice_sn(:,k,l)/norm_n
!* Defintion of Schmid matrix
      forall (i=1:3,j=1:3) lattice_Sslip(i,j,k,l)=lattice_sd(i,k,l)*lattice_sn(j,k,l)
!* Vectorization of normalized Schmid matrix
      forall (i=1:6) lattice_Sslip_v(i,k,l) = nrmMandel(i)/2.0_pReal * &
                  (lattice_Sslip(mapMandel(1,i),mapMandel(2,i),k,l)+lattice_Sslip(mapMandel(2,i),mapMandel(1,i),k,l))
   enddo

!* Iteration over the twin systems
   do k=1,lattice_MaxNtwinOfStructure(l)
!* Definition of transverse direction tt for the frame (td,tt,tn)
      lattice_tt(1,k,l)=lattice_tn(2,k,l)*lattice_td(3,k,l)-lattice_tn(3,k,l)*lattice_td(2,k,l)
	  lattice_tt(2,k,l)=lattice_tn(3,k,l)*lattice_td(1,k,l)-lattice_tn(1,k,l)*lattice_td(3,k,l)
	  lattice_tt(3,k,l)=lattice_tn(1,k,l)*lattice_td(2,k,l)-lattice_tn(2,k,l)*lattice_td(1,k,l)
	  norm_d=dsqrt(lattice_td(1,k,l)**2+lattice_td(2,k,l)**2+lattice_td(3,k,l)**2)
      norm_t=dsqrt(lattice_tt(1,k,l)**2+lattice_tt(2,k,l)**2+lattice_tt(3,k,l)**2)
      norm_n=dsqrt(lattice_tn(1,k,l)**2+lattice_tn(2,k,l)**2+lattice_tn(3,k,l)**2)
      lattice_td(:,k,l)=lattice_td(:,k,l)/norm_d
	  lattice_tt(:,k,l)=lattice_tt(:,k,l)/norm_t
	  lattice_tn(:,k,l)=lattice_tn(:,k,l)/norm_n
!* Defintion of Schmid matrix and transformation matrices
      lattice_Qtwin(:,:,k,l)=-math_I3
      forall (i=1:3,j=1:3)
	         lattice_Stwin(i,j,k,l)=lattice_td(i,k,l)*lattice_tn(j,k,l)
			 lattice_Qtwin(i,j,k,l)=lattice_Qtwin(i,j,k,l)+2*lattice_tn(i,k,l)*lattice_tn(j,k,l)
      endforall
!* Vectorization of normalized Schmid matrix
      lattice_Stwin_v(1,k,l)=lattice_Stwin(1,1,k,l)
      lattice_Stwin_v(2,k,l)=lattice_Stwin(2,2,k,l)
      lattice_Stwin_v(3,k,l)=lattice_Stwin(3,3,k,l)
	  !* be compatible with Mandel notation of Tstar
      lattice_Stwin_v(4,k,l)=(lattice_Stwin(1,2,k,l)+lattice_Stwin(2,1,k,l))/dsqrt(2.0_pReal)
      lattice_Stwin_v(5,k,l)=(lattice_Stwin(2,3,k,l)+lattice_Stwin(3,2,k,l))/dsqrt(2.0_pReal)
      lattice_Stwin_v(6,k,l)=(lattice_Stwin(1,3,k,l)+lattice_Stwin(3,1,k,l))/dsqrt(2.0_pReal)
   enddo
enddo

!*** printout schmid matrix (0nly Hexagonal structure)to check if the conversion is correctly done.

!* define the output location
!open(7, FILE='slip.prn')
!open(8, FILE='twin.prn')
!
!do k = 1,24
!       write(7,*) k
!        write(7,*) lattice_Sslip(1,1,k,3),lattice_Sslip(1,2,k,3),lattice_Sslip(1,3,k,3)
!        write(7,*) lattice_Sslip(2,1,k,3),lattice_Sslip(2,2,k,3),lattice_Sslip(2,3,k,3)
!        write(7,*) lattice_Sslip(3,1,k,3),lattice_Sslip(3,2,k,3),lattice_Sslip(3,3,k,3)
!        write(7,*)
!        write(8,*) k
!        write(8,*) lattice_Stwin(1,1,k,3),lattice_Stwin(2,2,k,3),lattice_Stwin(3,3,k,3)
!        write(8,*) lattice_Stwin(2,1,k,3),lattice_Stwin(2,2,k,3),lattice_Stwin(2,3,k,3)
!        write(8,*) lattice_Stwin(3,1,k,3),lattice_Stwin(3,2,k,3),lattice_Stwin(3,3,k,3)
!        write(8,*)
!enddo

end subroutine

END MODULE

	   
         
