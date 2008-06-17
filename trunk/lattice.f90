
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
reshape((/12,48,12/),(/lattice_MaxLatticeStructure/))
!* Total number of twin systems per lattice structure
!* (has to be changed according the definition of twin systems)
integer(pInt), dimension(lattice_MaxLatticeStructure), parameter :: lattice_MaxNtwinOfStructure = &
reshape((/12,0,0/),(/lattice_MaxLatticeStructure/))
!* Maximum number of slip systems over lattice structures
integer(pInt), parameter :: lattice_MaxMaxNslipOfStructure = 48
!* Maximum number of twin systems over lattice structures
integer(pInt), parameter :: lattice_MaxMaxNtwinOfStructure = 12
!* Slip direction, slip normales and Schmid matrices
real(pReal), dimension(3,3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_Sslip
real(pReal), dimension(6,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_Sslip_v
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_sn
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_sd
real(pReal), dimension(3,lattice_MaxMaxNslipOfStructure,lattice_MaxLatticeStructure) :: lattice_st
!* twin direction, twin normales, Schmid matrices and transformation matrices
real(pReal), dimension(3,3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_Stwin
real(pReal), dimension(6,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_Stwin_v
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_tn
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_td
real(pReal), dimension(3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_tt
real(pReal), dimension(3,3,lattice_MaxMaxNtwinOfStructure,lattice_MaxLatticeStructure) :: lattice_Qtwin
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
!* 1- (0 0 0 1)[-2  1  1  0]
!* 2- (0 0 0 1)[ 1 -2  1  0]
!* 3- (0 0 0 1)[ 1  1 -2  0]
!* Plane (hkil)->(hkl)
!* Direction [uvtw]->[(u-t) (v-t) w]
!* Automatical transformation from Bravais to Miller
!* not done for the moment
!* Sort?
data lattice_sd(:, 1,3)/-1, 0, 0/ ; data lattice_sn(:, 1,3)/ 0, 0, 1/
data lattice_sd(:, 2,3)/ 0,-1, 0/ ; data lattice_sn(:, 2,3)/ 0, 0, 1/
data lattice_sd(:, 3,3)/ 1, 1, 0/ ; data lattice_sn(:, 3,3)/ 0, 0, 1/
!* 1st type prismatic systems {1010}<1120>  (independent of c/a-ratio)
!* 1- ( 0  1 -1  0)[-2  1  1  0]
!* 2- ( 1  0 -1  0)[ 1 -2  1  0]
!* 3- (-1  1  0  0)[ 1  1 -2  0]
!* Sort?
data lattice_sd(:, 4,3)/-1, 0, 0/ ; data lattice_sn(:, 4,3)/ 0, 1, 0/
data lattice_sd(:, 5,3)/ 0,-1, 0/ ; data lattice_sn(:, 5,3)/ 1, 0, 0/
data lattice_sd(:, 6,3)/ 1, 1, 0/ ; data lattice_sn(:, 6,3)/-1, 1, 0/
!* 1st type 1st order pyramidal systems {1011}<1120>
!* plane normales depend on the c/a-ratio
!* 1- ( 0 -1  1  1)[-2  1  1  0]
!* 2- ( 0  1 -1  1)[-2  1  1  0]
!* 3- (-1  0  1  1)[ 1 -2  1  0]
!* 4- ( 1  0 -1  1)[ 1 -2  1  0]
!* 5- (-1  1  0  1)[ 1  1 -2  0]
!* 6- ( 1 -1  0  1)[ 1  1 -2  0]
!* Sort?
data lattice_sd(:, 7,3)/-1, 0, 0/ ; data lattice_sn(:, 7,3)/ 0,-1, 1/
data lattice_sd(:, 8,3)/ 0,-1, 0/ ; data lattice_sn(:, 8,3)/ 0, 1, 1/
data lattice_sd(:, 9,3)/ 1, 1, 0/ ; data lattice_sn(:, 9,3)/-1, 0, 1/
data lattice_sd(:,10,3)/-1, 0, 0/ ; data lattice_sn(:,10,3)/ 1, 0, 1/
data lattice_sd(:,11,3)/ 0,-1, 0/ ; data lattice_sn(:,11,3)/-1, 1, 1/
data lattice_sd(:,12,3)/ 1, 1, 0/ ; data lattice_sn(:,12,3)/ 1,-1, 1/

!*** Twin systems for HCP structures (3) ***
!* System {1012}<1011>
!* Sort?
!* MISSING: not implemented yet

!*** Slip-Slip interactions for HCP structures (3) ***
data lattice_SlipIntType( 1,1:lattice_MaxNslipOfStructure(3),3)/1,2,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 2,1:lattice_MaxNslipOfStructure(3),3)/2,1,2,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 3,1:lattice_MaxNslipOfStructure(3),3)/2,2,1,2,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 4,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,1,2,2,2,2,2,2,2,2/
data lattice_SlipIntType( 5,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,1,2,2,2,2,2,2,2/
data lattice_SlipIntType( 6,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,1,2,2,2,2,2,2/
data lattice_SlipIntType( 7,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,1,2,2,2,2,2/
data lattice_SlipIntType( 8,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,1,2,2,2,2/
data lattice_SlipIntType( 9,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,1,2,2,2/
data lattice_SlipIntType(10,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,1,2,2/
data lattice_SlipIntType(11,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,1,2/
data lattice_SlipIntType(12,1:lattice_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,1/

!*** slip-twin interactions for HCP structures (3) ***
! MISSING: not implemented yet

!*** Twin-twin interactions for HCP structures (3) ***
! MISSING: not implemented yet


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

end subroutine

END MODULE

	   
         
