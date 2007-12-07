
!************************************
!*         Module: CRYSTAL          *
!************************************
!* contains:                        *
!* - Crystal structure definition   *
!* - Slip system definition         *
!* - Schmid matrices calculation    *
!************************************

MODULE crystal

!*** Include other modules ***
use prec, only: pReal,pInt
implicit none

!************************************
!*      Crystal structures          *
!************************************
!* Number of crystal structures (1-FCC,2-BCC,3-HCP)
integer(pInt), parameter :: crystal_MaxCrystalStructure = 3
!* Total number of slip systems per crystal structure
!* (has to be changed according the definition of slip systems)
integer(pInt), dimension(crystal_MaxCrystalStructure), parameter :: crystal_MaxNslipOfStructure = &
reshape((/12,48,12/),(/crystal_MaxCrystalStructure/))
!* Total number of twin systems per crystal structure
!* (has to be changed according the definition of twin systems)
integer(pInt), dimension(crystal_MaxCrystalStructure), parameter :: crystal_MaxNtwinOfStructure = &
reshape((/12,12,6/),(/crystal_MaxCrystalStructure/))
!* Maximum number of slip systems over crystal structures
integer(pInt), parameter :: crystal_MaxMaxNslipOfStructure = 48
!* Maximum number of twin systems over crystal structures
integer(pInt), parameter :: crystal_MaxMaxNtwinOfStructure = 12
!* Slip direction, slip normales and Schmid matrices
real(pReal), dimension(3,3,crystal_MaxMaxNslipOfStructure,crystal_MaxCrystalStructure) :: crystal_Sslip
real(pReal), dimension(6,crystal_MaxMaxNslipOfStructure,crystal_MaxCrystalStructure) :: crystal_Sslip_v
real(pReal), dimension(3,crystal_MaxMaxNslipOfStructure,crystal_MaxCrystalStructure) :: crystal_sn
real(pReal), dimension(3,crystal_MaxMaxNslipOfStructure,crystal_MaxCrystalStructure) :: crystal_sd
real(pReal), dimension(3,crystal_MaxMaxNslipOfStructure,crystal_MaxCrystalStructure) :: crystal_st
!* twin direction, twin normales, Schmid matrices and transformation matrices
real(pReal), dimension(3,3,crystal_MaxMaxNtwinOfStructure,crystal_MaxCrystalStructure) :: crystal_Stwin
real(pReal), dimension(6,crystal_MaxMaxNtwinOfStructure,crystal_MaxCrystalStructure) :: crystal_Stwin_v
real(pReal), dimension(3,crystal_MaxMaxNtwinOfStructure,crystal_MaxCrystalStructure) :: crystal_tn
real(pReal), dimension(3,crystal_MaxMaxNtwinOfStructure,crystal_MaxCrystalStructure) :: crystal_td
real(pReal), dimension(3,crystal_MaxMaxNtwinOfStructure,crystal_MaxCrystalStructure) :: crystal_tt
real(pReal), dimension(3,3,crystal_MaxMaxNtwinOfStructure,crystal_MaxCrystalStructure) :: crystal_Qtwin

!* Slip_slip interaction matrices
integer(pInt), dimension(crystal_MaxMaxNslipOfStructure,crystal_MaxMaxNslipOfStructure,crystal_MaxCrystalStructure) :: &
crystal_SlipIntType

!*** Slip systems for FCC structures (1) ***
!* System {111}<110>  Sort according Eisenlohr&Hantcherli
data crystal_sd(:, 1,1)/ 0, 1,-1/ ; data crystal_sn(:, 1,1)/ 1, 1, 1/
data crystal_sd(:, 2,1)/-1, 0, 1/ ; data crystal_sn(:, 2,1)/ 1, 1, 1/
data crystal_sd(:, 3,1)/ 1,-1, 0/ ; data crystal_sn(:, 3,1)/ 1, 1, 1/
data crystal_sd(:, 4,1)/ 0,-1,-1/ ; data crystal_sn(:, 4,1)/-1,-1, 1/
data crystal_sd(:, 5,1)/ 1, 0, 1/ ; data crystal_sn(:, 5,1)/-1,-1, 1/
data crystal_sd(:, 6,1)/-1, 1, 0/ ; data crystal_sn(:, 6,1)/-1,-1, 1/
data crystal_sd(:, 7,1)/ 0,-1, 1/ ; data crystal_sn(:, 7,1)/ 1,-1,-1/
data crystal_sd(:, 8,1)/-1, 0,-1/ ; data crystal_sn(:, 8,1)/ 1,-1,-1/
data crystal_sd(:, 9,1)/ 1, 1, 0/ ; data crystal_sn(:, 9,1)/ 1,-1,-1/
data crystal_sd(:,10,1)/ 0, 1, 1/ ; data crystal_sn(:,10,1)/-1, 1,-1/
data crystal_sd(:,11,1)/ 1, 0,-1/ ; data crystal_sn(:,11,1)/-1, 1,-1/
data crystal_sd(:,12,1)/-1,-1, 0/ ; data crystal_sn(:,12,1)/-1, 1,-1/

!*** Twin systems for FCC structures (1) ***
!* System {111}<112>  Sort according Eisenlohr&Hantcherli
data crystal_td(:, 1,1)/-2, 1, 1/ ; data crystal_tn(:, 1,1)/ 1, 1, 1/
data crystal_td(:, 2,1)/ 1,-2, 1/ ; data crystal_tn(:, 2,1)/ 1, 1, 1/
data crystal_td(:, 3,1)/ 1, 1,-2/ ; data crystal_tn(:, 3,1)/ 1, 1, 1/
data crystal_td(:, 4,1)/ 2,-1, 1/ ; data crystal_tn(:, 4,1)/-1,-1, 1/
data crystal_td(:, 5,1)/-1, 2, 1/ ; data crystal_tn(:, 5,1)/-1,-1, 1/
data crystal_td(:, 6,1)/-1,-1,-2/ ; data crystal_tn(:, 6,1)/-1,-1, 1/
data crystal_td(:, 7,1)/-2,-1,-1/ ; data crystal_tn(:, 7,1)/ 1,-1,-1/
data crystal_td(:, 8,1)/ 1, 2,-1/ ; data crystal_tn(:, 8,1)/ 1,-1,-1/
data crystal_td(:, 9,1)/ 1,-1, 2/ ; data crystal_tn(:, 9,1)/ 1,-1,-1/
data crystal_td(:,10,1)/ 2, 1,-1/ ; data crystal_tn(:,10,1)/-1, 1,-1/
data crystal_td(:,11,1)/-1,-2,-1/ ; data crystal_tn(:,11,1)/-1, 1,-1/
data crystal_td(:,12,1)/-1, 1, 2/ ; data crystal_tn(:,12,1)/-1, 1,-1/

!*** Slip-Slip interactions for FCC structures (1) ***
data crystal_SlipIntType( 1,1:crystal_MaxNslipOfStructure(1),1)/1,2,2,4,6,5,3,5,5,4,5,6/
data crystal_SlipIntType( 2,1:crystal_MaxNslipOfStructure(1),1)/2,1,2,6,4,5,5,4,6,5,3,5/
data crystal_SlipIntType( 3,1:crystal_MaxNslipOfStructure(1),1)/2,2,1,5,5,3,5,6,4,6,5,4/
data crystal_SlipIntType( 4,1:crystal_MaxNslipOfStructure(1),1)/4,6,5,1,2,2,4,5,6,3,5,5/
data crystal_SlipIntType( 5,1:crystal_MaxNslipOfStructure(1),1)/6,4,5,2,1,2,5,3,5,5,4,6/
data crystal_SlipIntType( 6,1:crystal_MaxNslipOfStructure(1),1)/5,5,3,2,2,1,6,5,4,5,6,4/
data crystal_SlipIntType( 7,1:crystal_MaxNslipOfStructure(1),1)/3,5,5,4,5,6,1,2,2,4,6,5/
data crystal_SlipIntType( 8,1:crystal_MaxNslipOfStructure(1),1)/5,4,6,5,3,5,2,1,2,6,4,5/
data crystal_SlipIntType( 9,1:crystal_MaxNslipOfStructure(1),1)/5,6,4,6,5,4,2,2,1,5,5,3/
data crystal_SlipIntType(10,1:crystal_MaxNslipOfStructure(1),1)/4,5,6,3,5,5,4,6,5,1,2,2/
data crystal_SlipIntType(11,1:crystal_MaxNslipOfStructure(1),1)/5,3,5,5,4,6,6,4,5,2,1,2/
data crystal_SlipIntType(12,1:crystal_MaxNslipOfStructure(1),1)/6,5,4,5,6,4,5,5,3,2,2,1/

!*** Slip systems for BCC structures (2) ***
!* System {110}<111>
!* Sort?
data crystal_sd(:, 1,2)/ 1,-1, 1/ ; data crystal_sn(:, 1,2)/ 0, 1, 1/
data crystal_sd(:, 2,2)/-1,-1, 1/ ; data crystal_sn(:, 2,2)/ 0, 1, 1/
data crystal_sd(:, 3,2)/ 1, 1, 1/ ; data crystal_sn(:, 3,2)/ 0,-1, 1/
data crystal_sd(:, 4,2)/-1, 1, 1/ ; data crystal_sn(:, 4,2)/ 0,-1, 1/
data crystal_sd(:, 5,2)/-1, 1, 1/ ; data crystal_sn(:, 5,2)/ 1, 0, 1/
data crystal_sd(:, 6,2)/-1,-1, 1/ ; data crystal_sn(:, 6,2)/ 1, 0, 1/
data crystal_sd(:, 7,2)/ 1, 1, 1/ ; data crystal_sn(:, 7,2)/-1, 0, 1/
data crystal_sd(:, 8,2)/ 1,-1, 1/ ; data crystal_sn(:, 8,2)/-1, 0, 1/
data crystal_sd(:, 9,2)/-1, 1, 1/ ; data crystal_sn(:, 9,2)/ 1, 1, 0/
data crystal_sd(:,10,2)/-1, 1,-1/ ; data crystal_sn(:,10,2)/ 1, 1, 0/
data crystal_sd(:,11,2)/ 1, 1, 1/ ; data crystal_sn(:,11,2)/-1, 1, 0/
data crystal_sd(:,12,2)/ 1, 1,-1/ ; data crystal_sn(:,12,2)/-1, 1, 0/
!* System {112}<111>
!* Sort?
data crystal_sd(:,13,2)/-1, 1, 1/ ; data crystal_sn(:,13,2)/ 2, 1, 1/
data crystal_sd(:,14,2)/ 1, 1, 1/ ; data crystal_sn(:,14,2)/-2, 1, 1/
data crystal_sd(:,15,2)/ 1, 1,-1/ ; data crystal_sn(:,15,2)/ 2,-1, 1/
data crystal_sd(:,16,2)/ 1,-1, 1/ ; data crystal_sn(:,16,2)/ 2, 1,-1/
data crystal_sd(:,17,2)/ 1,-1, 1/ ; data crystal_sn(:,17,2)/ 1, 2, 1/
data crystal_sd(:,18,2)/ 1, 1,-1/ ; data crystal_sn(:,18,2)/-1, 2, 1/
data crystal_sd(:,19,2)/ 1, 1, 1/ ; data crystal_sn(:,19,2)/ 1,-2, 1/
data crystal_sd(:,20,2)/-1, 1, 1/ ; data crystal_sn(:,20,2)/ 1, 2,-1/
data crystal_sd(:,21,2)/ 1, 1,-1/ ; data crystal_sn(:,21,2)/ 1, 1, 2/
data crystal_sd(:,22,2)/ 1,-1, 1/ ; data crystal_sn(:,22,2)/-1, 1, 2/
data crystal_sd(:,23,2)/-1, 1, 1/ ; data crystal_sn(:,23,2)/ 1,-1, 2/
data crystal_sd(:,24,2)/ 1, 1, 1/ ; data crystal_sn(:,24,2)/ 1, 1,-2/
!* System {123}<111>
!* Sort?
data crystal_sd(:,25,2)/ 1, 1,-1/ ; data crystal_sn(:,25,2)/ 1, 2, 3/
data crystal_sd(:,26,2)/ 1,-1, 1/ ; data crystal_sn(:,26,2)/-1, 2, 3/
data crystal_sd(:,27,2)/-1, 1, 1/ ; data crystal_sn(:,27,2)/ 1,-2, 3/
data crystal_sd(:,28,2)/ 1, 1, 1/ ; data crystal_sn(:,28,2)/ 1, 2,-3/
data crystal_sd(:,29,2)/ 1,-1, 1/ ; data crystal_sn(:,29,2)/ 1, 3, 2/
data crystal_sd(:,30,2)/ 1, 1,-1/ ; data crystal_sn(:,30,2)/-1, 3, 2/
data crystal_sd(:,31,2)/ 1, 1, 1/ ; data crystal_sn(:,31,2)/ 1,-3, 2/
data crystal_sd(:,32,2)/-1, 1, 1/ ; data crystal_sn(:,32,2)/ 1, 3,-2/
data crystal_sd(:,33,2)/ 1, 1,-1/ ; data crystal_sn(:,33,2)/ 2, 1, 3/
data crystal_sd(:,34,2)/ 1,-1, 1/ ; data crystal_sn(:,34,2)/-2, 1, 3/
data crystal_sd(:,35,2)/-1, 1, 1/ ; data crystal_sn(:,35,2)/ 2,-1, 3/
data crystal_sd(:,36,2)/ 1, 1, 1/ ; data crystal_sn(:,36,2)/ 2, 1,-3/
data crystal_sd(:,37,2)/ 1,-1, 1/ ; data crystal_sn(:,37,2)/ 2, 3, 1/
data crystal_sd(:,38,2)/ 1, 1,-1/ ; data crystal_sn(:,38,2)/-2, 3, 1/
data crystal_sd(:,39,2)/ 1, 1, 1/ ; data crystal_sn(:,39,2)/ 2,-3, 1/
data crystal_sd(:,40,2)/-1, 1, 1/ ; data crystal_sn(:,40,2)/ 2, 3,-1/
data crystal_sd(:,41,2)/-1, 1, 1/ ; data crystal_sn(:,41,2)/ 3, 1, 2/
data crystal_sd(:,42,2)/ 1, 1, 1/ ; data crystal_sn(:,42,2)/-3, 1, 2/
data crystal_sd(:,43,2)/ 1, 1,-1/ ; data crystal_sn(:,43,2)/ 3,-1, 2/
data crystal_sd(:,44,2)/ 1,-1, 1/ ; data crystal_sn(:,44,2)/ 3, 1,-2/
data crystal_sd(:,45,2)/-1, 1, 1/ ; data crystal_sn(:,45,2)/ 3, 2, 1/
data crystal_sd(:,46,2)/ 1, 1, 1/ ; data crystal_sn(:,46,2)/-3, 2, 1/
data crystal_sd(:,47,2)/ 1, 1,-1/ ; data crystal_sn(:,47,2)/ 3,-2, 1/
data crystal_sd(:,48,2)/ 1,-1, 1/ ; data crystal_sn(:,48,2)/ 3, 2,-1/

!*** Twin systems for BCC structures (2) ***
!* System {112}<111>
!* Sort?
!* Not implemented yet

!*** Slip-Slip interactions for BCC structures (2) ***
data crystal_SlipIntType( 1,:,2)/1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 2,:,2)/2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 3,:,2)/2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 4,:,2)/2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 5,:,2)/2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 6,:,2)/2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 7,:,2)/2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 8,:,2)/2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 9,:,2)/2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(10,:,2)/2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(11,:,2)/2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(12,:,2)/2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(13,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(14,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(15,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(16,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(17,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(18,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(19,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(20,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(21,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(22,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(23,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(24,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(25,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(26,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(27,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(28,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(29,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(30,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(31,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(32,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(33,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(34,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(35,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(36,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(37,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(38,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(39,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(40,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2/
data crystal_SlipIntType(41,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2/
data crystal_SlipIntType(42,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2/
data crystal_SlipIntType(43,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2/
data crystal_SlipIntType(44,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2/
data crystal_SlipIntType(45,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2/
data crystal_SlipIntType(46,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2/
data crystal_SlipIntType(47,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2/
data crystal_SlipIntType(48,:,2)/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1/

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
data crystal_sd(:, 1,3)/-1, 0, 0/ ; data crystal_sn(:, 1,3)/ 0, 0, 1/
data crystal_sd(:, 2,3)/ 0,-1, 0/ ; data crystal_sn(:, 2,3)/ 0, 0, 1/
data crystal_sd(:, 3,3)/ 1, 1, 0/ ; data crystal_sn(:, 3,3)/ 0, 0, 1/
!* 1st type prismatic systems {1010}<1120>  (independent of c/a-ratio)
!* 1- ( 0  1 -1  0)[-2  1  1  0]
!* 2- ( 1  0 -1  0)[ 1 -2  1  0]
!* 3- (-1  1  0  0)[ 1  1 -2  0]
!* Sort?
data crystal_sd(:, 4,3)/-1, 0, 0/ ; data crystal_sn(:, 4,3)/ 0, 1, 0/
data crystal_sd(:, 5,3)/ 0,-1, 0/ ; data crystal_sn(:, 5,3)/ 1, 0, 0/
data crystal_sd(:, 6,3)/ 1, 1, 0/ ; data crystal_sn(:, 6,3)/-1, 1, 0/
!* 1st type 1st order pyramidal systems {1011}<1120>
!* plane normales depend on the c/a-ratio
!* 1- ( 0 -1  1  1)[-2  1  1  0]
!* 2- ( 0  1 -1  1)[-2  1  1  0]
!* 3- (-1  0  1  1)[ 1 -2  1  0]
!* 4- ( 1  0 -1  1)[ 1 -2  1  0]
!* 5- (-1  1  0  1)[ 1  1 -2  0]
!* 6- ( 1 -1  0  1)[ 1  1 -2  0]
!* Sort?
data crystal_sd(:, 7,3)/-1, 0, 0/ ; data crystal_sn(:, 7,3)/ 0,-1, 1/
data crystal_sd(:, 8,3)/ 0,-1, 0/ ; data crystal_sn(:, 8,3)/ 0, 1, 1/
data crystal_sd(:, 9,3)/ 1, 1, 0/ ; data crystal_sn(:, 9,3)/-1, 0, 1/
data crystal_sd(:,10,3)/-1, 0, 0/ ; data crystal_sn(:,10,3)/ 1, 0, 1/
data crystal_sd(:,11,3)/ 0,-1, 0/ ; data crystal_sn(:,11,3)/-1, 1, 1/
data crystal_sd(:,12,3)/ 1, 1, 0/ ; data crystal_sn(:,12,3)/ 1,-1, 1/

!*** Twin systems for HCP structures (2) ***
!* System {1012}<1011>
!* Sort?
!* Not implemented yet

!*** Slip-Slip interactions for HCP structures (3) ***
data crystal_SlipIntType( 1,1:crystal_MaxNslipOfStructure(3),3)/1,2,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 2,1:crystal_MaxNslipOfStructure(3),3)/2,1,2,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 3,1:crystal_MaxNslipOfStructure(3),3)/2,2,1,2,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 4,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,1,2,2,2,2,2,2,2,2/
data crystal_SlipIntType( 5,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,1,2,2,2,2,2,2,2/
data crystal_SlipIntType( 6,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,1,2,2,2,2,2,2/
data crystal_SlipIntType( 7,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,1,2,2,2,2,2/
data crystal_SlipIntType( 8,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,1,2,2,2,2/
data crystal_SlipIntType( 9,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,1,2,2,2/
data crystal_SlipIntType(10,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,1,2,2/
data crystal_SlipIntType(11,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,1,2/
data crystal_SlipIntType(12,1:crystal_MaxNslipOfStructure(3),3)/2,2,2,2,2,2,2,2,2,2,2,1/


CONTAINS
!****************************************
!* - crystal_Init
!* - crystal_SchmidMatrices
!****************************************


subroutine crystal_Init()
!**************************************
!*      Module initialization         *
!**************************************
call crystal_SchmidMatrices()
end subroutine


subroutine crystal_SchmidMatrices()
!**************************************
!*   Calculation of Schmid matrices   *
!**************************************
use prec, only: pReal,pInt
use math, only: math_identity2nd
implicit none

!* Definition of variables
integer(pInt) i,j,k,l
real(pReal) invNorm

!* Iteration over the crystal structures
do l=1,crystal_MaxCrystalStructure
!* Iteration over the slip systems
   do k=1,crystal_MaxNslipOfStructure(l)
!* Definition of transverse direction st for the frame (sd,st,sn)
      crystal_st(1,k,l)=crystal_sn(2,k,l)*crystal_sd(3,k,l)-crystal_sn(3,k,l)*crystal_sd(2,k,l)
	  crystal_st(2,k,l)=crystal_sn(3,k,l)*crystal_sd(1,k,l)-crystal_sn(1,k,l)*crystal_sd(3,k,l)
	  crystal_st(3,k,l)=crystal_sn(1,k,l)*crystal_sd(2,k,l)-crystal_sn(2,k,l)*crystal_sd(1,k,l)
!* Defintion of Schmid matrix
      forall (i=1:3,j=1:3)
	         crystal_Sslip(i,j,k,l)=crystal_sd(i,k,l)*crystal_sn(j,k,l)
      endforall
!* Normalization of Schmid matrix
      invNorm=dsqrt(1.0_pReal/((crystal_sn(1,k,l)**2+crystal_sn(2,k,l)**2+crystal_sn(3,k,l)**2)*&
	          (crystal_sd(1,k,l)**2+crystal_sd(2,k,l)**2+crystal_sd(3,k,l)**2)))
      crystal_Sslip(:,:,k,l)=crystal_Sslip(:,:,k,l)*invNorm
!* Vectorization of normalized Schmid matrix
      crystal_Sslip_v(1,k,l)=crystal_Sslip(1,1,k,l)
      crystal_Sslip_v(2,k,l)=crystal_Sslip(2,2,k,l)
      crystal_Sslip_v(3,k,l)=crystal_Sslip(3,3,k,l)
	  !* be compatible with Mandel notation of Tstar
      crystal_Sslip_v(4,k,l)=(crystal_Sslip(1,2,k,l)+crystal_Sslip(2,1,k,l))/dsqrt(2.0_pReal)
      crystal_Sslip_v(5,k,l)=(crystal_Sslip(2,3,k,l)+crystal_Sslip(3,2,k,l))/dsqrt(2.0_pReal)
      crystal_Sslip_v(6,k,l)=(crystal_Sslip(1,3,k,l)+crystal_Sslip(3,1,k,l))/dsqrt(2.0_pReal)
   enddo

!* Iteration over the twin systems
   do k=1,crystal_MaxNslipOfStructure(l)
!* Definition of transverse direction tt for the frame (td,tt,tn)
      crystal_tt(1,k,l)=crystal_tn(2,k,l)*crystal_td(3,k,l)-crystal_tn(3,k,l)*crystal_td(2,k,l)
	  crystal_tt(2,k,l)=crystal_tn(3,k,l)*crystal_td(1,k,l)-crystal_tn(1,k,l)*crystal_td(3,k,l)
	  crystal_tt(3,k,l)=crystal_tn(1,k,l)*crystal_td(2,k,l)-crystal_tn(2,k,l)*crystal_td(1,k,l)
!* Defintion of Schmid matrix and transformation matrices
      crystal_Qtwin(:,:,k,l)=-math_identity2nd(3)
      forall (i=1:3,j=1:3)
	         crystal_Stwin(i,j,k,l)=crystal_td(i,k,l)*crystal_tn(j,k,l)
			 crystal_Qtwin(i,j,k,l)=crystal_Qtwin(i,j,k,l)+2*crystal_tn(i,k,l)*crystal_tn(j,k,l)
      endforall
!* Normalization of Schmid matrix
      invNorm=dsqrt(1.0_pReal/((crystal_tn(1,k,l)**2+crystal_tn(2,k,l)**2+crystal_tn(3,k,l)**2)*&
	          (crystal_td(1,k,l)**2+crystal_td(2,k,l)**2+crystal_td(3,k,l)**2)))
      crystal_Stwin(:,:,k,l)=crystal_Stwin(:,:,k,l)*invNorm
!* Vectorization of normalized Schmid matrix
      crystal_Stwin_v(1,k,l)=crystal_Stwin(1,1,k,l)
      crystal_Stwin_v(2,k,l)=crystal_Stwin(2,2,k,l)
      crystal_Stwin_v(3,k,l)=crystal_Stwin(3,3,k,l)
	  !* be compatible with Mandel notation of Tstar
      crystal_Stwin_v(4,k,l)=(crystal_Stwin(1,2,k,l)+crystal_Stwin(2,1,k,l))/dsqrt(2.0_pReal)
      crystal_Stwin_v(5,k,l)=(crystal_Stwin(2,3,k,l)+crystal_Stwin(3,2,k,l))/dsqrt(2.0_pReal)
      crystal_Stwin_v(6,k,l)=(crystal_Stwin(1,3,k,l)+crystal_Stwin(3,1,k,l))/dsqrt(2.0_pReal)
   enddo
enddo

end subroutine

END MODULE

	   
         
