
!##############################################################
 MODULE prec
!##############################################################

 implicit none

!    *** Precision of real and integer variables ***
 integer, parameter :: pReal = 8
 integer, parameter :: pInt  = 4

!    *** Strain increment considered significant ***
 real(pReal), parameter :: relevantStrain = 1.0e-7_pReal

!    *** Numerical parameters ***
 integer(pInt), parameter :: ijaco        = 1_pInt       ! frequency of FEM Jacobi update
 integer(pInt), parameter :: nCutback     = 10_pInt      ! cutbacks in time-step integration
 integer(pInt), parameter :: nReg         = 1_pInt       ! regularization attempts for Jacobi inversion
 real(pReal),   parameter :: pert_Fg      = 1.0e-5_pReal ! strain perturbation for FEM Jacobi
 integer(pInt), parameter :: nOuter       = 10_pInt      ! outer loop limit
 integer(pInt), parameter :: nInner       = 1000_pInt    ! inner loop limit
 real(pReal),   parameter :: reltol_Outer = 1.0e-4_pReal ! relative tolerance in outer loop (state)
 real(pReal),   parameter :: reltol_Inner = 1.0e-6_pReal ! relative tolerance in inner loop (Lp)
 real(pReal),   parameter :: abstol_Inner = 1.0e-8_pReal ! absolute tolerance in inner loop (Lp)

 END MODULE prec
