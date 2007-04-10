
!##############################################################
 MODULE prec
!##############################################################

 implicit none
!    *** Precision of real and integer variables ***
 integer, parameter :: pReal = 8
 integer, parameter :: pInt  = 4
!    *** Numerical parameters ***
!    *** How frequently the jacobian is recalculated ***
 integer (pInt), parameter :: ijaco = 5_pInt
!    *** Maximum number of internal cutbacks in time step ***
 integer(pInt), parameter :: ncut=7_pInt
!    *** Maximum number of regularization attempts for Jacobi inversion ***
 integer(pInt), parameter :: nreg=1_pInt
!    *** Perturbation of strain array for numerical calculation of FEM Jacobi matrix ***
 real(pReal), parameter :: pert_e=1.0e-5_pReal  
!    *** Maximum number of iterations in outer (statevariables) loop ***
 integer(pInt), parameter :: nouter    = 50_pInt
!    *** Convergence criteria for outer (statevariables) loop ***
 real(pReal),   parameter :: tol_outer = 1.0e-4_pReal
!    *** Maximum number of iterations in inner (stress) loop ***
 integer(pInt), parameter :: ninner    = 2000_pInt
!    *** Convergence criteria for inner (stress) loop ***
 real(pReal),   parameter :: tol_inner = 1.0e-3_pReal   
!    *** Factor for maximum stress correction in inner (stress) loop ***
 real(pReal),   parameter :: crite     = 1.0e-1_pReal 

 END MODULE prec
