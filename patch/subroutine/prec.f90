
!##############################################################
 MODULE prec
!##############################################################

 implicit none

!    *** Precision of real and integer variables ***
 integer, parameter :: pReal = selected_real_kind(15,300)  ! 15 significant digits, up to 1e+-300
 integer, parameter :: pInt  = selected_int_kind(9)        ! up to +- 1e9
 integer, parameter :: pLongInt  = 8                       ! should be 64bit


 type :: p_vec
     real(pReal), dimension(:), pointer :: p
 end type p_vec
 
!    *** Strain increment considered significant ***
 real(pReal), parameter :: relevantStrain = 1.0e-7_pReal

!    *** Numerical parameters ***
 integer(pInt), parameter :: ijaco        = 1_pInt       ! frequency of FEM Jacobi update
 integer(pInt), parameter :: nCutback     = 20_pInt      ! max cutbacks accounted for in debug distribution
 integer(pInt), parameter :: nReg         = 1_pInt       ! regularization attempts for Jacobi inversion
 real(pReal),   parameter :: pert_Fg      = 1.0e-6_pReal ! strain perturbation for FEM Jacobi
 integer(pInt), parameter :: nOuter       = 20_pInt      ! outer loop limit 20
 integer(pInt), parameter :: nInner       = 200_pInt     ! inner loop limit 200
 real(pReal),   parameter :: reltol_Outer = 1.0e-5_pReal ! relative tolerance in outer loop (state)
 real(pReal),   parameter :: reltol_Inner = 1.0e-6_pReal ! relative tolerance in inner loop (Lp)
 real(pReal),   parameter :: abstol_Inner = 1.0e-8_pReal ! absolute tolerance in inner loop (Lp)
!
 real(pReal),   parameter :: resToler = 1.0e-4_pReal     ! relative tolerance of residual in GIA iteration
 real(pReal),   parameter :: resAbsol = 1.0e+2_pReal     ! absolute tolerance of residual in GIA iteration (corresponds to ~1 Pa)
 real(pReal),   parameter :: resBound = 1.0e+1_pReal     ! relative maximum value (upper bound) for GIA residual
 integer(pInt), parameter :: NRiterMax = 24_pInt         ! maximum number of GIA iteration
 real(pReal),   parameter :: subStepMin = 1.0e-3_pReal   ! minimum (relative) size of sub-step allowed during cutback

 END MODULE prec
