!$DAMASK:copyright.f90
!##############################################################
!* $Id$
!##############################################################

MODULE prec
 implicit none 
!    *** Precision of real and integer variables for python interfacing***
 integer, parameter :: pReal = 8
 integer, parameter :: pInt  = 4
 real(pReal), parameter :: DAMASK_NaN = real(Z'7FF0000000000001',pReal)
 real(pReal), parameter :: tol_math_check = 1.0e-8_pReal
END MODULE prec

MODULE debug
 use prec, only: pInt
 implicit none
 integer(pInt), parameter :: debug_verbosity = 1_pInt
END MODULE debug

MODULE numerics
 use prec, only: pInt, pReal
 implicit none
 real(pReal), parameter :: fftw_timelimit = -1.0_pReal
 integer(pInt), parameter :: fftw_planner_flag = 32_pInt
 integer(pInt), parameter :: fixedSeed = 1_pInt
END MODULE numerics

MODULE IO
 CONTAINS
 subroutine IO_error(error_ID,e,i,g,ext_msg)

 use prec, only: pInt
 implicit none
 integer(pInt), intent(in) :: error_ID
 integer(pInt), optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 character(len=1024) msg

 select case (error_ID)
   case default
     print*, 'Error messages not supported when interfacing to Python'
 end select
 end subroutine IO_error

END MODULE IO
