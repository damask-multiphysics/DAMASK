!* $Id$
!##############################################################
 MODULE prec
!##############################################################

 implicit none
 
!    *** Precision of real and integer variables ***
 integer, parameter :: pReal = selected_real_kind(15,300)     ! 15 significant digits, up to 1e+-300
 integer, parameter :: pInt  = selected_int_kind(9)           ! up to +- 1e9
 integer, parameter :: pLongInt  = 8                          ! should be 64bit
 real(pReal), parameter :: tol_math_check = 1.0e-8_pReal
 real(pReal), parameter :: tol_gravityNodePos = 1.0e-100_pReal
 
 type :: p_vec
     real(pReal), dimension(:), pointer :: p
 end type p_vec

CONTAINS

subroutine prec_init
 implicit none

 write(6,*)
 write(6,*) '<<<+-  prec init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 return
end subroutine


 END MODULE prec
