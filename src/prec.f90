!--------------------------------------------------------------------------------------------------
!> @author   Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    setting precision for real and int type
!--------------------------------------------------------------------------------------------------
module prec
 implicit none
 private 
#if (FLOAT==8)
 integer,     parameter, public :: pReal = 8                                                        !< floating point double precision (was selected_real_kind(15,300), number with 15 significant digits, up to 1e+-300)
#else
 NO SUITABLE PRECISION FOR REAL SELECTED, STOPPING COMPILATION
#endif

#if (INT==4)
 integer,     parameter, public :: pInt  = 4                                                        !< integer representation 32 bit (was selected_int_kind(9), number with at least up to +- 1e9)
#elif (INT==8)
 integer,     parameter, public :: pInt  = 8                                                        !< integer representation 64 bit (was selected_int_kind(12), number with at least up to +- 1e12)
#else
 NO SUITABLE PRECISION FOR INTEGER SELECTED, STOPPING COMPILATION
#endif

 integer,     parameter, public :: pLongInt  = 8                                                    !< integer representation 64 bit (was selected_int_kind(12), number with at least up to +- 1e12)
 real(pReal), parameter, public :: tol_math_check = 1.0e-8_pReal                                    !< tolerance for internal math self-checks (rotation)

 integer(pInt), allocatable, dimension(:) :: realloc_lhs_test
 
 type, public :: p_vec                                                                              !< variable length datatype used for storage of state
   real(pReal), dimension(:), pointer :: p
 end type p_vec

 type, public :: p_intvec
   integer(pInt), dimension(:), pointer :: p
 end type p_intvec

!http://stackoverflow.com/questions/3948210/can-i-have-a-pointer-to-an-item-in-an-allocatable-array
 type, public :: tState
   integer(pInt) :: &
     sizeState        = 0_pInt, &                                                                   !< size of state
     sizeDotState     = 0_pInt, &                                                                   !< size of dot state, i.e. state(1:sizeDot) follows time evolution by dotState rates
     offsetDeltaState = 0_pInt, &                                                                   !< offset of delta state
     sizeDeltaState   = 0_pInt, &                                                                   !< size of delta state, i.e. state(offset+1:offset+sizeDot) follows time evolution by deltaState increments
     sizePostResults  = 0_pInt                                                                      !< size of output data
   real(pReal), pointer,     dimension(:), contiguous :: &
     atolState
   real(pReal), pointer,     dimension(:,:), contiguous :: &                                        ! a pointer is needed here because we might point to state/doState. However, they will never point to something, but are rather allocated and, hence, contiguous 
     state0, &
     state, &                                                                                       !< state
     dotState, &                                                                                    !< rate of state change
     deltaState                                                                                     !< increment of state change
   real(pReal), allocatable, dimension(:,:) :: &
     partionedState0, &
     subState0, &
     previousDotState, &                                                                            !< state rate of previous xxxx
     previousDotState2, &                                                                           !< state rate two xxxx ago
     RK4dotState
   real(pReal), allocatable, dimension(:,:,:) :: &
     RKCK45dotState
 end type

 type, extends(tState), public :: tPlasticState
   integer(pInt) :: &
     nSlip = 0_pInt , &
     nTwin = 0_pInt, &
     nTrans = 0_pInt
   logical :: & 
     nonlocal = .false.
   real(pReal), pointer,     dimension(:,:), contiguous :: &
     slipRate, &                                                                                    !< slip rate
     accumulatedSlip                                                                                !< accumulated plastic slip
 end type

 type, public :: tSourceState
   type(tState), dimension(:), allocatable :: p                                                     !< tState for each active source mechanism in a phase
 end type
 
 type, public :: tHomogMapping
   integer(pInt), pointer, dimension(:,:) :: p                                  
 end type 

 type, public :: tPhaseMapping
   integer(pInt), pointer, dimension(:,:,:) :: p
 end type 

#ifdef FEM
 type, public :: tOutputData
   integer(pInt) :: &
     sizeIpCells = 0_pInt , &
     sizeResults = 0_pInt
   real(pReal), allocatable, dimension(:,:) :: &
     output                                                                                         !< output data
 end type 
#endif

 public :: &
   prec_init, &
   dEq, &
   dEq0, &
   cEq, &
   dNeq, &
   dNeq0, &
   cNeq
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief reporting precision
!--------------------------------------------------------------------------------------------------
subroutine prec_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif

 implicit none
 external :: &
   quit

 write(6,'(/,a)') ' <<<+-  prec init  -+>>>'
#include "compilation_info.f90"
 write(6,'(a,i3)')    ' Bytes for pReal:    ',pReal
 write(6,'(a,i3)')    ' Bytes for pInt:     ',pInt
 write(6,'(a,i3)')    ' Bytes for pLongInt: ',pLongInt

 realloc_lhs_test = [1_pInt,2_pInt]
 if (realloc_lhs_test(2)/=2_pInt) call quit(9000)
 
end subroutine prec_init


!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for float with double precision
! replaces "==" but for certain (relative) tolerance. Counterpart to dNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq(a,b,tol)

 implicit none
 real(pReal), intent(in)           :: a,b
 real(pReal), intent(in), optional :: tol
 real(pReal), parameter            :: eps = 2.220446049250313E-16                                   ! DBL_EPSILON in C

 dEq = merge(.True., .False.,abs(a-b) <= merge(tol,eps,present(tol))*maxval(abs([a,b])))
end function dEq


!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for float with double precision
! replaces "!=" but for certain (relative) tolerance. Counterpart to dEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative NOT
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq(a,b,tol)

 implicit none
 real(pReal), intent(in)           :: a,b
 real(pReal), intent(in), optional :: tol
 real(pReal), parameter            :: eps = 2.220446049250313E-16                                   ! DBL_EPSILON in C

 dNeq = merge(.False., .True.,abs(a-b) <= merge(tol,eps,present(tol))*maxval(abs([a,b])))
end function dNeq


!--------------------------------------------------------------------------------------------------
!> @brief equality to 0 comparison for float with double precision
! replaces "==0" but everything not representable as a normal number is treated as 0. Counterpart to dNeq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq0(a,tol)

 implicit none
 real(pReal), intent(in)           :: a
 real(pReal), intent(in), optional :: tol
 real(pReal), parameter            :: eps = 2.2250738585072014E-308                                 ! smallest non-denormalized number

 dEq0 = merge(.True., .False.,abs(a) <= merge(tol,eps,present(tol)))
end function dEq0


!--------------------------------------------------------------------------------------------------
!> @brief inequality to 0 comparison for float with double precision
! replaces "!=0" but everything not representable as a normal number is treated as 0. Counterpart to dEq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq0(a,tol)

 implicit none
 real(pReal), intent(in)           :: a
 real(pReal), intent(in), optional :: tol
 real(pReal), parameter            :: eps = 2.2250738585072014E-308                                 ! smallest non-denormalized number

 dNeq0 = merge(.False., .True.,abs(a) <= merge(tol,eps,present(tol)))
end function dNeq0


!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for complex with double precision
! replaces "==" but for certain (relative) tolerance. Counterpart to cNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cEq(a,b,tol)

 implicit none
 complex(pReal), intent(in)           :: a,b
 real(pReal),    intent(in), optional :: tol
 real(pReal),    parameter            :: eps = 2.220446049250313E-16                                ! DBL_EPSILON in C

 cEq = merge(.True., .False.,abs(a-b) <= merge(tol,eps,present(tol))*maxval(abs([a,b])))
end function cEq


!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for complex with double precision
! replaces "!=" but for certain (relative) tolerance. Counterpart to cEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cNeq(a,b,tol)

 implicit none
 complex(pReal), intent(in)           :: a,b
 real(pReal),    intent(in), optional :: tol
 real(pReal),    parameter            :: eps = 2.220446049250313E-16                                ! DBL_EPSILON in C

 cNeq = merge(.False., .True.,abs(a-b) <= merge(tol,eps,present(tol))*maxval(abs([a,b])))
end function cNeq

end module prec
