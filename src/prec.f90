!--------------------------------------------------------------------------------------------------
!> @author   Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    setting precision for real and int type
!--------------------------------------------------------------------------------------------------
module prec
  use, intrinsic :: IEEE_arithmetic, only:&
    IEEE_selected_real_kind

  implicit none
  private
  ! https://software.intel.com/en-us/blogs/2017/03/27/doctor-fortran-in-it-takes-all-kinds
#ifdef Abaqus
  integer,     parameter, public :: pReal      = selected_real_kind(15,307)                         !< number with 15 significant digits, up to 1e+-307 (typically 64 bit)
#else
  integer,     parameter, public :: pReal      = IEEE_selected_real_kind(15,307)                    !< number with 15 significant digits, up to 1e+-307 (typically 64 bit)
#endif
#if(INT==8)
  integer,     parameter, public :: pInt       = selected_int_kind(18)                              !< number with at least up to +-1e18 (typically 64 bit)
#else
  integer,     parameter, public :: pInt       = selected_int_kind(9)                               !< number with at least up to +-1e9 (typically 32 bit)
#endif
  integer,     parameter, public :: pLongInt   = selected_int_kind(18)                              !< number with at least up to +-1e18 (typically 64 bit)
  integer,     parameter, public :: pStringLen = 256                                                !< default string length

  real(pReal), parameter, public :: tol_math_check = 1.0e-8_pReal                                   !< tolerance for internal math self-checks (rotation)


  type, public :: group_float                                                                       !< variable length datatype used for storage of state
    real(pReal), dimension(:), pointer :: p
  end type group_float

  type, public :: group_int
    integer, dimension(:), pointer :: p
  end type group_int

  ! http://stackoverflow.com/questions/3948210/can-i-have-a-pointer-to-an-item-in-an-allocatable-array
  type, public :: tState
    integer :: &
      sizeState        = 0, &                                                                       !< size of state
      sizeDotState     = 0, &                                                                       !< size of dot state, i.e. state(1:sizeDot) follows time evolution by dotState rates
      offsetDeltaState = 0, &                                                                       !< index offset of delta state
      sizeDeltaState   = 0, &                                                                       !< size of delta state, i.e. state(offset+1:offset+sizeDelta) follows time evolution by deltaState increments
      sizePostResults  = 0                                                                          !< size of output data
    real(pReal), pointer,     dimension(:), contiguous :: &
      atolState
    real(pReal), pointer,     dimension(:,:), contiguous :: &                                       ! a pointer is needed here because we might point to state/doState. However, they will never point to something, but are rather allocated and, hence, contiguous
      state0, &
      state, &                                                                                      !< state
      dotState, &                                                                                   !< rate of state change
      deltaState                                                                                    !< increment of state change
    real(pReal), allocatable, dimension(:,:) :: &
      partionedState0, &
      subState0, &
      previousDotState, &                                                                           !< state rate of previous xxxx
      previousDotState2, &                                                                          !< state rate two xxxx ago
      RK4dotState
    real(pReal), allocatable, dimension(:,:,:) :: &
      RKCK45dotState
  end type

  type, extends(tState), public :: tPlasticState
    integer :: &
      nSlip  = 0, &
      nTwin  = 0, &
      nTrans = 0
    logical :: &
      nonlocal = .false.
    real(pReal), pointer,     dimension(:,:) :: &
      slipRate, &                                                                                   !< slip rate
      accumulatedSlip                                                                               !< accumulated plastic slip
  end type

  type, public :: tSourceState
    type(tState), dimension(:), allocatable :: p                                                    !< tState for each active source mechanism in a phase
  end type

  type, public :: tHomogMapping
    integer, pointer, dimension(:,:) :: p
  end type

  real(pReal), private, parameter :: PREAL_EPSILON = epsilon(0.0_pReal)                             !< minimum positive number such that 1.0 + EPSILON /= 1.0.
  real(pReal), private, parameter :: PREAL_MIN     = tiny(0.0_pReal)                                !< smallest normalized floating point number

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

  integer, allocatable, dimension(:) :: realloc_lhs_test

  external :: &
    quit

  write(6,'(/,a)') ' <<<+-  prec init  -+>>>'

  write(6,'(a,i3)')    ' Size of integer in bit: ',bit_size(0)
  write(6,'(a,i19)')   '   Maximum value:        ',huge(0)
  write(6,'(/,a,i3)')  ' Size of float in bit:   ',storage_size(0.0_pReal)
  write(6,'(a,e10.3)') '   Maximum value:        ',huge(0.0_pReal)
  write(6,'(a,e10.3)') '   Minimum value:        ',tiny(0.0_pReal)
  write(6,'(a,i3)')    '   Decimal precision:    ',precision(0.0_pReal)

  realloc_lhs_test = [1,2]
  if (realloc_lhs_test(2)/=2) call quit(9000)

end subroutine prec_init


!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for float with double precision
! replaces "==" but for certain (relative) tolerance. Counterpart to dNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq(a,b,tol)

  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol
  real(pReal)                       :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_EPSILON * maxval(abs([a,b]))
  endif

  dEq = merge(.True.,.False.,abs(a-b) < eps)

end function dEq


!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for float with double precision
! replaces "!=" but for certain (relative) tolerance. Counterpart to dEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative NOT
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq(a,b,tol)

  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol
  real(pReal)                       :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_EPSILON * maxval(abs([a,b]))
  endif

  dNeq = merge(.False.,.True.,abs(a-b) <= eps)

end function dNeq


!--------------------------------------------------------------------------------------------------
!> @brief equality to 0 comparison for float with double precision
! replaces "==0" but everything not representable as a normal number is treated as 0. Counterpart to dNeq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq0(a,tol)

  real(pReal), intent(in)           :: a
  real(pReal), intent(in), optional :: tol
  real(pReal)                       :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_MIN * 10.0_pReal
  endif

  dEq0 = merge(.True.,.False.,abs(a) < eps)

end function dEq0


!--------------------------------------------------------------------------------------------------
!> @brief inequality to 0 comparison for float with double precision
! replaces "!=0" but everything not representable as a normal number is treated as 0. Counterpart to dEq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq0(a,tol)

  real(pReal), intent(in)           :: a
  real(pReal), intent(in), optional :: tol
  real(pReal)                       :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_MIN * 10.0_pReal
  endif

  dNeq0 = merge(.False.,.True.,abs(a) <= eps)

end function dNeq0


!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for complex with double precision
! replaces "==" but for certain (relative) tolerance. Counterpart to cNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cEq(a,b,tol)

  complex(pReal), intent(in)           :: a,b
  real(pReal),    intent(in), optional :: tol
  real(pReal)                          :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_EPSILON * maxval(abs([a,b]))
  endif

  cEq = merge(.True.,.False.,abs(a-b) < eps)

end function cEq


!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for complex with double precision
! replaces "!=" but for certain (relative) tolerance. Counterpart to cEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cNeq(a,b,tol)

  complex(pReal), intent(in)           :: a,b
  real(pReal),    intent(in), optional :: tol
  real(pReal)                          :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_EPSILON * maxval(abs([a,b]))
  endif

  cNeq = merge(.False.,.True.,abs(a-b) <= eps)

end function cNeq

end module prec
