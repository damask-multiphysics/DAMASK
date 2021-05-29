!--------------------------------------------------------------------------------------------------
!> @author   Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    setting precision for real and int type
!--------------------------------------------------------------------------------------------------
module prec
  use, intrinsic :: IEEE_arithmetic
  use, intrinsic :: ISO_C_binding

  implicit none
  public

  ! https://software.intel.com/en-us/blogs/2017/03/27/doctor-fortran-in-it-takes-all-kinds
  integer,     parameter :: pReal      = IEEE_selected_real_kind(15,307)                            !< number with 15 significant digits, up to 1e+-307 (typically 64 bit)
  integer,     parameter :: pI32       = selected_int_kind(9)                                       !< number with at least up to +-1e9 (typically 32 bit)
  integer,     parameter :: pI64       = selected_int_kind(18)                                      !< number with at least up to +-1e18 (typically 64 bit)
#if(INT==8)
  integer,     parameter :: pInt       = pI64
#else
  integer,     parameter :: pInt       = pI32
#endif
  integer,     parameter :: pStringLen = 256                                                        !< default string length
  integer,     parameter :: pPathLen   = 4096                                                       !< maximum length of a path name on linux

  real(pReal), parameter :: tol_math_check = 1.0e-8_pReal                                           !< tolerance for internal math self-checks (rotation)


  type :: group_float                                                                               !< variable length datatype used for storage of state
    real(pReal), dimension(:), pointer :: p
  end type group_float

  type :: tState
    integer :: &
      sizeState        = 0, &                                                                       !< size of state
      sizeDotState     = 0, &                                                                       !< size of dot state, i.e. state(1:sizeDot) follows time evolution by dotState rates
      offsetDeltaState = 0, &                                                                       !< index offset of delta state
      sizeDeltaState   = 0                                                                          !< size of delta state, i.e. state(offset+1:offset+sizeDelta) follows time evolution by deltaState increments
    ! http://stackoverflow.com/questions/3948210
    real(pReal), pointer,     dimension(:), contiguous :: &
      atol
    real(pReal), pointer,     dimension(:,:), contiguous :: &                                       ! a pointer is needed here because we might point to state/doState. However, they will never point to something, but are rather allocated and, hence, contiguous
      state0, &
      state, &                                                                                      !< state
      dotState, &                                                                                   !< rate of state change
      deltaState                                                                                    !< increment of state change
    real(pReal), allocatable, dimension(:,:) :: &
      subState0
  end type

  type, extends(tState) :: tPlasticState
    logical :: &
      nonlocal = .false.
    real(pReal), pointer, dimension(:,:) :: &
      slipRate                                                                                      !< slip rate
  end type

  type :: tSourceState
    type(tState), dimension(:), allocatable :: p                                                    !< tState for each active source mechanism in a phase
  end type

  real(pReal), private, parameter :: PREAL_EPSILON = epsilon(0.0_pReal)                             !< minimum positive number such that 1.0 + EPSILON /= 1.0.
  real(pReal), private, parameter :: PREAL_MIN     = tiny(0.0_pReal)                                !< smallest normalized floating point number

  integer,                   dimension(0), parameter :: &
    emptyIntArray    = [integer::]
  real(pReal),               dimension(0), parameter :: &
    emptyRealArray   = [real(pReal)::]
  character(len=pStringLen), dimension(0), parameter :: &
    emptyStringArray = [character(len=pStringLen)::]

  private :: &
    selfTest

contains


!--------------------------------------------------------------------------------------------------
!> @brief Report precision and do self test.
!--------------------------------------------------------------------------------------------------
subroutine prec_init

  print'(/,a)', ' <<<+-  prec init  -+>>>'

  print'(a,i3)',    ' Size of integer in bit: ',bit_size(0)
  print'(a,i19)',   '   Maximum value:        ',huge(0)
  print'(/,a,i3)',  ' Size of float in bit:   ',storage_size(0.0_pReal)
  print'(a,e10.3)', '   Maximum value:        ',huge(0.0_pReal)
  print'(a,e10.3)', '   Minimum value:        ',PREAL_MIN
  print'(a,e10.3)', '   Epsilon value:        ',PREAL_EPSILON
  print'(a,i3)',    '   Decimal precision:    ',precision(0.0_pReal)

  call selfTest

end subroutine prec_init


!--------------------------------------------------------------------------------------------------
!> @brief Test floating point numbers with double precision for equality.
! replaces "==" but for certain (relative) tolerance. Counterpart to dNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq(a,b,tol)

  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol


  if (present(tol)) then
    dEq = abs(a-b) <= tol
  else
    dEq = abs(a-b) <= PREAL_EPSILON * maxval(abs([a,b]))
  endif

end function dEq


!--------------------------------------------------------------------------------------------------
!> @brief Test floating point numbers with double precision for inequality.
! replaces "!=" but for certain (relative) tolerance. Counterpart to dEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative NOT
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq(a,b,tol)

  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol


  dNeq = .not. dEq(a,b,tol)

end function dNeq


!--------------------------------------------------------------------------------------------------
!> @brief Test floating point number with double precision for equality to 0.
! replaces "==0" but everything not representable as a normal number is treated as 0. Counterpart to dNeq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq0(a,tol)

  real(pReal), intent(in)           :: a
  real(pReal), intent(in), optional :: tol


  if (present(tol)) then
    dEq0 = abs(a) <= tol
  else
    dEq0 = abs(a) <= PREAL_MIN * 10.0_pReal
  endif

end function dEq0


!--------------------------------------------------------------------------------------------------
!> @brief Test floating point number with double precision for inequality to 0.
! replaces "!=0" but everything not representable as a normal number is treated as 0. Counterpart to dEq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq0(a,tol)

  real(pReal), intent(in)           :: a
  real(pReal), intent(in), optional :: tol


  dNeq0 = .not. dEq0(a,tol)

end function dNeq0


!--------------------------------------------------------------------------------------------------
!> @brief Test complex floating point numbers with double precision for equality.
! replaces "==" but for certain (relative) tolerance. Counterpart to cNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cEq(a,b,tol)

  complex(pReal), intent(in)           :: a,b
  real(pReal),    intent(in), optional :: tol


  if (present(tol)) then
    cEq = abs(a-b) <= tol
  else
    cEq = abs(a-b) <= PREAL_EPSILON * maxval(abs([a,b]))
  endif

end function cEq


!--------------------------------------------------------------------------------------------------
!> @brief Test complex floating point numbers with double precision for inequality.
! replaces "!=" but for certain (relative) tolerance. Counterpart to cEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cNeq(a,b,tol)

  complex(pReal), intent(in)           :: a,b
  real(pReal),    intent(in), optional :: tol


  cNeq = .not. cEq(a,b,tol)

end function cNeq


!--------------------------------------------------------------------------------------------------
!> @brief Decode byte array (C_SIGNED_CHAR) as C_FLOAT array (4 byte float).
!--------------------------------------------------------------------------------------------------
pure function prec_bytesToC_FLOAT(bytes)

  integer(C_SIGNED_CHAR), dimension(:), intent(in) :: bytes                                         !< byte-wise representation of a C_FLOAT array
  real(C_FLOAT), dimension(size(bytes,kind=pI64)/(storage_size(0._C_FLOAT,pI64)/8_pI64)) :: &
    prec_bytesToC_FLOAT


  prec_bytesToC_FLOAT = transfer(bytes,prec_bytesToC_FLOAT,size(prec_bytesToC_FLOAT))

end function prec_bytesToC_FLOAT


!--------------------------------------------------------------------------------------------------
!> @brief Decode byte array (C_SIGNED_CHAR) as C_DOUBLE array (8 byte float).
!--------------------------------------------------------------------------------------------------
pure function prec_bytesToC_DOUBLE(bytes)

  integer(C_SIGNED_CHAR), dimension(:), intent(in) :: bytes                                         !< byte-wise representation of a C_DOUBLE array
  real(C_DOUBLE), dimension(size(bytes,kind=pI64)/(storage_size(0._C_DOUBLE,pI64)/8_pI64)) :: &
    prec_bytesToC_DOUBLE


  prec_bytesToC_DOUBLE = transfer(bytes,prec_bytesToC_DOUBLE,size(prec_bytesToC_DOUBLE))

end function prec_bytesToC_DOUBLE


!--------------------------------------------------------------------------------------------------
!> @brief Decode byte array (C_SIGNED_CHAR) as C_INT32_T array (4 byte signed integer).
!--------------------------------------------------------------------------------------------------
pure function prec_bytesToC_INT32_T(bytes)

  integer(C_SIGNED_CHAR), dimension(:), intent(in) :: bytes                                         !< byte-wise representation of a C_INT32_T array
  integer(C_INT32_T), dimension(size(bytes,kind=pI64)/(storage_size(0_C_INT32_T,pI64)/8_pI64)) :: &
    prec_bytesToC_INT32_T


  prec_bytesToC_INT32_T = transfer(bytes,prec_bytesToC_INT32_T,size(prec_bytesToC_INT32_T))

end function prec_bytesToC_INT32_T


!--------------------------------------------------------------------------------------------------
!> @brief Decode byte array (C_SIGNED_CHAR) as C_INT64_T array (8 byte signed integer).
!--------------------------------------------------------------------------------------------------
pure function prec_bytesToC_INT64_T(bytes)

  integer(C_SIGNED_CHAR), dimension(:), intent(in) :: bytes                                         !< byte-wise representation of a C_INT64_T array
  integer(C_INT64_T), dimension(size(bytes,kind=pI64)/(storage_size(0_C_INT64_T,pI64)/8_pI64)) :: &
     prec_bytesToC_INT64_T


  prec_bytesToC_INT64_T = transfer(bytes,prec_bytesToC_INT64_T,size(prec_bytesToC_INT64_T))

end function prec_bytesToC_INT64_T


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some prec functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer, allocatable, dimension(:) :: realloc_lhs_test
  real(pReal),   dimension(1) :: f
  integer(pInt), dimension(1) :: i
  real(pReal),   dimension(2) :: r


  realloc_lhs_test = [1,2]
  if (any(realloc_lhs_test/=[1,2]))        error stop 'LHS allocation'

  call random_number(r)
  r = r/minval(r)
  if(.not. all(dEq(r,r+PREAL_EPSILON)))    error stop 'dEq'
  if(dEq(r(1),r(2)) .and. dNeq(r(1),r(2))) error stop 'dNeq'
  if(.not. all(dEq0(r-(r+PREAL_MIN))))     error stop 'dEq0'

  ! https://www.binaryconvert.com
  ! https://www.rapidtables.com/convert/number/binary-to-decimal.html
  f = real(prec_bytesToC_FLOAT(int([-65,+11,-102,+75],C_SIGNED_CHAR)),pReal)
  if(dNeq(f(1),20191102.0_pReal,0.0_pReal)) error stop 'prec_bytesToC_FLOAT'

  f = real(prec_bytesToC_DOUBLE(int([0,0,0,-32,+119,+65,+115,65],C_SIGNED_CHAR)),pReal)
  if(dNeq(f(1),20191102.0_pReal,0.0_pReal)) error stop 'prec_bytesToC_DOUBLE'

  i = int(prec_bytesToC_INT32_T(int([+126,+23,+52,+1],C_SIGNED_CHAR)),pInt)
  if(i(1) /= 20191102_pInt)                 error stop 'prec_bytesToC_INT32_T'

  i = int(prec_bytesToC_INT64_T(int([+126,+23,+52,+1,0,0,0,0],C_SIGNED_CHAR)),pInt)
  if(i(1) /= 20191102_pInt)                 error stop 'prec_bytesToC_INT64_T'

end subroutine selfTest

end module prec
