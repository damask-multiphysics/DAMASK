!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @author Philip Eisenlohr, Michigan State University
!> @brief Miscellaneous tools.
!--------------------------------------------------------------------------------------------------
module misc
  use prec

  implicit none(type,external)
  private

  interface misc_optional
    module procedure misc_optional_bool
    module procedure misc_optional_integer
    module procedure misc_optional_real
    module procedure misc_optional_str
  end interface misc_optional

  public :: &
    misc_init, &
    misc_optional

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine misc_init()

  print'(/,1x,a)', '<<<+-  misc init  -+>>>'

  call misc_selfTest()

end subroutine misc_init


!--------------------------------------------------------------------------------------------------
!> @brief Return bool value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure function misc_optional_bool(given,default) result(var)

  logical, intent(in), optional :: given
  logical, intent(in)           :: default
  logical                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_bool


!--------------------------------------------------------------------------------------------------
!> @brief Return integer value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure function misc_optional_integer(given,default) result(var)

  integer, intent(in), optional :: given
  integer, intent(in)           :: default
  integer                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_integer


!--------------------------------------------------------------------------------------------------
!> @brief Return real value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure function misc_optional_real(given,default) result(var)

  real(pReal), intent(in), optional :: given
  real(pReal), intent(in)           :: default
  real(pReal)                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_real


!--------------------------------------------------------------------------------------------------
!> @brief Return string value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure function misc_optional_str(given,default) result(var)

  character(len=*), intent(in), optional :: given
  character(len=*), intent(in)           :: default
  character(len=:), allocatable          :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_str


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some misc functions.
!--------------------------------------------------------------------------------------------------
subroutine misc_selfTest()

  real(pReal) :: r

  call random_number(r)
  if (test_str('DAMASK') /= 'DAMASK')                        error stop 'optional_str, present'
  if (test_str() /= 'default')                               error stop 'optional_str, not present'
  if (misc_optional(default='default') /= 'default')         error stop 'optional_str, default only'
  if (test_int(20191102) /= 20191102)                        error stop 'optional_int, present'
  if (test_int() /= 42)                                      error stop 'optional_int, not present'
  if (misc_optional(default=20191102) /= 20191102)           error stop 'optional_int, default only'
  if (dNeq(test_real(r),r))                                  error stop 'optional_real, present'
  if (dNeq(test_real(),0.0_pReal))                           error stop 'optional_real, not present'
  if (dNeq(misc_optional(default=r),r))                      error stop 'optional_real, default only'
  if (test_bool(r<0.5_pReal) .neqv. r<0.5_pReal)             error stop 'optional_bool, present'
  if (.not. test_bool())                                     error stop 'optional_bool, not present'
  if (misc_optional(default=r>0.5_pReal) .neqv. r>0.5_pReal) error stop 'optional_bool, default only'

contains

  function test_str(str_in) result(str_out)

    character(len=:), allocatable           :: str_out
    character(len=*), intent(in), optional :: str_in


    str_out = misc_optional_str(str_in,'default')

  end function test_str


  function test_int(int_in) result(int_out)

    integer                       :: int_out
    integer, intent(in), optional :: int_in


    int_out = misc_optional_integer(int_in,42)

  end function test_int


  function test_real(real_in) result(real_out)

    real(pReal)                       :: real_out
    real(pReal), intent(in), optional :: real_in


    real_out = misc_optional_real(real_in,0.0_pReal)

  end function test_real


  function test_bool(bool_in) result(bool_out)

    logical                       :: bool_out
    logical, intent(in), optional :: bool_in


    bool_out = misc_optional_bool(bool_in,.true.)

  end function test_bool


end subroutine misc_selfTest

end module misc
