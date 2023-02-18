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
    module procedure misc_optional_string
  end interface misc_optional

  public :: &
    misc_optional

contains


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
pure function misc_optional_string(given,default) result(var)

  character(len=*), intent(in), optional :: given
  character(len=*), intent(in)           :: default
  character(len=:), allocatable          :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_string

end module misc
