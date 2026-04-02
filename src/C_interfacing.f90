! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
!> @author Martin Diehl, KU Leuven
!> @brief  Helpers to interface between C and Fortran strings
!--------------------------------------------------------------------------------------------------
module C_interfacing
  use, intrinsic :: ISO_C_binding

  use prec

  implicit none(type,external)
  private

  public :: &
    c_f_string
#if (defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) \
 || (defined(__GFORTRAN__) && __GNUC__ < 15) \
 || defined(__flang__)
  public :: &
    f_c_string
#endif

  interface c_f_string
    module procedure c_f_string_scalar
    module procedure c_f_string_array
  end interface

contains

!--------------------------------------------------------------------------------------------------
!> @brief Convert C string to Fortran string.
!> @details: C string is a fixed-size fortran array referencing a NULL terminated c string.
!            Due to the NULL-termination, C string has one more element than the Fortran string.
!--------------------------------------------------------------------------------------------------
pure function c_f_string_scalar(c_string) result(f_string)

  character(kind=C_CHAR,len=*), intent(in) :: c_string
  character(len=:), allocatable            :: f_string

  integer(pI64) :: i


  allocate(character(len=len(c_string,kind=pI64))::f_string)
  do i=1_pI64,len(f_string,pI64)
    if (c_string(i:i) /= C_NULL_CHAR) then
      f_string(i:i)=c_string(i:i)
    else
      f_string = f_string(:i-1_pI64)
      exit
    end if
  end do

end function c_f_string_scalar

!--------------------------------------------------------------------------------------------------
!> @brief Convert C string pointer to Fortran string.
!> @details: C string is an assumed-size fortran array referencing a NULL terminated c string.
!--------------------------------------------------------------------------------------------------
pure function c_f_string_array(c_string) result(f_string)

  character(kind=C_CHAR), intent(in), dimension(*) :: c_string
  character(len=:), allocatable                    :: f_string

  integer :: n, i


  n = 0
  do
    if (c_string(n+1) == C_NULL_CHAR) exit
    n = n + 1
  end do

  allocate(character(len=n) :: f_string)
  do i = 1, n
    f_string(i:i) = c_string(i)
  end do

end function c_f_string_array


#if (defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) \
 || (defined(__GFORTRAN__) && __GNUC__ < 15) \
 || defined(__flang__)
!--------------------------------------------------------------------------------------------------
!> @brief Fortran 2023 "f_c_string" (without optional argument).
!> @details: C string is NULL terminated and, hence, longer by one than the Fortran string.
!--------------------------------------------------------------------------------------------------
pure function f_c_string(f_string) result(c_string)

  character(len=*), intent(in) :: f_string
  character(kind=C_CHAR,len=len_trim(f_string,pI64)+1_pI64) :: c_string


  c_string = trim(f_string)//C_NULL_CHAR

end function f_c_string
#endif


end module C_interfacing
