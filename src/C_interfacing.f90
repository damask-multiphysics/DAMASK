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

#if (defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) \
 || (defined(__GFORTRAN__) && __GNUC__ < 15) \
 || defined(__flang__)
  public :: &
    f_c_string
#endif

contains

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
