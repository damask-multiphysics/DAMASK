! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Wrappers to C routines for system operations
!--------------------------------------------------------------------------------------------------
module OS
  use, intrinsic :: ISO_C_binding
  use, intrinsic :: ISO_fortran_env

  use prec
  use C_interfacing

  implicit none(type,external)
  private

  public :: &
    OS_init, &
    OS_selfTest, &
    OS_setCWD, &
    OS_getCWD, &
    OS_getHostName, &
    OS_getUserName
#ifdef OLD_STYLE_C_TO_FORTRAN_STRING
  public :: free_C
#endif

  interface

    function set_CWD_C(cwd) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR

      implicit none(type,external)
      integer(C_INT) :: set_CWD_C
      character(kind=C_CHAR), dimension(*), intent(in) :: cwd
    end function set_CWD_C
#ifndef OLD_STYLE_C_TO_FORTRAN_STRING
    subroutine get_CWD_C(cwd, stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR

      implicit none(type,external)
      character(kind=C_CHAR,len=:), allocatable, intent(out) :: cwd
      integer(C_INT),                            intent(out) :: stat
    end subroutine get_CWD_C

    subroutine get_hostname_C(hostname, stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR

      implicit none(type,external)
      character(kind=C_CHAR,len=:), allocatable, intent(out) :: hostname
      integer(C_INT),                            intent(out) :: stat
    end subroutine get_hostname_C

    subroutine get_username_C(username, stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR

      implicit none(type,external)
      character(kind=C_CHAR,len=:), allocatable, intent(out) :: username
      integer(C_INT),                            intent(out) :: stat
    end subroutine get_username_C
#else
    subroutine get_CWD_C(cwd, stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR
      use prec

      implicit none(type,external)
      character(kind=C_CHAR), dimension(4096+1), intent(out) :: cwd                                 ! NULL-terminated array
      integer(C_INT),                            intent(out) :: stat
    end subroutine get_CWD_C

    subroutine get_hostname_C(hostname, stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR
      use prec

      implicit none(type,external)
      character(kind=C_CHAR), dimension(pSTRLEN+1), intent(out) :: hostname                         ! NULL-terminated array
      integer(C_INT),                               intent(out) :: stat
    end subroutine get_hostname_C

    subroutine get_username_C(username, stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR
      use prec

      implicit none(type,external)
      character(kind=C_CHAR), dimension(pSTRLEN+1), intent(out) :: username                         ! NULL-terminated array
      integer(C_INT),                               intent(out) :: stat
    end subroutine get_username_C

    subroutine free_C(ptr) bind(C,name='free')
      use, intrinsic :: ISO_C_binding, only: C_PTR

      implicit none(type,external)
      type(C_PTR), value :: ptr
    end subroutine free_C
#endif

  end interface

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine OS_init()

  print'(/,1x,a)', '<<<+-  OS init  -+>>>'; flush(OUTPUT_UNIT)

  call OS_selfTest()

end subroutine OS_init


!--------------------------------------------------------------------------------------------------
!> @brief Set the current working directory.
!--------------------------------------------------------------------------------------------------
logical function OS_setCWD(path)

  character(len=*), intent(in) :: path


  OS_setCWD = set_CWD_C(f_c_string(path)) /= 0_C_INT

end function OS_setCWD


!--------------------------------------------------------------------------------------------------
!> @brief Get the current working directory.
!--------------------------------------------------------------------------------------------------
function OS_getCWD()

  character(kind=C_CHAR,len=:), allocatable :: OS_getCWD

  integer(C_INT) :: stat

#ifndef OLD_STYLE_C_TO_FORTRAN_STRING
  call get_CWD_C(OS_getCWD,stat)
  if (stat /= 0) error stop 'invalid working directory'
#else
  character(kind=C_CHAR,len=(4096+1)) :: CWD_Cstring


  call get_CWD_C(CWD_Cstring,stat)

  if (stat == 0) then
    OS_getCWD = c_f_string(CWD_Cstring)
  else
    error stop 'invalid working directory'
  end if
#endif
end function OS_getCWD


!--------------------------------------------------------------------------------------------------
!> @brief Get the host name.
!--------------------------------------------------------------------------------------------------
function OS_getHostName()

  character(kind=C_CHAR,len=:), allocatable :: OS_getHostName

  integer(C_INT) :: stat

#ifndef OLD_STYLE_C_TO_FORTRAN_STRING
  call get_hostname_C(OS_getHostName,stat)
  if (stat /= 0) OS_getHostName = 'n/a (Error!)'
#else
  character(kind=C_CHAR,len=(pSTRLEN+1)) :: hostname_Cstring


  call get_hostname_C(hostname_Cstring,stat)

  if (stat == 0) then
    OS_getHostName = c_f_string(hostname_Cstring)
  else
    OS_getHostName = 'n/a (Error!)'
  end if
#endif
end function OS_getHostName


!--------------------------------------------------------------------------------------------------
!> @brief Get the user name.
!--------------------------------------------------------------------------------------------------
function OS_getUserName()

  character(kind=C_CHAR,len=:), allocatable :: OS_getUserName

  integer(C_INT) :: stat

#ifndef OLD_STYLE_C_TO_FORTRAN_STRING
  call get_username_C(OS_getUserName,stat)
  if (stat /= 0) OS_getUserName = 'n/a (Error!)'
#else
  character(kind=C_CHAR,len=(pSTRLEN+1)) :: username_Cstring


  call get_username_C(username_Cstring,stat)

  if (stat == 0) then
    OS_getUserName = c_f_string(username_Cstring)
  else
    OS_getUserName = 'n/a (Error!)'
  end if
#endif
end function OS_getUserName


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some OS functions.
!--------------------------------------------------------------------------------------------------
subroutine OS_selfTest()

  real :: r
  real, dimension(:), allocatable :: rnd_real
  character(len=:), allocatable :: rnd_str
  integer :: i


  call random_number(r)
  allocate(rnd_real(30+int(r*50.)))
  call random_number(rnd_real)
  allocate(character(size(rnd_real))::rnd_str)

  do i = 1, size(rnd_real)
    rnd_str(i:i) = char(32 + int(rnd_real(i)*(127.-32.)))
  end do

  if (c_f_string(f_c_string(rnd_str)) /= rnd_str) error stop 'c_f_string/f_c_string'

end subroutine OS_selfTest


end module OS
