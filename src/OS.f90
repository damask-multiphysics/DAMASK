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

  interface

    function set_CWD_C(cwd) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_INT, C_CHAR

      implicit none(type,external)
      integer(C_INT) :: set_CWD_C
      character(kind=C_CHAR), dimension(*), intent(in) :: cwd
    end function set_CWD_C

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


  call get_CWD_C(OS_getCWD,stat)
  if (stat /= 0) error stop 'invalid working directory'

end function OS_getCWD


!--------------------------------------------------------------------------------------------------
!> @brief Get the host name.
!--------------------------------------------------------------------------------------------------
function OS_getHostName()

  character(kind=C_CHAR,len=:), allocatable :: OS_getHostName

  integer(C_INT) :: stat


  call get_hostname_C(OS_getHostName,stat)
  if (stat /= 0) OS_getHostName = 'n/a (Error!)'

end function OS_getHostName


!--------------------------------------------------------------------------------------------------
!> @brief Get the user name.
!--------------------------------------------------------------------------------------------------
function OS_getUserName()

  character(kind=C_CHAR,len=:), allocatable :: OS_getUserName

  integer(C_INT) :: stat


  call get_username_C(OS_getUserName,stat)
  if (stat /= 0) OS_getUserName = 'n/a (Error!)'

end function OS_getUserName


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some OS functions.
!--------------------------------------------------------------------------------------------------
subroutine OS_selfTest()

  character(len=:), allocatable :: cwd, hostname, username
  logical :: failed


  cwd = OS_getCWD()
  if (len_trim(cwd) == 0) error stop 'OS_getCWD returned empty string'
  hostname = OS_getHostName()
  if (len_trim(hostname) == 0) error stop 'OS_getHostName returned empty string'
  username = OS_getUserName()
  if (len_trim(username) == 0) error stop 'OS_getUserName returned empty string'
  failed = OS_setCWD(cwd)
  if (failed) error stop 'OS_setCWD failed to set to current directory'

end subroutine OS_selfTest


end module OS
