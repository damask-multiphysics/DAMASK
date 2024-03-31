!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Wrappers to C routines for system operations
!--------------------------------------------------------------------------------------------------
module system_routines
  use, intrinsic :: ISO_C_Binding
  use, intrinsic :: ISO_fortran_env

  use prec

  implicit none(type,external)
  private

  public :: &
    system_routines_init, &
    system_routines_selfTest, &
    setCWD, &
    getCWD, &
    getHostName, &
    getUserName, &
    signalint_C, &
    signalusr1_C, &
    signalusr2_C, &
    isatty, &
    f_c_string, &
    free_C


  interface


    function setCWD_C(cwd) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT, C_CHAR

      implicit none(type,external)
      integer(C_INT) :: setCWD_C
      character(kind=C_CHAR), dimension(*), intent(in) :: cwd
    end function setCWD_C

    subroutine getCWD_C(cwd, stat) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT, C_CHAR
      use prec

      implicit none(type,external)
      character(kind=C_CHAR), dimension(pPathLen+1), intent(out) :: cwd                             ! NULL-terminated array
      integer(C_INT),                                intent(out) :: stat
    end subroutine getCWD_C

    subroutine getHostName_C(hostname, stat) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT, C_CHAR
      use prec

      implicit none(type,external)
      character(kind=C_CHAR), dimension(pSTRLEN+1), intent(out) :: hostname                         ! NULL-terminated array
      integer(C_INT),                               intent(out) :: stat
    end subroutine getHostName_C

    subroutine getUserName_C(username, stat) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT, C_CHAR
      use prec

      implicit none(type,external)
      character(kind=C_CHAR), dimension(pSTRLEN+1), intent(out) :: username                         ! NULL-terminated array
      integer(C_INT),                               intent(out) :: stat
    end subroutine getUserName_C

    subroutine signalint_C(handler) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_FUNPTR

      implicit none(type,external)
      type(C_FUNPTR), intent(in), value :: handler
    end subroutine signalint_C

    subroutine signalusr1_C(handler) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_FUNPTR

      implicit none(type,external)
      type(C_FUNPTR), intent(in), value :: handler
    end subroutine signalusr1_C

    subroutine signalusr2_C(handler) bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_FUNPTR

      implicit none(type,external)
      type(C_FUNPTR), intent(in), value :: handler
    end subroutine signalusr2_C

    subroutine free_C(ptr) bind(C,name='free')
      use, intrinsic :: ISO_C_Binding, only: C_PTR

      implicit none(type,external)
      type(C_PTR), value :: ptr
    end subroutine free_C

    function stdout_isatty_C() bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT

      implicit none(type,external)
      integer(C_INT) :: stdout_isatty_C
    end function stdout_isatty_C

    function stderr_isatty_C() bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT

      implicit none(type,external)
      integer(C_INT) :: stderr_isatty_C
    end function stderr_isatty_C

    function stdin_isatty_C() bind(C)
      use, intrinsic :: ISO_C_Binding, only: C_INT

      implicit none(type,external)
      integer(C_INT) :: stdin_isatty_C
    end function stdin_isatty_C


  end interface

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine system_routines_init()

  print'(/,1x,a)', '<<<+-  system_routines init  -+>>>'; flush(OUTPUT_UNIT)

  call system_routines_selfTest()

end subroutine system_routines_init


!--------------------------------------------------------------------------------------------------
!> @brief Set the current working directory.
!--------------------------------------------------------------------------------------------------
logical function setCWD(path)

  character(len=*), intent(in) :: path


  setCWD = setCWD_C(f_c_string(path)) /= 0_C_INT

  call system_routines_selfTest()

end function setCWD


!--------------------------------------------------------------------------------------------------
!> @brief Get the current working directory.
!--------------------------------------------------------------------------------------------------
function getCWD()

  character(len=:), allocatable :: getCWD

  character(kind=C_CHAR,len=(pPathLen+1)) :: getCWD_Cstring
  integer(C_INT) :: stat


  call getCWD_C(getCWD_Cstring,stat)

  if (stat == 0) then
    getCWD = c_f_string(getCWD_Cstring)
  else
    error stop 'invalid working directory'
  end if

end function getCWD


!--------------------------------------------------------------------------------------------------
!> @brief Get the host name.
!--------------------------------------------------------------------------------------------------
function getHostName()

  character(len=:), allocatable :: getHostName

  character(kind=C_CHAR,len=(pSTRLEN+1)) :: getHostName_Cstring
  integer(C_INT) :: stat


  call getHostName_C(getHostName_Cstring,stat)

  if (stat == 0) then
    getHostName = c_f_string(getHostName_Cstring)
  else
    getHostName = 'n/a (Error!)'
  end if

end function getHostName


!--------------------------------------------------------------------------------------------------
!> @brief Get the user name.
!--------------------------------------------------------------------------------------------------
function getUserName()

  character(len=:), allocatable :: getUserName

  character(kind=C_CHAR,len=(pSTRLEN+1)) :: getUserName_Cstring
  integer(C_INT) :: stat


  call getUserName_C(getUserName_Cstring,stat)

  if (stat == 0) then
    getUserName = c_f_string(getUserName_Cstring)
  else
    getUserName = 'n/a (Error!)'
  end if

end function getUserName

!--------------------------------------------------------------------------------------------------
!> @brief Convert C string to Fortran string.
!> @details: C string is NULL terminated and, hence, longer by one than the Fortran string.
!--------------------------------------------------------------------------------------------------
pure function c_f_string(c_string) result(f_string)

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

end function c_f_string


#if __INTEL_LLVM_COMPILER < 20240100
!--------------------------------------------------------------------------------------------------
!> @brief Convert Fortran string to C string.
!> @details: C string is NULL terminated and, hence, longer by one than the Fortran string.
!--------------------------------------------------------------------------------------------------
pure function f_c_string(f_string) result(c_string)

  character(len=*), intent(in) :: f_string
  character(kind=C_CHAR,len=len_trim(f_string,pI64)+1_pI64) :: c_string


  c_string = trim(f_string)//C_NULL_CHAR

end function f_c_string
#endif


!--------------------------------------------------------------------------------------------------
!> @brief Test whether a file descriptor refers to a terminal.
!> @detail A terminal is neither a file nor a redirected STDOUT/STDERR/STDIN.
!--------------------------------------------------------------------------------------------------
logical function isatty(unit)

  integer, intent(in) :: unit


  select case(unit)
#ifndef LOGFILE
    case (OUTPUT_UNIT)
      isatty = stdout_isatty_C()==1
    case (ERROR_UNIT)
      isatty = stderr_isatty_C()==1
#endif
    case (INPUT_UNIT)
      isatty = stdin_isatty_C()==1
    case default
      isatty = .false.
  end select

end function isatty


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some system_routine functions.
!--------------------------------------------------------------------------------------------------
subroutine system_routines_selfTest()

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

end subroutine system_routines_selfTest


end module system_routines

