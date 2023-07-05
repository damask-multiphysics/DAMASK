!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Wrappers to C routines for system operations
!--------------------------------------------------------------------------------------------------
module system_routines
  use, intrinsic :: ISO_C_Binding

  use prec

  implicit none(type,external)
  private

  public :: &
    setCWD, &
    getCWD, &
    getHostName, &
    getUserName, &
    signalint_C, &
    signalusr1_C, &
    signalusr2_C, &
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

  end interface

contains


!--------------------------------------------------------------------------------------------------
!> @brief Set the current working directory.
!--------------------------------------------------------------------------------------------------
logical function setCWD(path)

  character(len=*), intent(in) :: path


  setCWD = setCWD_C(f_c_string(path)) /= 0_C_INT

end function setCWD


!--------------------------------------------------------------------------------------------------
!> @brief Get the current working directory.
!--------------------------------------------------------------------------------------------------
function getCWD()

  character(len=:), allocatable :: getCWD

  character(kind=C_CHAR), dimension(pPathLen+1) :: getCWD_Cstring
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

  character(kind=C_CHAR), dimension(pSTRLEN+1) :: getHostName_Cstring
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

  character(kind=C_CHAR), dimension(pSTRLEN+1) :: getUserName_Cstring
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

  character(kind=C_CHAR), dimension(:), intent(in) :: c_string
  character(len=:),       allocatable              :: f_string

  integer(pI64) :: i


  allocate(character(len=size(c_string,kind=pI64))::f_string)
  arrayToString: do i=1_pI64,len(f_string,pI64)
    if (c_string(i) /= C_NULL_CHAR) then
      f_string(i:i)=c_string(i)
    else
      f_string = f_string(:i-1_pI64)
      exit
    end if
  end do arrayToString

end function c_f_string


!--------------------------------------------------------------------------------------------------
!> @brief Convert Fortran string to C string.
!> @details: C string is NULL terminated and, hence, longer by one than the Fortran string.
!--------------------------------------------------------------------------------------------------
pure function f_c_string(f_string) result(c_string)

  character(len=*), intent(in) :: f_string
  character(kind=C_CHAR), dimension(len_trim(f_string,pI64)+1_pI64) :: c_string


  c_string = transfer(trim(f_string)//C_NULL_CHAR,c_string,size=size(c_string,kind=pI64))

end function f_c_string


end module system_routines

