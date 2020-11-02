!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Wrappers to C routines for system operations
!--------------------------------------------------------------------------------------------------
module system_routines
  use, intrinsic :: ISO_C_Binding

  use prec

  implicit none

  public :: &
    signalterm_C, &
    signalusr1_C, &
    signalusr2_C, &
    isDirectory, &
    getCWD, &
    getHostName, &
    setCWD

  interface

  function isDirectory_C(path) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_INT, &
      C_CHAR

    use prec

    integer(C_INT) :: isDirectory_C
    character(kind=C_CHAR), dimension(pPathLen), intent(in) :: path                                 ! C string is an array
  end function isDirectory_C

  subroutine getCurrentWorkDir_C(path, stat) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_INT, &
      C_CHAR

    use prec

    character(kind=C_CHAR), dimension(pPathLen), intent(out) :: path                                ! C string is an array
    integer(C_INT),                              intent(out) :: stat
  end subroutine getCurrentWorkDir_C

  subroutine getHostName_C(str, stat) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_INT, &
      C_CHAR

    use prec

    character(kind=C_CHAR), dimension(pStringLen), intent(out) :: str                               ! C string is an array
    integer(C_INT),                                intent(out) :: stat
  end subroutine getHostName_C

  function chdir_C(path) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_INT, &
      C_CHAR

    use prec

    integer(C_INT) :: chdir_C
    character(kind=C_CHAR), dimension(pPathLen), intent(in) :: path                                 ! C string is an array
  end function chdir_C

  subroutine signalterm_C(handler) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_FUNPTR

    type(C_FUNPTR), intent(in), value :: handler
  end subroutine signalterm_C

  subroutine signalusr1_C(handler) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_FUNPTR

    type(C_FUNPTR), intent(in), value :: handler
  end subroutine signalusr1_C

  subroutine signalusr2_C(handler) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_FUNPTR

    type(C_FUNPTR), intent(in), value :: handler
  end subroutine signalusr2_C

  end interface

contains

!--------------------------------------------------------------------------------------------------
!> @brief figures out if a given path is a directory (and not an ordinary file)
!--------------------------------------------------------------------------------------------------
logical function isDirectory(path)

  character(len=*), intent(in) :: path
  
  isDirectory=merge(.True.,.False.,isDirectory_C(f_c_string(path)) /= 0_C_INT)

end function isDirectory


!--------------------------------------------------------------------------------------------------
!> @brief gets the current working directory
!--------------------------------------------------------------------------------------------------
function getCWD()

  character(kind=C_CHAR), dimension(pPathLen) :: getCWD_Cstring
  character(len=:),       allocatable         :: getCWD
  integer(C_INT) :: stat

  call getCurrentWorkDir_C(getCWD_Cstring,stat)

  if(stat == 0) then
    getCWD = c_f_string(getCWD_Cstring)
  else
    getCWD = 'Error occured when getting currend working directory'
  endif

end function getCWD


!--------------------------------------------------------------------------------------------------
!> @brief gets the current host name
!--------------------------------------------------------------------------------------------------
function getHostName()

  character(kind=C_CHAR), dimension(pPathLen) :: getHostName_Cstring
  character(len=:),       allocatable         :: getHostName
  integer(C_INT) :: stat

  call getHostName_C(getHostName_Cstring,stat)

  if(stat == 0) then
    getHostName = c_f_string(getHostName_Cstring)
  else
    getHostName = 'Error occured when getting host name'
  endif

end function getHostName


!--------------------------------------------------------------------------------------------------
!> @brief changes the current working directory
!--------------------------------------------------------------------------------------------------
logical function setCWD(path)

  character(len=*), intent(in) :: path

  setCWD=merge(.True.,.False.,chdir_C(f_c_string(path)) /= 0_C_INT)

end function setCWD


!--------------------------------------------------------------------------------------------------
!> @brief convert C string to Fortran string
!> @details: C string is NULL terminated and, hence, longer by one than the Fortran string
!--------------------------------------------------------------------------------------------------
pure function c_f_string(c_string) result(f_string)

  character(kind=C_CHAR), dimension(:), intent(in) :: c_string
  character(len=:),       allocatable              :: f_string
  integer :: i

  allocate(character(len=size(c_string))::f_string)
  arrayToString: do i=1,len(f_string)
    if (c_string(i) /= C_NULL_CHAR) then
      f_string(i:i)=c_string(i)
    else
      f_string = f_string(:i-1)
      exit
    endif
  enddo arrayToString

end function c_f_string


!--------------------------------------------------------------------------------------------------
!> @brief convert Fortran string to C string
!> @details: C string is NULL terminated and, hence, longer by one than the Fortran string
!--------------------------------------------------------------------------------------------------
pure function f_c_string(f_string) result(c_string)

  character(len=*), intent(in)                       :: f_string
  character(kind=C_CHAR), dimension(len(f_string)+1) :: c_string
  integer :: i

  do i=1,len(f_string)
    c_string(i)=f_string(i:i)
  enddo
  c_string(i) = C_NULL_CHAR

end function f_c_string


end module system_routines

