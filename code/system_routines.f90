!--------------------------------------------------------------------------------------------------
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    provides wrappers to C routines
!--------------------------------------------------------------------------------------------------
module system_routines
 
 implicit none
 private
 
 public :: &
   isDirectory, &
   getCWD

interface

 function  isDirectory_C(path) BIND(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   integer(C_INT) :: isDirectory_C
   character(kind=C_CHAR),intent(in) :: path(*)
  end function isDirectory_C

 subroutine getCurrentWorkDir_C(str, stat) bind(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   character(kind=c_char), dimension(1024), intent(out)  :: str                                     ! C string is an array
   integer(C_INT),intent(out)         :: stat

  end subroutine getCurrentWorkDir_C

end interface


contains

!--------------------------------------------------------------------------------------------------
!> @brief figures out if a given path is a directory (and not an ordinary file)
!--------------------------------------------------------------------------------------------------
logical function isDirectory(path)
  use, intrinsic :: ISO_C_Binding, only: &
    C_INT

  implicit none
  character(len=*), intent(in) :: path

  isDirectory=merge(.True.,.False.,isDirectory_C(trim(path)) /= 0_C_INT)

end function isDirectory


!--------------------------------------------------------------------------------------------------
!> @brief gets the current working directory
!--------------------------------------------------------------------------------------------------
logical function getCWD(str)
  use, intrinsic :: ISO_C_Binding, only: &
    C_INT, &
    C_CHAR, &
    C_NULL_CHAR
  implicit none
  character(len=*), intent(out) :: str
  character(kind=c_char), dimension(1024) :: strFixedLength                                         ! C string is an array
  integer(C_INT) :: stat
  integer :: i

  str = repeat('',len(str))
  call getCurrentWorkDir_C(strFixedLength,stat)
  do i=1,1024                                                                                       ! copy array components until Null string is found
    if (strFixedLength(i) /= C_NULL_CHAR) then
      str(i:i)=strFixedLength(i)
    else
      exit
    endif
  enddo
  getCWD=merge(.True.,.False.,stat /= 0_C_INT)

end function getCWD

end module system_routines

