!--------------------------------------------------------------------------------------------------
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    provides wrappers to C routines
!--------------------------------------------------------------------------------------------------
module system_routines
 
 implicit none
 private
 
 public :: &
   isDirectory, &
   getCWD, &
   getHostName, &
   setCWD

interface

 function isDirectory_C(path) bind(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   integer(C_INT) :: isDirectory_C
   character(kind=C_CHAR), dimension(1024), intent(in) :: path                                      ! C string is an array
  end function isDirectory_C

 subroutine getCurrentWorkDir_C(str, stat) bind(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   character(kind=C_CHAR), dimension(1024), intent(out)  :: str                                     ! C string is an array
   integer(C_INT),intent(out)         :: stat
  end subroutine getCurrentWorkDir_C

 subroutine getHostName_C(str, stat) bind(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   character(kind=C_CHAR), dimension(1024), intent(out)  :: str                                     ! C string is an array
   integer(C_INT),intent(out)         :: stat
  end subroutine getHostName_C

 function chdir_C(path) bind(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   integer(C_INT) :: chdir_C
   character(kind=C_CHAR), dimension(1024), intent(in) :: path                                      ! C string is an array
 end function chdir_C

end interface


contains

!--------------------------------------------------------------------------------------------------
!> @brief figures out if a given path is a directory (and not an ordinary file)
!--------------------------------------------------------------------------------------------------
logical function isDirectory(path)
  use, intrinsic :: ISO_C_Binding, only: &
    C_INT, &
    C_CHAR, &
    C_NULL_CHAR

  implicit none
  character(len=*), intent(in) :: path
  character(kind=C_CHAR), dimension(1024) :: strFixedLength
  integer :: i 

  strFixedLength = repeat(C_NULL_CHAR,len(strFixedLength))
  do i=1,len(path)                                                                                  ! copy array components
    strFixedLength(i)=path(i:i)
  enddo
  isDirectory=merge(.True.,.False.,isDirectory_C(strFixedLength) /= 0_C_INT)

end function isDirectory


!--------------------------------------------------------------------------------------------------
!> @brief gets the current working directory
!--------------------------------------------------------------------------------------------------
character(len=1024) function getCWD()
  use, intrinsic :: ISO_C_Binding, only: &
    C_INT, &
    C_CHAR, &
    C_NULL_CHAR

  implicit none
  character(kind=C_CHAR), dimension(1024) :: charArray                                              ! C string is an array
  integer(C_INT) :: stat
  integer :: i

  call getCurrentWorkDir_C(charArray,stat)
  if (stat /= 0_C_INT) then
    getCWD = 'Error occured when getting currend working directory'
  else
    getCWD = repeat('',len(getCWD))
    arrayToString: do i=1,len(getCWD)
      if (charArray(i) /= C_NULL_CHAR) then
        getCWD(i:i)=charArray(i)
      else
        exit
      endif
    enddo arrayToString
  endif

end function getCWD


!--------------------------------------------------------------------------------------------------
!> @brief gets the current host name
!--------------------------------------------------------------------------------------------------
character(len=1024) function getHostName()
  use, intrinsic :: ISO_C_Binding, only: &
    C_INT, &
    C_CHAR, &
    C_NULL_CHAR

  implicit none
  character(kind=C_CHAR), dimension(1024) :: charArray                                              ! C string is an array
  integer(C_INT) :: stat
  integer :: i

  call getHostName_C(charArray,stat)
  if (stat /= 0_C_INT) then
    getHostName = 'Error occured when getting host name'
  else
    getHostName = repeat('',len(getHostName))
    arrayToString: do i=1,len(getHostName)
      if (charArray(i) /= C_NULL_CHAR) then
        getHostName(i:i)=charArray(i)
      else
        exit
      endif
    enddo arrayToString
  endif

end function getHostName

!--------------------------------------------------------------------------------------------------
!> @brief changes the current working directory
!--------------------------------------------------------------------------------------------------
logical function setCWD(path)
  use, intrinsic :: ISO_C_Binding, only: &
    C_INT, &
    C_CHAR, &
    C_NULL_CHAR

  implicit none
  character(len=*), intent(in) :: path
  character(kind=C_CHAR), dimension(1024) :: strFixedLength                                         ! C string is an array
  integer :: i

  strFixedLength = repeat(C_NULL_CHAR,len(strFixedLength))
  do i=1,len(path)                                                                                  ! copy array components
    strFixedLength(i)=path(i:i)
  enddo
  setCWD=merge(.True.,.False.,chdir_C(strFixedLength) /= 0_C_INT)

end function setCWD

end module system_routines

