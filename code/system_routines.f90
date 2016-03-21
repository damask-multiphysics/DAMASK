module system_routines
 
 implicit none
 private
 
 public :: &
   isDirectory, &
   getCWD2

interface

 function  isDirectory_C(path) BIND(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   integer(C_INT) :: isDirectory_C
   character(kind=C_CHAR),intent(in) :: path(*)
  end function isDirectory_C

 subroutine getCurrentWorkDir_C(str_out, stat) bind(C)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR
   character( kind=c_char ), dimension(*), intent(inout)  :: str_out
   integer(C_INT),intent(out)         :: stat

  end subroutine getCurrentWorkDir_C

end interface


contains

 logical function isDirectory(path)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT

   implicit none
   character(len=*), intent(in) :: path

   isDirectory=merge(.True.,.False.,isDirectory_C(trim(path)) /= 0_C_INT)

 end function isDirectory


 logical function getCWD2(str)
   use, intrinsic :: ISO_C_Binding, only: &
     C_INT, &
     C_CHAR, &
     C_NULL_CHAR

   implicit none
   character(len=*), intent(out) :: str
   character(len=1024) :: strFixedLength
   integer(C_INT) :: stat
 
   str = repeat(C_NULL_CHAR,1024)
   call getCurrentWorkDir_C(strFixedLength,stat)
   str = strFixedLength(1:scan(strFixedLength,C_NULL_CHAR,.True.)-1)
   getCWD2=merge(.True.,.False.,stat /= 0_C_INT)

 end function getCWD2

end module system_routines

