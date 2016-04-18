!********************************************************************
! quit subroutine to satisfy IO_error for core module
!
!********************************************************************
subroutine quit(stop_id)
 use prec, only: &
   pInt
   
 implicit none
 integer(pInt), intent(in) :: stop_id

end subroutine
