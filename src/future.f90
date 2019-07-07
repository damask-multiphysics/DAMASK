!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief New fortran functions for compiler versions that do not support them
!--------------------------------------------------------------------------------------------------
module future
  use prec

  implicit none
  public

contains

#if defined(__GFORTRAN__) || __INTEL_COMPILER < 1800
!--------------------------------------------------------------------------------------------------
!> @brief substitute for the findloc intrinsic (only for integer, dimension(:) at the moment)
!--------------------------------------------------------------------------------------------------
function findloc(a,v)

  integer, intent(in), dimension(:) :: a
  integer, intent(in) :: v
  integer :: i,j
  integer, allocatable, dimension(:) ::  findloc

  allocate(findloc(count(a==v)))
  j = 1
  do i = 1, size(a)
    if (a(i)==v) then
      findloc(j) = i
      j = j + 1
    endif
  enddo
end function findloc
#endif


#if defined(__PGI)
!--------------------------------------------------------------------------------------------------
!> @brief substitute for the norm2 intrinsic (only for real, dimension(3) at the moment)
!--------------------------------------------------------------------------------------------------
real(pReal) pure function norm2(v)
  
  real(pReal), intent(in), dimension(3) :: v
  
  norm2 = sqrt(sum(v**2))
 
end function norm2
#endif

end module future
