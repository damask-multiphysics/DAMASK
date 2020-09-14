!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Inflate zlib compressed data
!--------------------------------------------------------------------------------------------------
module zlib
  use prec

  implicit none
  private

  public :: &
    zlib_inflate

  interface

  subroutine inflate_C(s_deflated,s_inflated,deflated,inflated) bind(C)
    use, intrinsic :: ISO_C_Binding, only: &
      C_SIGNED_CHAR, C_INT64_T

    integer(C_INT64_T),                            intent(in)  :: s_deflated,s_inflated
    integer(C_SIGNED_CHAR), dimension(s_deflated), intent(in)  :: deflated
    integer(C_SIGNED_CHAR), dimension(s_inflated), intent(out) :: inflated

  end subroutine inflate_C

  end interface

contains

!--------------------------------------------------------------------------------------------------
!> @brief Inflate byte-wise representation
!--------------------------------------------------------------------------------------------------
function zlib_inflate(deflated,size_inflated)

  integer(C_SIGNED_CHAR), dimension(:), intent(in) :: deflated
  integer(pI64),                        intent(in) :: size_inflated

  integer(C_SIGNED_CHAR), dimension(size_inflated) :: zlib_inflate

  call inflate_C(size(deflated,kind=C_INT64_T),int(size_inflated,C_INT64_T),deflated,zlib_inflate)

end function zlib_inflate

end module zlib
