! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Inflate zlib compressed data
!--------------------------------------------------------------------------------------------------
module zlib
  use prec
  use IO
  use base64

  implicit none(type,external)
  private

  public :: &
    zlib_init, &
    zlib_inflate

  interface

    subroutine inflate_C(s_deflated,s_inflated,deflated,inflated,stat) bind(C)
      use, intrinsic :: ISO_C_binding, only: C_SIGNED_CHAR, C_INT64_T, C_INT
      implicit none(type,external)

      integer(C_INT64_T),                            intent(in)  :: s_deflated,s_inflated
      integer(C_SIGNED_CHAR), dimension(s_deflated), intent(in)  :: deflated                        ! ok for unsigned char
      integer(C_SIGNED_CHAR), dimension(s_inflated), intent(out) :: inflated                        ! ok for unsigned char
      integer(C_INT),                                intent(out) :: stat
    end subroutine inflate_C

  end interface

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine zlib_init()

  print'(/,1x,a)', '<<<+-  zlib init  -+>>>'; flush(IO_STDOUT)

  call selfTest()

end subroutine zlib_init


!--------------------------------------------------------------------------------------------------
!> @brief Inflate byte-wise representation.
!--------------------------------------------------------------------------------------------------
function zlib_inflate(deflated,size_inflated)

  integer(C_SIGNED_CHAR), dimension(:), intent(in) :: deflated
  integer(pI64),                        intent(in) :: size_inflated
  integer(C_SIGNED_CHAR), dimension(size_inflated) :: zlib_inflate

  integer(C_INT) :: stat


  call inflate_C(size(deflated,kind=C_INT64_T),int(size_inflated,C_INT64_T),deflated,zlib_inflate,stat)
  if (stat /= 0) error stop 'inflate failed'

end function zlib_inflate


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of zlib inflate.
!> @details From Python: print(base64.b64encode(zlib.compress(b'DAMASK zlib FFTW PETSc fyaml',9))).
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  character(len=*), parameter :: deflated = 'eNpzcfR1DPZWqMrJTFJwcwsJVwhwDQlOVkirTMzNAQB5gwjS'
  character(len=*), parameter :: inflated = 'DAMASK zlib FFTW PETSc fyaml'
  integer(C_SIGNED_CHAR), dimension(:), allocatable :: bytes_deflated
  integer(C_SIGNED_CHAR), dimension(len(inflated)) :: bytes_inflated


  bytes_deflated = base64_to_bytes(deflated)
  bytes_inflated = zlib_inflate(bytes_deflated,len(inflated,kind=pI64))
  if (inflated /= transfer(bytes_inflated,inflated)) error stop 'zlib_inflate'

end subroutine selfTest

end module zlib
