!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Wrap FFTW3 into a module.
!--------------------------------------------------------------------------------------------------
module FFTW3
  use, intrinsic :: ISO_C_binding

  implicit none(type,external)
  public

  include 'fftw3-mpi.f03'

end module FFTW3
