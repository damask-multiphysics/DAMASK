module HDF5_io
  use prec
  use IO
  use hdf5

contains

subroutine HDF5_init(filename, total_inc, total_time)
  integer(pInt), intent(in) :: total_inc
  real(pReal),   intent(in) :: total_time

  write(6,*) 'pretend to write something'

end subroutine HDF5_init

end module HDF5_io