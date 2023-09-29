module test_HDF5_utilities
  use prec
  use HDF5
  use HDF5_utilities

  implicit none(type,external)

  private
  public :: test_HDF5_utilities_run

  contains

subroutine test_HDF5_utilities_run()

  call read_write()

end subroutine test_HDF5_utilities_run


subroutine read_write()

  integer(HID_T) :: f
  real(pREAL), dimension(3) :: d_in,d_out


  call random_number(d_in)

  f = HDF5_openFile('test.hdf5','w')

  call HDF5_write(d_in,f,'test')
  call HDF5_read(d_out,f,'test')

  if (any(d_in /= d_out)) error stop 'test_read_write'

end subroutine read_write

end module test_HDF5_utilities
