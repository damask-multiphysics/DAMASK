module test_HDF5_utilities
  use prec
  use HDF5
  use HDF5_utilities

  implicit none(type,external)

  private
  public :: HDF5_utilities_test

  contains

subroutine HDF5_utilities_test()

  print*, 'begin test HDF5_utilities'
  call test_read_write()
  print*, 'end test HDF5_utilities'

end subroutine HDF5_utilities_test


subroutine test_read_write()

  integer(HID_T) :: f
  real(pREAL), dimension(3) :: d_in,d_out


  call random_number(d_in)

  f = HDF5_openFile('test.hdf5','w')

  call HDF5_write(d_in,f,'test')
  call HDF5_read(d_out,f,'test')

  if (any(d_in /= d_out)) error stop 'test_read_write'

end subroutine test_read_write

end module test_HDF5_utilities
