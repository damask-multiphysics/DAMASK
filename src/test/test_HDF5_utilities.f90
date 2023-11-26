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
  real(pREAL), dimension(3) :: d1_in,d1_out
  real(pREAL), dimension(3,3) :: d2_in,d2_out
  real(pREAL), dimension(3,3,3) :: d3_in,d3_out
  real(pREAL), dimension(3,3,3,3) :: d4_in,d4_out
  real(pREAL), dimension(3,3,3,3,3) :: d5_in,d5_out


  call random_number(d1_in)
  call random_number(d2_in)
  call random_number(d3_in)
  call random_number(d4_in)
  call random_number(d5_in)

  f = HDF5_openFile('test.hdf5','w')

  call HDF5_write(d1_in,f,'d1')
  call HDF5_write(d2_in,f,'d2')
  call HDF5_write(d3_in,f,'d3')
  call HDF5_write(d4_in,f,'d4')
  call HDF5_write(d5_in,f,'d5')

  call HDF5_read(d1_out,f,'d1')
  call HDF5_read(d2_out,f,'d2')
  call HDF5_read(d3_out,f,'d3')
  call HDF5_read(d4_out,f,'d4')
  call HDF5_read(d5_out,f,'d5')

  if (any(d1_in /= d1_out)) error stop 'test_read_write(w)/d1'
  if (any(d2_in /= d2_out)) error stop 'test_read_write(w)/d2'
  if (any(d3_in /= d3_out)) error stop 'test_read_write(w)/d3'
  if (any(d4_in /= d4_out)) error stop 'test_read_write(w)/d4'
  if (any(d5_in /= d5_out)) error stop 'test_read_write(w)/d5'

  call HDF5_closeFile(f)


  f = HDF5_openFile('test.hdf5','r')

  call HDF5_read(d1_out,f,'d1')
  call HDF5_read(d2_out,f,'d2')
  call HDF5_read(d3_out,f,'d3')
  call HDF5_read(d4_out,f,'d4')
  call HDF5_read(d5_out,f,'d5')

  if (any(d1_in /= d1_out)) error stop 'test_read_write(r)/d1'
  if (any(d2_in /= d2_out)) error stop 'test_read_write(r)/d2'
  if (any(d3_in /= d3_out)) error stop 'test_read_write(r)/d3'
  if (any(d4_in /= d4_out)) error stop 'test_read_write(r)/d4'
  if (any(d5_in /= d5_out)) error stop 'test_read_write(r)/d5'

  call HDF5_closeFile(f)


end subroutine read_write

end module test_HDF5_utilities
