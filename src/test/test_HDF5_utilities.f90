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

  integer(HID_T) :: f, create_list
  integer :: hdferr, order

  real(pREAL), dimension(3) :: real_d1_in,real_d1_out
  real(pREAL), dimension(3,3) :: real_d2_in,real_d2_out
  real(pREAL), dimension(3,3,3) :: real_d3_in,real_d3_out
  real(pREAL), dimension(3,3,3,3) :: real_d4_in,real_d4_out
  real(pREAL), dimension(3,3,3,3,3) :: real_d5_in,real_d5_out

  integer, dimension(3) :: int_d1_in,int_d1_out
  integer, dimension(3,3) :: int_d2_in,int_d2_out
  integer, dimension(3,3,3) :: int_d3_in,int_d3_out
  integer, dimension(3,3,3,3) :: int_d4_in,int_d4_out
  integer, dimension(3,3,3,3,3) :: int_d5_in,int_d5_out


  call random_number(real_d1_in)
  call random_number(real_d2_in)
  call random_number(real_d3_in)
  call random_number(real_d4_in)
  call random_number(real_d5_in)

  int_d1_in = int(real_d1_in*2048._pREAL)
  int_d2_in = int(real_d2_in*2048._pREAL)
  int_d3_in = int(real_d3_in*2048._pREAL)
  int_d4_in = int(real_d4_in*2048._pREAL)
  int_d5_in = int(real_d5_in*2048._pREAL)


  f = HDF5_openFile('test.hdf5','w')


  call HDF5_write(real_d1_in,f,'real_d1')
  call HDF5_write(real_d2_in,f,'real_d2')
  call HDF5_write(real_d3_in,f,'real_d3')
  call HDF5_write(real_d4_in,f,'real_d4')
  call HDF5_write(real_d5_in,f,'real_d5')

  call HDF5_write(int_d1_in,f,'int_d1')
  call HDF5_write(int_d2_in,f,'int_d2')
  call HDF5_write(int_d3_in,f,'int_d3')
  call HDF5_write(int_d4_in,f,'int_d4')
  call HDF5_write(int_d5_in,f,'int_d5')


  call HDF5_read(real_d1_out,f,'real_d1')
  call HDF5_read(real_d2_out,f,'real_d2')
  call HDF5_read(real_d3_out,f,'real_d3')
  call HDF5_read(real_d4_out,f,'real_d4')
  call HDF5_read(real_d5_out,f,'real_d5')

  call HDF5_read(int_d1_out,f,'int_d1')
  call HDF5_read(int_d2_out,f,'int_d2')
  call HDF5_read(int_d3_out,f,'int_d3')
  call HDF5_read(int_d4_out,f,'int_d4')
  call HDF5_read(int_d5_out,f,'int_d5')


  if (any(real_d1_in /= real_d1_out)) error stop 'test_read_write(w)/real_d1'
  if (any(real_d2_in /= real_d2_out)) error stop 'test_read_write(w)/real_d2'
  if (any(real_d3_in /= real_d3_out)) error stop 'test_read_write(w)/real_d3'
  if (any(real_d4_in /= real_d4_out)) error stop 'test_read_write(w)/real_d4'
  if (any(real_d5_in /= real_d5_out)) error stop 'test_read_write(w)/real_d5'

  if (any(int_d1_in /= int_d1_out)) error stop 'test_read_write(w)/int_d1'
  if (any(int_d2_in /= int_d2_out)) error stop 'test_read_write(w)/int_d2'
  if (any(int_d3_in /= int_d3_out)) error stop 'test_read_write(w)/int_d3'
  if (any(int_d4_in /= int_d4_out)) error stop 'test_read_write(w)/int_d4'
  if (any(int_d5_in /= int_d5_out)) error stop 'test_read_write(w)/int_d5'


  call HDF5_closeFile(f)

  f = HDF5_openFile('test.hdf5','r')
  call H5Fget_create_plist_f(f,create_list,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Pget_link_creation_order_f(create_list,order,hdferr)
  call HDF5_chkerr(hdferr)
  ! https://github.com/HDFGroup/hdf5/issues/5183
  !if (iand(order,H5P_CRT_ORDER_INDEXED_F) /= H5P_CRT_ORDER_INDEXED_F) error stop 'CRT_ORDER_INDEXED'
  !if (iand(order,H5P_CRT_ORDER_TRACKED_F) /= H5P_CRT_ORDER_TRACKED_F) error stop 'CRT_ORDER_TRACKED'
  call H5Pclose_f(create_list,hdferr)
  call HDF5_chkerr(hdferr)

  call HDF5_read(real_d1_out,f,'real_d1')
  call HDF5_read(real_d2_out,f,'real_d2')
  call HDF5_read(real_d3_out,f,'real_d3')
  call HDF5_read(real_d4_out,f,'real_d4')
  call HDF5_read(real_d5_out,f,'real_d5')

  call HDF5_read(int_d1_out,f,'int_d1')
  call HDF5_read(int_d2_out,f,'int_d2')
  call HDF5_read(int_d3_out,f,'int_d3')
  call HDF5_read(int_d4_out,f,'int_d4')
  call HDF5_read(int_d5_out,f,'int_d5')


  if (any(real_d1_in /= real_d1_out)) error stop 'test_read_write(r)/real_d1'
  if (any(real_d2_in /= real_d2_out)) error stop 'test_read_write(r)/real_d2'
  if (any(real_d3_in /= real_d3_out)) error stop 'test_read_write(r)/real_d3'
  if (any(real_d4_in /= real_d4_out)) error stop 'test_read_write(r)/real_d4'
  if (any(real_d5_in /= real_d5_out)) error stop 'test_read_write(r)/real_d5'

  if (any(int_d1_in /= int_d1_out)) error stop 'test_read_write(r)/int_d1'
  if (any(int_d2_in /= int_d2_out)) error stop 'test_read_write(r)/int_d2'
  if (any(int_d3_in /= int_d3_out)) error stop 'test_read_write(r)/int_d3'
  if (any(int_d4_in /= int_d4_out)) error stop 'test_read_write(r)/int_d4'
  if (any(int_d5_in /= int_d5_out)) error stop 'test_read_write(r)/int_d5'


  call HDF5_closeFile(f)


end subroutine read_write

end module test_HDF5_utilities
