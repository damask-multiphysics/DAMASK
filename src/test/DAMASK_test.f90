program DAMASK_test

  use parallelization
  use HDF5_utilities
  use IO

  use test_prec
  use test_crystal
  use test_IO
  use test_rotations
  use test_misc
  use test_tables
  use test_HDF5_utilities

  external :: quit

  character(len=*), parameter :: &
    tab = '19', &
    ok  = achar(27)//'[32mok'//achar(27)//'[0m'

  call parallelization_init()
  call HDF5_utilities_init()

  write(IO_STDOUT,fmt='(/,1x,a,/)') achar(27)//'[1m'//'testing'//achar(27)//'[0m'

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'prec','...'
  call test_prec_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'tables','...'
  call test_tables_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'crystal','...'
  call test_crystal_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'rotations','...'
  call test_rotations_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'IO','...'
  call test_IO_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'misc','...'
  call test_misc_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  write(IO_STDOUT,fmt='(3x,a,T'//tab//',a)', advance='no') 'HDF5_utilities','...'
  call test_HDF5_utilities_run()
  write(IO_STDOUT,fmt='(1x,a)') ok

  call quit(0)

end program DAMASK_test
