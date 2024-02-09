program DAMASK_test

  use parallelization
  use HDF5_utilities
  use IO

  use test_prec
  use test_system_routines
  use test_misc
  use test_math
  use test_polynomials
  use test_tables
  use test_crystal
  use test_rotations
  use test_IO
  use test_types
  use test_HDF5_utilities

  external :: quit

  character(len=*), parameter :: &
    ok  = achar(27)//'[32mok'//achar(27)//'[0m', &
    fmt = '(3x,a,T20,a,1x)'

  call parallelization_init()
  call HDF5_utilities_init()

  write(IO_STDOUT,fmt='(/,1x,a,/)') achar(27)//'[1m'//'testing'//achar(27)//'[0m'

  write(IO_STDOUT,fmt=fmt, advance='no') 'prec','...'
  call test_prec_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'misc','...'
  call test_misc_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'system_routines','...'
  call test_system_routines_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'math','...'
  call test_math_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'polynomials','...'
  call test_polynomials_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'tables','...'
  call test_tables_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'crystal','...'
  call test_crystal_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'rotations','...'
  call test_rotations_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'IO','...'
  call test_IO_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'types','...'
  call test_types_run()
  write(IO_STDOUT,fmt='(a)') ok

  write(IO_STDOUT,fmt=fmt, advance='no') 'HDF5_utilities','...'
  call test_HDF5_utilities_run()
  write(IO_STDOUT,fmt='(a)') ok

  call quit(0)

end program DAMASK_test
