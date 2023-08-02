program DAMASK_test

  use parallelization
  use HDF5_utilities

  use test_prec
  use test_crystal
  use test_IO
  use test_rotations
  use test_misc
  use test_tables
  use test_HDF5_utilities

  external :: quit

  call parallelization_init()

  print('(/,a)'), 'begin test prec'
  call test_prec_run()
  print('(a)'), 'begin test prec'

  print('(/,a)'), 'begin test tables'
  call test_tables_run()
  print('(a)'), 'end test tables'

  print('(/,a)'), 'begin test crystal'
  call test_crystal_run()
  print('(a)'), 'end test crystal'

  print('(/,a)'), 'begin test rotations'
  call test_rotations_run()
  print('(a)'), 'end test rotations'

  print('(/,a)'), 'begin test IO'
  call test_IO_run()
  print('(a)'), 'end test IO'

  print('(/,a)'), 'begin test misc'
  call test_misc_run()
  print('(a)'), 'end test misc'

  call HDF5_utilities_init()

  print('(/,a)'), 'begin test HDF5_utilities'
  call test_HDF5_utilities_run()
  print('(a)'), 'begin test HDF5_utilities'

  call quit(0)

end program DAMASK_test
