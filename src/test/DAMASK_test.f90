program DAMASK_test

  use parallelization
  use HDF5_utilities

  use test_prec
  use test_tables
  use test_crystal
  use test_HDF5_utilities

  print('(/,a)'), 'begin test prec'
  call test_prec_run()
  print('(a)'), 'begin test prec'

  print('(/,a)'), 'begin test tables'
  call test_tables_run()
  print('(a)'), 'end test tables'

  print('(/,a)'), 'begin test crystal'
  call test_crystal_run()
  print('(a)'), 'end test crystal'

  call parallelization_init()
  call HDF5_utilities_init()

  print('(/,a)'), 'begin test HDF5_utilities'
  call test_HDF5_utilities_run()
  print('(a)'), 'begin test HDF5_utilities'

end program DAMASK_test
