program DAMASK_test
  use parallelization
  use HDF5_utilities

  use test_prec
  use test_HDF5_utilities

  call prec_test()

  call parallelization_init()
  call HDF5_utilities_init()

  call HDF5_utilities_test()

end program DAMASK_test
