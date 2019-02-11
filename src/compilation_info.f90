! https://github.com/jeffhammond/HPCInfo/blob/master/docs/Preprocessor-Macros.md
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
  write(6,*) 'Compiled with ', compiler_version()
  write(6,*) 'With options  ', compiler_options()
#elif defined(__INTEL_COMPILER)
  write(6,'(a,i4.4,a,i8.8)') ' Compiled with Intel fortran version ', __INTEL_COMPILER,&
                                                     ', build date ', __INTEL_COMPILER_BUILD_DATE
#elif defined(__PGI)
  write(6,'(a,i4.4,a,i8.8)') ' Compiled with PGI fortran version ', __PGIC__,&
                                                               '.', __PGIC_MINOR__
#endif
write(6,*) 'Compiled on ', __DATE__,' at ',__TIME__
write(6,*)
flush(6)
