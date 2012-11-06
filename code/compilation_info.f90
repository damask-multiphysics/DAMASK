!$Id$
#ifdef __GFORTRAN__
#if __GNUC__<=4 &&  __GNUC_MINOR__<=5
    write(6,'(a,i1.1,a,i1.1,a,i2.2))') ' Compiled with GCC version ', __GNUC__,'.',__GNUC_MINOR__,&
                                                                           '.',__GNUC_PATCHLEVEL__
#else
    write(6,*) 'Compiled with ', compiler_version() !not supported by GFORTRAN 4.5 and ifort 12
    write(6,*) 'With options  ', compiler_options()
#endif
#endif 
#ifdef __INTEL_COMPILER
  write(6,'(a,i4.4,a,i8.8)') ' Compiled with Intel fortran version ', __INTEL_COMPILER,&
                                                        ', build date ', __INTEL_COMPILER_BUILD_DATE
#endif
write(6,*) 'Compiled on ', __DATE__,' at ',__TIME__
write(6,*)
flush(6)