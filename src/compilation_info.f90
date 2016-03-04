!##############################################################
!$Id$
#ifdef __GFORTRAN__
  write(6,*) 'Compiled with ', compiler_version() !not supported by and ifort <= 15 (and old gfortran)
  write(6,*) 'With options  ', compiler_options()
#endif 
#ifdef __INTEL_COMPILER
  write(6,'(a,i4.4,a,i8.8)') ' Compiled with Intel fortran version ', __INTEL_COMPILER,&
                                                     ', build date ', __INTEL_COMPILER_BUILD_DATE
#endif
write(6,*) 'Compiled on ', __DATE__,' at ',__TIME__
write(6,*)
flush(6)
