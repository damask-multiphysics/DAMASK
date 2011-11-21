########################################################################################
# Makefile to compile the Material subroutine for BVP solution using spectral method
########################################################################################
# Be sure to remove all files compiled with different options by using "make clean"
# 
# Uses OpenMP to parallelise the material subroutines (set number of cores with "export DAMASK_NUM_THREADS=n" to n)
# Uses linux threads to parallelise fftw3
# Besides of the f90 files written at MPIE, the two library files of fftw3 "libfftw3_threads.a" "libfftw3.a" are also needed
# Install fftw3 (v3.3 is tested) with "./configure --enable-threads --enable-sse2 --enable-shared" and "make"; "make install" is not needed
# as long as the two library files "libfftw3_threads.a" "libfftw3.a" are copied to the /../lib directory.
# Need the AMD Core Math Library to be installed (v 4.4 is tested)
########################################################################################
# OPTIONS = standard (alternative): meaning
#-------------------------------------------------------------
# PRECISION = double (single): floating point precision
# F90 = ifort (gfortran): compiler, choose Intel or GNU
# PORTABLE = TRUE (FALSE): decision, if executable is optimized for the machine on which it was built.
# OPTIMIZATION = DEFENSIVE (OFF,AGGRESSIVE,ULTRA): Optimization mode: O2, O0, O3 + further options for most files, 03 + further options for all files
# OPENMP = TRUE (FALSE): OpenMP multiprocessor support
# ACMLPATH =/opt/acml4.4.0/ifort64/lib (...) Path to ACML Library, choose according to your system
# ACMLPATH =/opt/acml4.4.0/ifort64_mp/lib (...) Path to ACML Library with multicore support, choose according to your system
# FFTWOPTIONS =./../lib/libfftw3.a include/libfftw3_threads.a -lpthread (...) Path to FFTW library files with Linux threads (multicore) support
# FFTWOPTIONS =./../lib/libfftw3.a (...) Path to FFTW library files without Linux threads (multicore) support
# FFTWOPTIONS is different for single and double precision. Choose the options to use OpenMP instead of pthreads support or change the directory
# PREFIX: specify an arbitrary prefix
# SUFFIX: specify an arbitrary suffix
# COMPILERNAME = overwrite name of Compiler, e.g. using mpich-g90 instead of ifort
########################################################################################
# Here are some usefull debugging switches. Switch on by uncommenting the #SUFFIX line at the end of this section:
# information on http://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors/
# check if an array index is too small (<1) or too large!
DEBUG1 =-check bounds -g
#will cause a lot of warnings because we create a bunch of temporary arrays
DEBUG2 =-check arg_temp_created
#check from time to time
DEBUG3 =-fp-stack-check -g -traceback -gen-interfaces -warn interfaces
#should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Problably it helps also to unlimit other limits
DEBUG4 =-heap-arrays
#checks for standard
DEBUG5 =stand std03/std95
#SUFFIX =$(DEBUG1) $(DEBUG2) $(DEBUG3) 
########################################################################################

ifndef F90
F90 =ifort
endif

ifndef OPTIMIZATION
OPTIMIZATION =DEFENSIVE
endif
MAXOPTI =$(OPTIMIZATION)

ifeq ($(OPTIMIZATION),AGGRESSIVE)
MAXOPTI=DEFENSIVE
endif

ifeq ($(OPTIMIZATION),ULTRA)
MAXOPTI=AGGRESSIVE
OPTIMIZATION=AGGRESSIVE
endif

ifeq ($(PORTABLE),FALSE)
PORTABLE_SWITCH =-msse3
endif

ifndef OPENMP
OPENMP=ON
endif

ifeq ($(OPENMP),ON)
OPENMP_FLAG_ifort =-openmp -openmp-report0 -parallel
OPENMP_FLAG_gfortran =-fopenmp
ifndef ACMLPATH
ACMLPATH =/opt/acml4.4.0/ifort64_mp/lib
endif
ifndef FFTWOPTIONS
ifeq ($(PRECISION),single)
FFTWOPTIONS =./../lib/libfftw3f_threads.a ./../lib/libfftw3f.a -lpthread
else
FFTWOPTIONS =./../lib/libfftw3_threads.a ./../lib/libfftw3.a -lpthread
endif
endif
BLAS_ifort =-L $(ACMLPATH) -lacml_mp
BLAS_gfortran =
else
ifndef ACMLPATH
ACMLPATH =/opt/acml4.4.0/ifort64/lib
endif
ifndef FFTWOPTIONS
ifeq ($(PRECISION),single)
FFTWOPTIONS =./../lib/libfftw3f.a
else
FFTWOPTIONS =./../lib/libfftw3.a
endif
endif
BLAS_ifort =-L $(ACMLPATH) -lacml
BLAS_gfortran =
endif


OPTIMIZATION_OFF_ifort =-O0
OPTIMIZATION_OFF_gfortran =-O0
OPTIMIZATION_DEFENSIVE_ifort =-O2
OPTIMIZATION_DEFENSIVE_gfortran =-O2
OPTIMIZATION_AGGRESSIVE_ifort =-O3 $(PORTABLE_SWITCH) -ip -static -fp-model fast=2 -no-prec-div
OPTIMIZATION_AGGRESSIVE_gfortran =-O3 $(PORTABLE_SWITCH) -march=opteron -ffast-math -funroll-loops -ftree-vectorize -ftree-loop-linear 

COMPILE_OPTIONS_ifort =-fpp -diag-disable 8291,8290
COMPILE_OPTIONS_gfortran =-xf95-cpp-input -ffree-line-length-none -fno-range-check

COMPILE =${OPENMP_FLAG_${F90}} ${COMPILE_OPTIONS_${F90}} ${OPTIMIZATION_${OPTIMIZATION}_${F90}} -c
COMPILE_MAXOPTI =${OPENMP_FLAG_${F90}} ${COMPILE_OPTIONS_${F90}} ${OPTIMIZATION_${MAXOPTI}_${F90}} -c

ifndef COMPILERNAME
COMPILERNAME=$(F90)
endif

ifeq ($(PRECISION),single)
DAMASK_spectral_single.exe: DAMASK_spectral_single.o CPFEM.a
	$(PREFIX)	$(COMPILERNAME) ${OPENMP_FLAG_${F90}} -o DAMASK_spectral_single.exe DAMASK_spectral_single.o CPFEM.a $(FFTWOPTIONS)\
  constitutive.a advanced.a basics.a ${BLAS_${F90}}
DAMASK_spectral_single.o: DAMASK_spectral_single.f90 CPFEM.o 
	$(PREFIX)	$(COMPILERNAME) $(COMPILE_MAXOPTI) DAMASK_spectral_single.f90 $(SUFFIX)
else
DAMASK_spectral.exe: DAMASK_spectral.o CPFEM.a
	$(PREFIX)	$(COMPILERNAME) ${OPENMP_FLAG_${F90}} -o DAMASK_spectral.exe DAMASK_spectral.o CPFEM.a $(FFTWOPTIONS)\
  constitutive.a advanced.a basics.a ${BLAS_${F90}}
DAMASK_spectral.o: DAMASK_spectral.f90 CPFEM.o 
	$(PREFIX)	$(COMPILERNAME) $(COMPILE_MAXOPTI) DAMASK_spectral.f90 $(SUFFIX)
endif

CPFEM.a: CPFEM.o
	ar rc CPFEM.a homogenization.o homogenization_RGC.o homogenization_isostrain.o crystallite.o CPFEM.o constitutive.o

CPFEM.o: CPFEM.f90 homogenization.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) CPFEM.f90 $(SUFFIX)
homogenization.o: homogenization.f90 homogenization_isostrain.o homogenization_RGC.o crystallite.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) homogenization.f90 $(SUFFIX)
homogenization_RGC.o: homogenization_RGC.f90 constitutive.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) homogenization_RGC.f90 $(SUFFIX)
homogenization_isostrain.o: homogenization_isostrain.f90 basics.a advanced.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) homogenization_isostrain.f90 $(SUFFIX)
crystallite.o: crystallite.f90 constitutive.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) crystallite.f90 $(SUFFIX)



constitutive.a: constitutive.o
	ar rc constitutive.a constitutive.o constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o basics.a advanced.a

constitutive.o: constitutive.f90 constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) constitutive.f90 $(SUFFIX)

constitutive_titanmod.o: constitutive_titanmod.f90 basics.a advanced.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) constitutive_titanmod.f90 $(SUFFIX)

constitutive_nonlocal.o: constitutive_nonlocal.f90 basics.a advanced.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) constitutive_nonlocal.f90 $(SUFFIX)

constitutive_dislotwin.o: constitutive_dislotwin.f90 basics.a advanced.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) constitutive_dislotwin.f90 $(SUFFIX)

constitutive_j2.o: constitutive_j2.f90 basics.a advanced.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) constitutive_j2.f90 $(SUFFIX)

constitutive_phenopowerlaw.o: constitutive_phenopowerlaw.f90 basics.a advanced.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) constitutive_phenopowerlaw.f90 $(SUFFIX)



advanced.a: lattice.o
	ar rc advanced.a FEsolving.o mesh.o material.o lattice.o


lattice.o: lattice.f90 material.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) lattice.f90 $(SUFFIX)
material.o: material.f90 mesh.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) material.f90 $(SUFFIX)
mesh.o: mesh.f90 FEsolving.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) mesh.f90 $(SUFFIX)
FEsolving.o: FEsolving.f90 basics.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) FEsolving.f90 $(SUFFIX)

ifeq ($(PRECISION),single)
basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o DAMASK_spectral_interface.o prec_single.o
else
basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o DAMASK_spectral_interface.o prec.o
endif

debug.o: debug.f90 numerics.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) debug.f90 $(SUFFIX)
math.o: math.f90 numerics.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) math.f90 $(SUFFIX)
numerics.o: numerics.f90 IO.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) numerics.f90 $(SUFFIX)
IO.o: IO.f90 DAMASK_spectral_interface.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) IO.f90 $(SUFFIX)
  
ifeq ($(PRECISION),single)
DAMASK_spectral_interface.o: DAMASK_spectral_interface.f90 prec_single.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) DAMASK_spectral_interface.f90 $(SUFFIX)
prec_single.o: prec_single.f90
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) prec_single.f90 $(SUFFIX)
else
DAMASK_spectral_interface.o: DAMASK_spectral_interface.f90 prec.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) DAMASK_spectral_interface.f90 $(SUFFIX)
prec.o: prec.f90
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) prec.f90 $(SUFFIX)
endif


clean: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	rm -rf *.exe