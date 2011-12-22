########################################################################################
# Makefile to compile the Material subroutine for BVP solution using spectral method
########################################################################################
# Be sure to remove all files compiled with different options by using "make clean"
# 
# Uses OpenMP to parallelize the material subroutines (set number of threads with "export DAMASK_NUM_THREADS=n" to n)
#
# Install fftw3 (v3.3 is tested):
# + patch api/f77funcs.h, line 92 in the FFTW source with:
#  FFTW_VOIDFUNC F77(set_timelimit,SET_TIMELIMIT)(double *t)
#  {
#      X(set_timelimit)(*t);
#  }
# + run
#   ./configure --enable-threads --enable-sse2 --enable-shared [-enable-float]
#   make
#   make install
# + specify in the "pathinfo:FFTW" where FFTW was installed.
#   We essentially look for two library files "lib/libfftw3_threads.a" and "lib/libfftw3.a", so you can copy those, for instance,
#   into DAMASK_ROOT/lib/fftw/lib/ and specify "./fftw/" as pathinfo:FFTW
# Use --enable-float in above configure for single precision...
# Uses linux threads to parallelize fftw3
#
# Instead of the AMD Core Math Library a standard "liblapack.a/dylib/etc." can be used by leaving pathinfo:ACML blank
########################################################################################
# OPTIONS = standard (alternative): meaning
#-------------------------------------------------------------
# PRECISION = double (single): floating point precision
# F90 = ifort (gfortran): compiler, choose Intel or GNU
# PORTABLE = TRUE (FALSE): decision, if executable is optimized for the machine on which it was built.
# OPTIMIZATION = DEFENSIVE (OFF,AGGRESSIVE,ULTRA): Optimization mode: O2, O0, O3 + further options for most files, 03 + further options for all files
# OPENMP = TRUE (FALSE): OpenMP multiprocessor support
# FFTWROOT = pathinfo:FFTW (will be adjusted by setup_code.py)
# ACMLROOT = pathinfo:ACML (will be adjusted by setup_code.py)
# ACMLPATH =$(ACMLROOT)/"compilerdir"/lib (...) Path to ACML, choose according to your system
# ACMLPATH =$(ACMLROOT/"compilerdir"_mp/lib (...) Path to ACML with multicore support, choose according to your system
# "compilerdir" is "intel64" for ifort
# PREFIX: specify an arbitrary prefix
# SUFFIX: specify an arbitrary suffix
# COMPILERNAME = overwrite name of Compiler, e.g. using mpich-g90 instead of ifort
########################################################################################
# Here are some useful debugging switches. Switch on by uncommenting the #SUFFIX line at the end of this section:
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
DEBUG5 =-stand std03/std95

#SUFFIX =$(DEBUG1) $(DEBUG2) $(DEBUG3)
########################################################################################

#auto values will be set by setup_code.py
FFTWROOT := 
IKMLROOT := 
ACMLROOT := 
LAPACKROOT := 

F90 ?= ifort
COMPILERNAME ?= $(F90)
OPENMP ?= ON
OPTIMIZATION ?= DEFENSIVE

ifeq "$(OPTIMIZATION)" "OFF"
OPTI    := OFF
MAXOPTI := OFF
endif
ifeq "$(OPTIMIZATION)" "DEFENSIVE"
OPTI    := DEFENSIVE
MAXOPTI := DEFENSIVE
endif
ifeq "$(OPTIMIZATION)" "AGGRESSIVE"
OPTI    := AGGRESSIVE
MAXOPTI := DEFENSIVE
endif
ifeq "$(OPTIMIZATION)" "ULTRA"
OPTI    := AGGRESSIVE
MAXOPTI := AGGRESSIVE
endif

ifndef OPTI
OPTI    := DEFENSIVE
MAXOPTI := DEFENSIVE
endif

ifeq  "$(PORTABLE)" "FALSE"
PORTABLE_SWITCH = -msse3
endif

ifeq "$(PRECISION)" "single"
FFTWPREC = f
else
FFTWPREC =
endif


# settings for multicore support
ifeq "$(OPENMP)" "ON"
OPENMP_FLAG_ifort = -openmp -openmp-report0 -parallel
OPENMP_FLAG_gfortran = -fopenmp
ACML_ARCH =_mp
LIBRARIES += -lfftw3$(FFTWPREC)_threads -lpthread
endif

LIBRARIES += -lfftw3$(FFTWPREC)
LIB_DIRS += -L$(FFTWROOT)/lib

ifdef IKMLROOT
LIBRARIES += -mkl
else
ifdef ACMLROOT
LIB_DIRS += -L$(ACMLROOT)/$(F90)64$(ACML_ARCH)/lib
LIBRARIES += -lacml$(ACML_ARCH)
else
ifdef LAPACKROOT
LIB_DIRS += -L$(LAPACKROOT)/lib -L$(LAPACKROOT)/lib64
LIBRARIES += -llapack
endif
endif
endif


OPTIMIZATION_OFF_ifort           :=-O0
OPTIMIZATION_OFF_gfortran        :=-O0
OPTIMIZATION_DEFENSIVE_ifort     :=-O2
OPTIMIZATION_DEFENSIVE_gfortran  :=-O2
OPTIMIZATION_AGGRESSIVE_ifort    :=-O3 $(PORTABLE_SWITCH) -ip -static -fp-model fast=2 -no-prec-div
OPTIMIZATION_AGGRESSIVE_gfortran :=-O3 $(PORTABLE_SWITCH) -ffast-math -funroll-loops -ftree-vectorize 


COMPILE_OPTIONS_ifort    := -fpp -diag-disable 8291,8290
COMPILE_OPTIONS_gfortran := -xf95-cpp-input -fno-range-check

COMPILE         = $(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(OPTIMIZATION_$(OPTI)_$(F90)) -c
COMPILE_MAXOPTI = $(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90)) -c


ifeq "$(PRECISION)" "single"
DAMASK_spectral_single.exe: DAMASK_spectral_single.o CPFEM.a
	$(PREFIX)	$(COMPILERNAME) $(OPENMP_FLAG_$(F90)) -o DAMASK_spectral_single.exe DAMASK_spectral_single.o CPFEM.a \
  constitutive.a advanced.a basics.a $(LIB_DIRS) $(LIBRARIES)
DAMASK_spectral_single.o: DAMASK_spectral_single.f90 CPFEM.o 
	$(PREFIX)	$(COMPILERNAME) $(COMPILE_MAXOPTI) DAMASK_spectral_single.f90 $(SUFFIX)
else
DAMASK_spectral.exe: DAMASK_spectral.o CPFEM.a
	$(PREFIX)	$(COMPILERNAME) ${OPENMP_FLAG_${F90}} -o DAMASK_spectral.exe DAMASK_spectral.o CPFEM.a \
  constitutive.a advanced.a basics.a $(LIB_DIRS) $(LIBRARIES)
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

ifeq "$(PRECISION)" "single"
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
  
ifeq "$(PRECISION)" "single"
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


tidy: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a

clean: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	rm -rf *.exe

