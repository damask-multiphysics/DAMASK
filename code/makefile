# Makefile to compile the Material subroutine for BVP solution using spectral method
#
# use switch on make to determine PRECISION, e.g make PRECISION=single
# default is PRECISION=double
# be sure to remove all librarys with different PRECISION (make clean)
# 
# Uses openmp to parallelise the material subroutines (set number of cores with "export MPIE_NUM_THREADS=n" to n)
# Uses linux threads to parallelise fftw3 (should also be possible with openmp)
# Besides of the f90 files written at MPIE, the two library files of fftw3 "libfftw3_threads.a" "libfftw3.a" are also needed
# Install fftw3 (v3.2.2 is tested) with "./configure --enable-threads --enable-float"  and "make", "make install" is not needed
# as long as the two library files are copied to the source code directory.

#standart values
PRECISION :=double  # precision
F90 :=ifort         # compiler (Intel or gfortran, default intel)
VERSION :=10        # version of intel compiler. More aggressive optimization if VERSION =12
PORTABLE :=TRUE     # decision, if executable is optimized for the machine on which it was built

ifeq ($(PORTABLE), FALSE)
PORTABLE_SWITCH = -xHost
endif

ifeq ($(F90), ifort)
OPENMP_FLAG = -openmp
COMPILE_OPTIONS = -fpp -diag-disable 8291   #for preprocessor, switch of messages on formatting of output
OPTIMIZATION_AGGRESSIVE = -O3 -static $(PORTABLE_SWITCH)
endif

ifeq ($(F90), gfortran)
OPENMP_FLAG := -fopenmp
OPTIMIZATION_AGGRESSIVE = -O3
endif

ifeq ($(VERSION), 12)
OPTIMIZATION_DEFENSIVE = -O3 -static $(PORTABLE_SWITCH)
else
OPTIMIZATION_DEFENSIVE = -O2
endif

COMPILE = $(OPENMP_FLAG) $(COMPILE_OPTIONS) $(OPTIMIZATION_AGGRESSIVE) -c
COMPILE_HEAP = $(COMPILE) -heap-arrays 500000000
COMPILE_HEAP_DEFENSIVE = $(OPENMP_FLAG) $(COMPILE_OPTIONS) $(OPTIMIZATION_DEFENSIVE) -c -heap-arrays 500000000

   
ifeq ($(PRECISION),single)
DAMASK_spectral_single.exe: DAMASK_spectral_single.o CPFEM.a
	$(F90) $(OPENMP_FLAG) -o DAMASK_spectral_single.exe DAMASK_spectral_single.o CPFEM.a include/libfftw3f_threads.a include/libfftw3f.a constitutive.a advanced.a basics.a -lpthread
DAMASK_spectral_single.o: DAMASK_spectral_single.f90 CPFEM.o 
	$(F90) $(COMPILE_HEAP_DEFENSIVE) DAMASK_spectral_single.f90 

else
DAMASK_spectral.exe: DAMASK_spectral.o CPFEM.a
	$(F90) $(OPENMP_FLAG) -o DAMASK_spectral.exe DAMASK_spectral.o CPFEM.a include/libfftw3_threads.a include/libfftw3.a constitutive.a advanced.a basics.a -lpthread
DAMASK_spectral.o: DAMASK_spectral.f90 CPFEM.o 
	$(F90) $(COMPILE_HEAP_DEFENSIVE) DAMASK_spectral.f90 
endif

CPFEM.a: CPFEM.o
	ar rc CPFEM.a homogenization.o homogenization_RGC.o homogenization_isostrain.o crystallite.o CPFEM.o constitutive.o

CPFEM.o: CPFEM.f90 homogenization.o
	$(F90) $(COMPILE_HEAP) CPFEM.f90
homogenization.o: homogenization.f90 homogenization_isostrain.o homogenization_RGC.o crystallite.o
	$(F90) $(COMPILE_HEAP) homogenization.f90
homogenization_RGC.o: homogenization_RGC.f90 constitutive.a
	$(F90) $(COMPILE_HEAP) homogenization_RGC.f90
homogenization_isostrain.o: homogenization_isostrain.f90 basics.a advanced.a
	$(F90) $(COMPILE_HEAP) homogenization_isostrain.f90
crystallite.o: crystallite.f90 constitutive.a
	$(F90) $(COMPILE_HEAP) crystallite.f90



constitutive.a: constitutive.o
	ar rc constitutive.a constitutive.o constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o basics.a advanced.a

constitutive.o: constitutive.f90 constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o
	$(F90) $(COMPILE_HEAP) constitutive.f90

constitutive_titanmod.o: constitutive_titanmod.f90 basics.a advanced.a
	$(F90) $(COMPILE_HEAP) constitutive_titanmod.f90

constitutive_nonlocal.o: constitutive_nonlocal.f90 basics.a advanced.a
	$(F90) $(COMPILE_HEAP) constitutive_nonlocal.f90

constitutive_dislotwin.o: constitutive_dislotwin.f90 basics.a advanced.a
	$(F90) $(COMPILE_HEAP) constitutive_dislotwin.f90

constitutive_j2.o: constitutive_j2.f90 basics.a advanced.a
	$(F90) $(COMPILE_HEAP) constitutive_j2.f90

constitutive_phenopowerlaw.o: constitutive_phenopowerlaw.f90 basics.a advanced.a
	$(F90) $(COMPILE_HEAP) constitutive_phenopowerlaw.f90



advanced.a: lattice.o
	ar rc advanced.a FEsolving.o mesh.o material.o lattice.o


lattice.o: lattice.f90 material.o
	$(F90) $(COMPILE_HEAP) lattice.f90
material.o: material.f90 mesh.o
	$(F90) $(COMPILE_HEAP) material.f90
mesh.o: mesh.f90 FEsolving.o
	$(F90) $(COMPILE_HEAP) mesh.f90
FEsolving.o: FEsolving.f90 basics.a
	$(F90) $(COMPILE_HEAP) FEsolving.f90

ifeq ($(PRECISION),single)
basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o DAMASK_spectral_interface.o prec_single.o
else
basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o DAMASK_spectral_interface.o prec.o
endif

debug.o: debug.f90 numerics.o
	$(F90) $(COMPILE) debug.f90
math.o: math.f90 numerics.o
	$(F90) $(COMPILE) math.f90

numerics.o: numerics.f90 IO.o
	$(F90) $(COMPILE) numerics.f90
IO.o: IO.f90 DAMASK_spectral_interface.o
	$(F90) $(COMPILE) IO.f90

ifeq ($(PRECISION),single)
DAMASK_spectral_interface.o: DAMASK_spectral_interface.f90 prec_single.o
	$(F90) $(COMPILE) DAMASK_spectral_interface.f90
prec_single.o: prec_single.f90
	$(F90) $(COMPILE) prec_single.f90
else
DAMASK_spectral_interface.o: DAMASK_spectral_interface.f90 prec.o
	$(F90) $(COMPILE) DAMASK_spectral_interface.f90
prec.o: prec.f90
	$(F90) $(COMPILE) prec.f90
endif


clean: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
