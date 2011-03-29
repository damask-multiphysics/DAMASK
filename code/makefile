# Makefile to compile the Material subroutine for BVP solution using spectral method
#
# use switch on make to determine precision, e.g make precision=single
# default is precision=double
# be sure to remove all librarys with different precision (make clean)
# 
# Uses openmp to parallelise the material subroutines (set number of cores with "export MPIE_NUM_THREADS=n" to n)
# Uses linux threads to parallelise fftw3 (should also be possible with openmp)
# Besides of the f90 files written at MPIE, the two library files of fftw3 "libfftw3_threads.a" "libfftw3.a" are also needed
# Install fftw3 (v3.2.2 is tested) with "./configure --enable-threads --enable-float"  and "make", "make install" is not needed
# as long as the two library files are copied to the source code directory.

COMPILE_DOUBLE = -openmp -c -O3 -fpp -heap-arrays 500000000
COMPILE_SINGLE = -openmp -c -O3 -fpp


precision :=double
ifeq ($(precision),single)
cpspectral_single.exe: mpie_spectral_single.o CPFEM.a
	ifort -openmp -o cpspectral_single.exe mpie_spectral_single.o CPFEM.a libfftw3f_threads.a libfftw3f.a constitutive.a advanced.a basics.a -lpthread
mpie_spectral_single.o: mpie_spectral_single.f90 CPFEM.o 
	ifort $(COMPILE_DOUBLE) mpie_spectral_single.f90 

else
cpspectral_double.exe: mpie_spectral.o CPFEM.a
	ifort -openmp -o cpspectral_double.exe mpie_spectral.o CPFEM.a libfftw3_threads.a libfftw3.a constitutive.a advanced.a basics.a -lpthread
mpie_spectral.o: mpie_spectral.f90 CPFEM.o 
	ifort $(COMPILE_DOUBLE) mpie_spectral.f90 
endif

CPFEM.a: CPFEM.o
	ar rc CPFEM.a homogenization.o homogenization_RGC.o homogenization_isostrain.o crystallite.o CPFEM.o constitutive.o

CPFEM.o: CPFEM.f90 homogenization.o
	ifort $(COMPILE_DOUBLE) CPFEM.f90
homogenization.o: homogenization.f90 homogenization_isostrain.o homogenization_RGC.o crystallite.o
	ifort $(COMPILE_DOUBLE) homogenization.f90
homogenization_RGC.o: homogenization_RGC.f90 constitutive.a
	ifort $(COMPILE_DOUBLE) homogenization_RGC.f90
homogenization_isostrain.o: homogenization_isostrain.f90 basics.a advanced.a
	ifort $(COMPILE_DOUBLE) homogenization_isostrain.f90
crystallite.o: crystallite.f90 constitutive.a
	ifort $(COMPILE_DOUBLE) crystallite.f90



constitutive.a: constitutive.o
	ar rc constitutive.a constitutive.o constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o basics.a advanced.a

constitutive.o: constitutive.f90 constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o
	ifort $(COMPILE_DOUBLE) constitutive.f90

constitutive_titanmod.o: constitutive_titanmod.f90 basics.a advanced.a
	ifort $(COMPILE_DOUBLE) constitutive_titanmod.f90

constitutive_nonlocal.o: constitutive_nonlocal.f90 basics.a advanced.a
	ifort $(COMPILE_DOUBLE) constitutive_nonlocal.f90

constitutive_dislotwin.o: constitutive_dislotwin.f90 basics.a advanced.a
	ifort $(COMPILE_DOUBLE) constitutive_dislotwin.f90

constitutive_j2.o: constitutive_j2.f90 basics.a advanced.a
	ifort $(COMPILE_DOUBLE) constitutive_j2.f90

constitutive_phenopowerlaw.o: constitutive_phenopowerlaw.f90 basics.a advanced.a
	ifort $(COMPILE_DOUBLE) constitutive_phenopowerlaw.f90



advanced.a: lattice.o
	ar rc advanced.a FEsolving.o mesh.o material.o lattice.o


lattice.o: lattice.f90 material.o
	ifort $(COMPILE_DOUBLE) lattice.f90
material.o: material.f90 mesh.o
	ifort $(COMPILE_DOUBLE) material.f90
mesh.o: mesh.f90 FEsolving.o
	ifort $(COMPILE_DOUBLE) mesh.f90
FEsolving.o: FEsolving.f90 basics.a
	ifort $(COMPILE_DOUBLE) FEsolving.f90

ifeq ($(precision),single)
basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o mpie_spectral_interface.o prec_single.o
else
basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o mpie_spectral_interface.o prec.o
endif

debug.o: debug.f90 numerics.o
	ifort $(COMPILE_SINGLE) debug.f90
math.o: math.f90 numerics.o
	ifort $(COMPILE_SINGLE) math.f90

numerics.o: numerics.f90 IO.o
	ifort $(COMPILE_SINGLE) numerics.f90
IO.o: IO.f90 mpie_spectral_interface.o
	ifort $(COMPILE_SINGLE) IO.f90

ifeq ($(precision),single)
mpie_spectral_interface.o: mpie_spectral_interface.f90 prec_single.o
	ifort $(COMPILE_SINGLE) mpie_spectral_interface.f90
prec_single.o: prec_single.f90
	ifort $(COMPILE_SINGLE) prec_single.f90
else
mpie_spectral_interface.o: mpie_spectral_interface.f90 prec.o
	ifort $(COMPILE_SINGLE) mpie_spectral_interface.f90
prec.o: prec.f90
	ifort $(COMPILE_SINGLE) prec.f90
endif


clean: 
	rm -rf *.o 
	rm -rf *.mod
	rm -rf *.exe 
	rm -rf basics.a
	rm -rf advanced.a
	rm -rf constitutive.a
	rm -rf CPFEM.a
