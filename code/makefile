cpspectral.out: mpie_spectral.o CPFEM.a
	ifort -o cpspectral.out mpie_spectral.o CPFEM.a fft.o libfftw3.a constitutive.a advanced.a basics.a

mpie_spectral.o: mpie_spectral.f90 CPFEM.o 
	ifort -c -heap-arrays 500000000 mpie_spectral.f90
	
CPFEM.a: CPFEM.o
	ar rc CPFEM.a homogenization.o homogenization_RGC.o homogenization_isostrain.o crystallite.o CPFEM.o constitutive.o

CPFEM.o: CPFEM.f90 homogenization.o
	ifort -c -heap-arrays 500000000 CPFEM.f90
homogenization.o: homogenization.f90 homogenization_isostrain.o homogenization_RGC.o crystallite.o
	ifort -c -heap-arrays 500000000 homogenization.f90
homogenization_RGC.o: homogenization_RGC.f90 constitutive.a
	ifort -c -heap-arrays 500000000 homogenization_RGC.f90
homogenization_isostrain.o: homogenization_isostrain.f90 basics.a advanced.a
	ifort -c -heap-arrays 500000000 homogenization_isostrain.f90
crystallite.o: crystallite.f90 constitutive.a
	ifort -c -heap-arrays 500000000 crystallite.f90



constitutive.a: constitutive.o
	ar rc constitutive.a constitutive.o constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o basics.a advanced.a

constitutive.o: constitutive.f90 constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o
	ifort -c -heap-arrays 500000000 constitutive.f90

constitutive_titanmod.o: constitutive_titanmod.f90 basics.a advanced.a
	ifort -c -heap-arrays 500000000 constitutive_titanmod.f90

constitutive_nonlocal.o: constitutive_nonlocal.f90 basics.a advanced.a
	ifort -c -heap-arrays 500000000 constitutive_nonlocal.f90

constitutive_dislotwin.o: constitutive_dislotwin.f90 basics.a advanced.a
	ifort -c -heap-arrays 500000000 constitutive_dislotwin.f90

constitutive_j2.o: constitutive_j2.f90 basics.a advanced.a
	ifort -c -heap-arrays 500000000 constitutive_j2.f90

constitutive_phenopowerlaw.o: constitutive_phenopowerlaw.f90 basics.a advanced.a
	ifort -c -heap-arrays 500000000 constitutive_phenopowerlaw.f90



advanced.a: lattice.o
	ar rc advanced.a FEsolving.o mesh.o material.o lattice.o

lattice.o: lattice.f90 material.o
	ifort -c -heap-arrays 500000000 lattice.f90
material.o: material.f90 mesh.o
	ifort -c -heap-arrays 500000000 material.f90
mesh.o: mesh.f90 FEsolving.o
	ifort -c -heap-arrays 500000000 mesh.f90
FEsolving.o: FEsolving.f90 basics.a
	ifort -c -heap-arrays 500000000 FEsolving.f90



basics.a: debug.o math.o
	ar rc basics.a debug.o math.o numerics.o IO.o mpie_spectral_interface.o prec.o

debug.o: debug.f90 numerics.o
	ifort -c debug.f90

math.o: math.f90 numerics.o
	ifort -c math.f90

numerics.o: numerics.f90 IO.o
	ifort -c numerics.f90
IO.o: IO.f90 mpie_spectral_interface.o
	ifort -c IO.f90
mpie_spectral_interface.o: mpie_spectral_interface.f90 prec.o
	ifort -c mpie_spectral_interface.f90
prec.o: prec.f90
	ifort -c prec.f90
