########################################################################################
# Makefile to compile the Material subroutine for BVP solution using spectral method
########################################################################################
# Be sure to remove all files compiled with different options by using "make clean"
# 
# Uses OpenMP to parallelize the material subroutines (set number of threads with "export DAMASK_NUM_THREADS=n" to n)
#
# Install fftw3 (v3.3 is tested):
# + run
#   ./configure --enable-threads --enable-sse2 --enable-shared [-enable-float]
#   make
#   make install
# + specify in the "pathinfo:FFTW" where FFTW was installed.
#   We essentially look for two library files "lib/libfftw3_threads" and "lib/libfftw3", so you can copy those, for instance,
#   into DAMASK_ROOT/lib/fftw/lib/ and specify "./fftw/" as pathinfo:FFTW
# Use --enable-float in above configure for single precision...
# Uses linux threads to parallelize fftw3
#
# Instead of the AMD Core Math Library a standard "liblapack.a/dylib/etc." can be used by leaving pathinfo:ACML and pathinfo:IKML blank
########################################################################################
# OPTIONS = standard (alternative): meaning
#-------------------------------------------------------------
# F90            = ifort (gfortran): compiler, choose Intel or GNU
# COMPILERNAME   = overwrite name of Compiler, e.g. using mpich-g90 instead of ifort
# PORTABLE       = TRUE (FALSE): decision, if executable is optimized for the machine on which it was built.
# OPTIMIZATION   = DEFENSIVE (OFF,AGGRESSIVE,ULTRA): Optimization mode: O2, O0, O3 + further options for most files, O3 + further options for all files
# OPENMP         = TRUE (FALSE): OpenMP multiprocessor support
# FFTWROOT       = pathinfo:FFTW   (will be adjusted by setup_code.py - required in pathinfo)
# IKMLROOT       = pathinfo:IKML   (will be adjusted by setup_code.py if present in pathinfo)
# ACMLROOT       = pathinfo:ACML   (will be adjusted by setup_code.py if present in pathinfo)
# LAPACKROOT     = pathinfo:LAPACK (will be adjusted by setup_code.py if present in pathinfo)
# PREFIX         = arbitrary prefix
# SUFFIX         = arbitrary suffix
# STANDARD_CHECK = checking for Fortran 2008, compiler dependend
########################################################################################
# Here are some useful debugging switches for ifort. Switch on by uncommenting the #SUFFIX line at the end of this section:
# information on http://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors/

# check if an array index is too small (<1) or too large!
DEBUG1 =-check bounds -g

#will cause a lot of warnings because we create a bunch of temporary arrays
DEBUG2 =-check arg_temp_created

#check from time to time
DEBUG3 =-fp-stack-check -g -traceback -gen-interfaces -warn interfaces

#should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
DEBUG4 =-heap-arrays

#additional warnings
DEBUG5 =-warn all
# or one or more of those: alignments, declarations,general, ignore_loc, uncalled, unuses, usage

#set precision (check for missing _pInt and _pReal)
DEBUG6= -real-size 32  -integer-size 16
#or one of those 16/32/64/128 (= 2,4,8,16 bytes)
#SUFFIX =$(DEBUG1) $(DEBUG2) $(DEBUG3) $(DEBUG4) $(DEBUG5) $(DEBUG6)

# Here are some useful debugging switches for gfortran
# fcheck-bounds: eqv to DEBUG1 of ifort


########################################################################################

#auto values will be set by setup_code.py
FFTWROOT := $(DAMASK_ROOT)/lib/fftw
IKMLROOT :=
ACMLROOT := /opt/acml4.4.0
LAPACKROOT :=

F90 ?= ifort
COMPILERNAME ?= $(F90)
OPENMP ?= ON
OPTIMIZATION ?= DEFENSIVE

ifeq "$(F90)" "ifort"
ARCHIVE_COMMAND :=xiar
else
ARCHIVE_COMMAND :=ar
endif

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
PORTABLE_SWITCH =-msse3
endif


# settings for multicore support
ifeq "$(OPENMP)" "ON"
OPENMP_FLAG_ifort =-openmp -openmp-report0 -parallel
OPENMP_FLAG_gfortran =-fopenmp
ACML_ARCH =_mp
LIBRARIES +=-lfftw3_threads -lpthread
endif

LIBRARIES +=-lfftw3
LIB_DIRS +=-L$(FFTWROOT)/lib

ifdef IKMLROOT
LIBRARIES +=-mkl
else
ifdef ACMLROOT
LIB_DIRS +=-L$(ACMLROOT)/$(F90)64$(ACML_ARCH)/lib
LIBRARIES +=-lacml$(ACML_ARCH)
else
ifdef LAPACKROOT
LIB_DIRS +=-L$(LAPACKROOT)/lib64 -L$(LAPACKROOT)/lib
LIBRARIES +=-llapack 
endif
endif
endif

ifdef STANDARD_CHECK
STANDARD_CHECK_ifort =$(STANDARD_CHECK)
STANDARD_CHECK_gfortran =$(STANDARD_CHECK)
endif
STANDARD_CHECK_ifort ?=-stand f08 -standard-semantics
STANDARD_CHECK_gfortran ?=-std=f2008


OPTIMIZATION_OFF_ifort           :=-O0
OPTIMIZATION_OFF_gfortran        :=-O0
OPTIMIZATION_DEFENSIVE_ifort     :=-O2
OPTIMIZATION_DEFENSIVE_gfortran  :=-O2
OPTIMIZATION_AGGRESSIVE_ifort    :=-O3 $(PORTABLE_SWITCH) -ip -static -fp-model fast=2 -no-prec-div
OPTIMIZATION_AGGRESSIVE_gfortran :=-O3 $(PORTABLE_SWITCH) -ffast-math -funroll-loops -ftree-vectorize 


COMPILE_OPTIONS_ifort    :=-fpp\
                           -diag-enable sc3\
                           -diag-disable 8291,8290,5268\
                           -warn declarations\
                           -warn general\
                           -warn usage

#alignments:       Determines whether warnings occur for data that is not naturally aligned.
#declarations:     Determines whether warnings occur for any undeclared names.
#errors:           Determines whether warnings are changed to errors.
#general:          Determines whether warning messages and informational messages are issued by the compiler.
#ignore_loc:       Determines whether warnings occur when %LOC is stripped from an actual argument.
#interfaces:       Determines whether the compiler checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks.
#stderrors:        Determines whether warnings about Fortran standard violations are changed to errors.
#truncated_source: Determines whether warnings occur when source exceeds the maximum column width in fixed-format files.
#uncalled:         Determines whether warnings occur when a statement function is never called
#unused:           Determines whether warnings occur for declared variables that are never used.
#usage:            Determines whether warnings occur for questionable programming practices.                         
                           
#-fpp:                 preprocessor
#-diag-disable:        disables warnings, where                           
#   warning ID 9291:
#   warning ID 8290:
#   warning ID 5268: The text exceeds right hand column allowed on the line (we have only comments there)
COMPILE_OPTIONS_gfortran :=-xf95-cpp-input\
                           -ffree-line-length-132\
                           -fno-range-check\
                           -fimplicit-none\
                           -pedantic\
                           -Warray-bounds\
                           -Wunused-parameter\
                           -Wampersand\
                           -Wno-tabs\
                           -Wcharacter-truncation\
                           -Wintrinsic-shadow\
                           -Waliasing\
                           -Wconversion\
                           -Wsurprising\
                           -Wunused-value\
                           -Wunderflow
                            
#-xf95-cpp-input:         preprocessor
#-ffree-line-length-132:  restrict line length to the standard 132 characters
#-fno-range-check:        disables checking if result can be represented by variable. Needs to be set to enable DAMASK_NaN
#-fimplicit-none:         assume "implicit-none" even if not present in source
#-pedantic:               more strict on standard, enables some of the warnings below
#-Warray-bounds:          checks if array reference is out of bounds at compile time. use -fcheck-bounds to also check during runtime
#-Wunused-parameter:      find usused variables with "parameter" attribute
#-Wampersand:             checks if a character expression is continued proberly by an ampersand at the end of the line and at the beginning of the new line
#-Wno-tabs:               do not allow tabs in source
#-Wcharacter-truncation:  warn if character expressions (strings) are truncated
#-Wintrinsic-shadow:      warn if a user-defined procedure or module procedure has the same name as an intrinsic
#-Waliasing:              warn about possible aliasing of dummy arguments. Specifically, it warns if the same actual argument is associated with a dummy argument with "INTENT(IN)" and a dummy argument with "INTENT(OUT)" in a call with an explicit interface. 
#-Wconversion:            warn about implicit conversions between different type
#-Wsurprising:            warn when "suspicious" code constructs are encountered. While technically legal these usually indicate that an error has been made. 
#-Wunused-value:
#-Wunderflow:             produce a warning when numerical constant expressions are encountered, which yield an UNDERFLOW during compilation



#MORE OPTIONS
   
# only for gfortran 4.6:
                            #-Wsuggest-attribute=const
                            #-Wsuggest-attribute=noreturn
                            #-Wsuggest-attribute=pure  
                            
# too many warnings because we have comments beyond character 132:
                            #-Wline-truncation
# warnings because of "flush" is not longer in the standard, but still an intrinsic fuction of the compilers:
                            #-Wintrinsic-std
# warnings because we have many temporary arrays (performance issue?):
                            #-Warray-temporaries
                            
# -Wimplicit-interface
# -pedantic-errors
# -fmodule-private


COMPILE         =$(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(STANDARD_CHECK_$(F90)) $(OPTIMIZATION_$(OPTI)_$(F90)) -c
COMPILE_MAXOPTI =$(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(STANDARD_CHECK_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90)) -c


DAMASK_spectral.exe: DAMASK_spectral.o CPFEM.a
	$(PREFIX)	$(COMPILERNAME) ${OPENMP_FLAG_${F90}} -o DAMASK_spectral.exe DAMASK_spectral.o CPFEM.a \
  constitutive.a advanced.a basics.a $(LIB_DIRS) $(LIBRARIES)

DAMASK_spectral.o: DAMASK_spectral.f90 CPFEM.o 
	$(PREFIX)	$(COMPILERNAME) $(COMPILE_MAXOPTI) DAMASK_spectral.f90 $(SUFFIX)



CPFEM.a: CPFEM.o
	$(ARCHIVE_COMMAND) rc CPFEM.a homogenization.o homogenization_RGC.o homogenization_isostrain.o crystallite.o CPFEM.o constitutive.o

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
	$(ARCHIVE_COMMAND) rc constitutive.a constitutive.o constitutive_titanmod.o constitutive_nonlocal.o constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o basics.a advanced.a

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
	$(ARCHIVE_COMMAND) rc advanced.a FEsolving.o mesh.o material.o lattice.o

lattice.o: lattice.f90 material.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) lattice.f90 $(SUFFIX)

material.o: material.f90 mesh.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) material.f90 $(SUFFIX)

mesh.o: mesh.f90 FEsolving.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) mesh.f90 $(SUFFIX)

FEsolving.o: FEsolving.f90 basics.a
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) FEsolving.f90 $(SUFFIX)



basics.a: math.o
	$(ARCHIVE_COMMAND) rc basics.a math.o debug.o numerics.o IO.o DAMASK_spectral_interface.o prec.o

math.o: math.f90 debug.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) math.f90 $(SUFFIX)
  
debug.o: debug.f90 numerics.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) debug.f90 $(SUFFIX)

numerics.o: numerics.f90 IO.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) numerics.f90 $(SUFFIX)

IO.o: IO.f90 DAMASK_spectral_interface.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) IO.f90 $(SUFFIX)
  
DAMASK_spectral_interface.o: DAMASK_spectral_interface.f90 prec.o
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) DAMASK_spectral_interface.f90 $(SUFFIX)

prec.o: prec.f90
	$(PREFIX)	$(COMPILERNAME) $(COMPILE) prec.f90 $(SUFFIX)


tidy: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a

clean: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	rm -rf *.exe

