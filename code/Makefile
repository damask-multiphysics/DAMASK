SHELL = /bin/sh
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

#auto values will be set by setup_code.py
FFTWROOT :=$(DAMASK_ROOT)/lib/fftw
#IKMLROOT := 
ACMLROOT :=/opt/acml4.4.0
#LAPACKROOT :=/usr

F90 ?=ifort
COMPILERNAME ?= $(F90)

ifeq "$(FASTBUILD)" "YES"
OPENMP :=OFF
OPTIMIZATION :=OFF
endif

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

ifneq "$(FASTBUILD)" "YES"
STANDARD_CHECK_ifort ?=-stand f08 -standard-semantics -warn stderrors
STANDARD_CHECK_gfortran ?=-std=f2008
endif

OPTIMIZATION_OFF_ifort           :=-O0 -no-ip
OPTIMIZATION_OFF_gfortran        :=-O0
OPTIMIZATION_DEFENSIVE_ifort     :=-O2
OPTIMIZATION_DEFENSIVE_gfortran  :=-O2
OPTIMIZATION_AGGRESSIVE_ifort    :=-O3 $(PORTABLE_SWITCH) -ipo -static -no-prec-div -fp-model fast=2 
OPTIMIZATION_AGGRESSIVE_gfortran :=-O3 $(PORTABLE_SWITCH) -ffast-math -funroll-loops -ftree-vectorize 


COMPILE_OPTIONS_ifort    :=-fpp\
                           -assume byterecl
ifneq "$(FASTBUILD)" "YES"
COMPILE_OPTIONS_ifort    :=$(COMPILE_OPTIONS_ifort)\
                           -implicitnone\
                           -diag-enable sc3\
                           -diag-disable 5268\
                           -warn declarations\
                           -warn general\
                           -warn usage\
                           -warn interfaces\
                           -warn ignore_loc\
                           -warn alignments\
                           -warn unused\
                           -warn errors
endif

#-fpp:                     preprocessor
#-fimplicit-none:          assume "implicit-none" even if not present in source
#-diag-disable:            disables warnings, where                           
#   warning ID 5268:       the text exceeds right hand column allowed on the line (we have only comments there)
#-warn:                    enables warnings, where
#   declarations:            any undeclared names
#   general:                 warning messages and informational messages are issued by the compiler
#   usage:                   questionable programming practices
#   interfaces:              checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks
#   ignore_loc:              %LOC is stripped from an actual argument
#   alignments:              data that is not naturally aligned
#   unused:                  declared variables that are never used
#   errors:                  warnings are changed to errors
#   stderrors:               warnings about Fortran standard violations are changed to errors (STANDARD_CHECK)
#
###################################################################################################
#MORE OPTIONS FOR DEBUGGING DURING COMPILING
#-warn:                    enables warnings, where
#   truncated_source:        Determines whether warnings occur when source exceeds the maximum column width in fixed-format files. (too many warnings because we have comments beyond character 132)
#   uncalled:                Determines whether warnings occur when a statement function is never called
#   all:
#
#OPTIONS FOR DEGUBBING DURING RUNTIME
# information on http://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors/
#-g:                       Generate symbolic debugging information in the object file
#-traceback:               Generate extra information in the object file to provide source file traceback information when a severe error occurs at run time.
#-gen-interfaces:          Generate an interface block for each routine. http://software.intel.com/en-us/blogs/2012/01/05/doctor-fortran-gets-explicit-again/  
#-fp-stack-check:          Generate extra code after every function call to ensure that the floating-point (FP) stack is in the expected state.
#-check:                   checks at runtime, where  
#   bounds:                  check if an array index is too small (<1) or too large!
#   arg_temp_created:        will cause a lot of warnings because we create a bunch of temporary arrays
#   format:                  Checking for the data type of an item being formatted for output.
#   output_conversion:       Checking for the fit of data items within a designated format descriptor field.
#   pointers:                Checking for certain disassociated or uninitialized pointers or unallocated allocatable objects.
#   uninit:                  Checking for uninitialized variables.
#-heap-arrays:            should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
#
#OPTIONS FOR TYPE DEBUGGING
#-real-size 32:           set precision to one of those 32/64/128 (= 4/8/16 bytes) for standard real (=8 for pReal)
#-integer-size 16:        set precision to one of those 16/32/64 (= 2/4/8 bytes) for standard integer (=4 for pInt)
###################################################################################################


COMPILE_OPTIONS_gfortran :=-xf95-cpp-input\
                           -fno-range-check
ifneq "$(FASTBUILD)" "YES"
COMPILE_OPTIONS_gfortran :=$(COMPILE_OPTIONS_gfortran)\
                           -xf95-cpp-input\
                           -ffree-line-length-132\
                           -fno-range-check\
                           -fimplicit-none\
                           -fall-intrinsics\
                           -pedantic\
                           -Warray-bounds\
                           -Wampersand\
                           -Wno-tabs\
                           -Wcharacter-truncation\
                           -Wintrinsic-shadow\
                           -Waliasing\
                           -Wconversion\
                           -Wsurprising\
                           -Wunderflow\
                           -Wswitch\
                           -Wstrict-overflow\
                           -Wattributes\
                           -Wunsafe-loop-optimizations\
                           -Wunused\
                           -Wall\
                           -Wextra\
                           -Wintrinsics-std
endif

#-xf95-cpp-input:              preprocessor
#-ffree-line-length-132:       restrict line length to the standard 132 characters
#-fno-range-check:             disables checking if result can be represented by variable. Needs to be set to enable DAMASK_NaN
#-fimplicit-none:              assume "implicit-none" even if not present in source
#-fall-intrinsics:        
#-pedantic:                    more strict on standard, enables some of the warnings below
#-Warray-bounds:               checks if array reference is out of bounds at compile time. use -fcheck-bounds to also check during runtime
#-Wampersand:                  checks if a character expression is continued proberly by an ampersand at the end of the line and at the beginning of the new line
#-Wno-tabs:                    do not allow tabs in source
#-Wcharacter-truncation:       warn if character expressions (strings) are truncated
#-Wintrinsic-shadow:           warn if a user-defined procedure or module procedure has the same name as an intrinsic
#-Waliasing:                   warn about possible aliasing of dummy arguments. Specifically, it warns if the same actual argument is associated with a dummy argument with "INTENT(IN)" and a dummy argument with "INTENT(OUT)" in a call with an explicit interface. 
#-Wconversion:                 warn about implicit conversions between different type
#-Wsurprising:                 warn when "suspicious" code constructs are encountered. While technically legal these usually indicate that an error has been made. 
#-Wunderflow:                  produce a warning when numerical constant expressions are encountered, which yield an UNDERFLOW during compilation
#-Wswitch:                     warn whenever a "switch" statement has an index of enumerated type and lacks a "case" for one or more of the named codes of that enumeration.  (The presence of a "default" label prevents this warning.) "case" labels outside the enumeration range also provokewarnings when this option is used (even if there is a "default" label)
#-Wstrict-overflow:
#-Wattributes:                 warn about inappropriate attribute usage                              
#-Wunsafe-loop-optimizations:  warn if the loop cannot be optimized due to nontrivial assumptions.
#-Wunused:   
#        -value:
#        -parameter:           find usused variables with "parameter" attribute
#-Wextra:
#-Wintrinsics-std:             warnings because of "flush" is not longer in the standard, but still an intrinsic fuction of the compilers:                        
###################################################################################################
#OPTIONS FOR GFORTRAN 4.6
#-Wsuggest-attribute=const:
#-Wsuggest-attribute=noreturn:
#-Wsuggest-attribute=pure:  
#-Wreal-q-constant:            Warn about real-literal-constants with 'q'  exponent-letter
#
#MORE OPTIONS FOR DEBUGGING DURING COMPILING                      
#-Wline-truncation:           too many warnings because we have comments beyond character 132
#-Warray-temporarieswarnings: because we have many temporary arrays (performance issue?):                     
#-Wimplicit-interface:
#-pedantic-errors:
#-fmodule-private:
#
#OPTIONS FOR DEGUBBING DURING RUNTIME
#-fcheck-bounds:              check if an array index is too small (<1) or too large!
#
#OPTIONS FOR TYPE DEBUGGING
#-fdefault-real-8:           set precision to 8 bytes for standard real (=8 for pReal). Will set size of double to 16 bytes as long as -fdefault-double-8 is not set
#-fdefault-integer-8:        set precision to 8 bytes for standard integer (=4 for pInt)
##################################################################################################

COMPILE         =$(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(STANDARD_CHECK_$(F90)) $(OPTIMIZATION_$(OPTI)_$(F90))
COMPILE_MAXOPTI =$(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(STANDARD_CHECK_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90))
###################################################################################################
COMPILED_FILES = prec.o DAMASK_spectral_interface.o IO.o numerics.o debug.o math.o \
                 FEsolving.o mesh.o material.o lattice.o \
                 constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o \
                 constitutive_titanmod.o constitutive_nonlocal.o constitutive.o \
                 homogenization_RGC.o homogenization_isostrain.o homogenization.o CPFEM.o crystallite.o

ifneq "$(SOLVER)" "AL"
DAMASK_spectral.exe: DAMASK_spectral.o 
	$(PREFIX) $(COMPILERNAME) $(OPENMP_FLAG_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90)) $(STANDARD_CHECK_$(F90)) \
  -o DAMASK_spectral.exe DAMASK_spectral.o \
  $(COMPILED_FILES) $(LIB_DIRS) $(LIBRARIES) $(SUFFIX)

DAMASK_spectral.o: DAMASK_spectral.f90 CPFEM.o
	$(PREFIX) $(COMPILERNAME) $(COMPILE_MAXOPTI) -c DAMASK_spectral.f90 $(SUFFIX)
else
DAMASK_spectral_AL.exe: DAMASK_spectral_AL.o 
	$(PREFIX) $(COMPILERNAME) $(OPENMP_FLAG_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90)) $(STANDARD_CHECK_$(F90)) \
  -o DAMASK_spectral_AL.exe DAMASK_spectral_AL.o \
  $(COMPILED_FILES) $(LIB_DIRS) $(LIBRARIES) $(SUFFIX)

DAMASK_spectral_AL.o: DAMASK_spectral_AL.f90 CPFEM.o
	$(PREFIX) $(COMPILERNAME) $(COMPILE_MAXOPTI) -c DAMASK_spectral_AL.f90 $(SUFFIX)
endif
CPFEM.o:                      CPFEM.f90\
                              homogenization.o

homogenization.o:             homogenization.f90\
                              homogenization_RGC.o \
                              homogenization_isostrain.o

homogenization_RGC.o:         homogenization_RGC.f90 \
                              crystallite.o

homogenization_isostrain.o:   homogenization_isostrain.f90 \
                              crystallite.o

crystallite.o:                crystallite.f90 \
                              constitutive.o
                            
constitutive.o:               constitutive.f90 \
                              constitutive_nonlocal.o \
                              constitutive_titanmod.o \
                              constitutive_dislotwin.o \
                              constitutive_phenopowerlaw.o \
                              constitutive_j2.o

constitutive_nonlocal.o:      constitutive_nonlocal.f90 \
                              lattice.o

constitutive_titanmod.o:      constitutive_titanmod.f90 \
                              lattice.o

constitutive_dislotwin.o:     constitutive_dislotwin.f90 \
                              lattice.o

constitutive_phenopowerlaw.o: constitutive_phenopowerlaw.f90 \
                              lattice.o

constitutive_j2.o:            constitutive_j2.f90 \
                              lattice.o

lattice.o:                    lattice.f90 \
                              material.o

material.o:                   material.f90 \
                              mesh.o

mesh.o:                       mesh.f90 \
                              FEsolving.o \
                              math.o

FEsolving.o:                  FEsolving.f90 \
                              debug.o

math.o:                       math.f90 \
                              debug.o

debug.o:                      debug.f90 \
                              numerics.o

numerics.o:                   numerics.f90 \
                              IO.o
                              
IO.o:                         IO.f90 \
                              DAMASK_spectral_interface.o

DAMASK_spectral_interface.o:  DAMASK_spectral_interface.f90 \
                              prec.o

prec.o:                       prec.f90

%.o : %.f90
	$(PREFIX) $(COMPILERNAME) $(COMPILE) -c $< $(SUFFIX)

tidy: 
	rm -rf *.o
	rm -rf *.mod

clean: 
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.exe
