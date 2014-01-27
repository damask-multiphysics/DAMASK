SHELL = /bin/sh
########################################################################################
# Makefile to compile the Material subroutine for BVP solution using spectral method
########################################################################################
# Be sure to remove all files compiled with different options by using "make clean"
# 
# Uses OpenMP to parallelize the material subroutines
# (set number of threads with "export DAMASK_NUM_THREADS=n" to n)
#
# Install fftw3 (v3.3 is tested):
# + execute
#   ./configure --enable-threads --enable-sse2 --enable-shared
#   make
#   make install
# + set FFTW_ROOT in your envinronment
#   We essentially look for two library files "lib/libfftw3_threads" and "lib/libfftw3", 
# Uses linux threads to parallelize fftw3
#
# Instead of the AMD Core Math Library or the Intel Kernel Math Library 
# a standard "lib(64)/liblapack.a/dylib/etc." can be used.
########################################################################################
# OPTIONS = standard (alternative): meaning
#-------------------------------------------------------------
# F90            = ifort (gfortran): compiler, choose Intel or GNU
# COMPILERNAME   = overwrite name of Compiler, e.g. using mpich-g90 instead of ifort
# PORTABLE       = TRUE (FALSE): decision, if executable is optimized for the machine on which it was built.
# OPTIMIZATION   = DEFENSIVE (OFF,AGGRESSIVE,ULTRA): Optimization mode: O2, O0, O3 + further options for most files, O3 + further options for all files
# OPENMP         = TRUE (FALSE): OpenMP multiprocessor support
# FFTW_ROOT       = root path
# IMKL_ROOT       = root path
# ACML_ROOT       = root path
# LAPACK_ROOT     = root path
# PREFIX         = arbitrary prefix (before compilername)
# OPTION         = arbitrary option (just before file to compile)
# SUFFIX         = arbitrary suffix (after file to compile)
# STANDARD_CHECK = checking for Fortran 2008, compiler dependend
########################################################################################
ifeq ($(wildcard "${HOME}/.damask/damask.conf"),)
	include ${HOME}/.damask/damask.conf
else
	include /etc/damask.conf
endif  


COMPILERNAME ?= $(F90)

INCLUDE_DIRS :=-I../lib
LIBRARIES :=-lfftw3
LIB_DIRS  :=-L$$FFTW_ROOT/lib

ifeq "$(FASTBUILD)" "YES"
OPENMP := OFF
OPTIMIZATION := OFF
else
OPENMP ?= ON
OPTIMIZATION ?= DEFENSIVE
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
IMKL_ARCH =_intel_thread
LIBRARIES +=-lfftw3_threads -lpthread
else
IMKL_ARCH =_sequential
endif

ifneq  "x$(IMKL_ROOT)" "x"
LIB_DIRS  :=-L$$IMKL_ROOT/lib/intel64
INCLUDE_DIRS +=-I$(IMKL_ROOT)/include
LIBRARIES +=-lmkl_intel_lp64 -lmkl_core -lmkl$(ACML_ARCH) -lm -lpthread -liomp5
else
ifneq  "x$(ACML_ROOT)" "x"
LIB_DIRS  +=-L$$ACML_ROOT/$(F90)64$(ACML_ARCH)/lib
LIBRARIES +=-lacml$(ACML_ARCH)
else
ifneq  "x$(LAPACK_ROOT)" "x"
LIB_DIRS  +=-L$$LAPACK_ROOT/lib64 -L$$LAPACK_ROOT/lib
LIBRARIES +=-llapack
endif
endif
endif




#hdf5
ifeq "$(HDF5)" "ON"
LIBRARIES +=-lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
LIB_DIRS  +=-L$$HDF5ROOT/lib
INCLUDE_DIRS +=-I$(HDF5ROOT)/include -DHDF
endif

ifdef PETSC_DIR
include ${PETSC_DIR}/conf/variables
INCLUDE_DIRS +=${PETSC_FC_INCLUDES} -DPETSc
LIBRARIES    +=${PETSC_WITH_EXTERNAL_LIB}
endif

ifdef STANDARD_CHECK
STANDARD_CHECK_ifort =$(STANDARD_CHECK)
STANDARD_CHECK_gfortran =$(STANDARD_CHECK)
endif

ifneq "$(FASTBUILD)" "YES"
STANDARD_CHECK_ifort ?=-stand f08 -standard-semantics -warn stderrors
STANDARD_CHECK_gfortran ?=-std=f2008 -pedantic-errors
endif
#-std=f2008ts:      for newer gfortran
#-pedantic:                    more strict on standard, enables some warnings
#   -pedantic-errors:          like pedantic, but errors instead of warnings
OPTIMIZATION_OFF_ifort           :=-O0 -no-ip
OPTIMIZATION_OFF_gfortran        :=-O0
OPTIMIZATION_DEFENSIVE_ifort     :=-O2
OPTIMIZATION_DEFENSIVE_gfortran  :=-O2
OPTIMIZATION_AGGRESSIVE_ifort    :=-O3 $(PORTABLE_SWITCH) -no-prec-div -fp-model fast=2 -ipo #-static causes trouble at the moment
OPTIMIZATION_AGGRESSIVE_gfortran :=-O3 $(PORTABLE_SWITCH) -ffast-math -funroll-loops -ftree-vectorize 


COMPILE_OPTIONS_ifort    :=-fpp\
                           -ftz\
                           -assume byterecl
ifneq "$(FASTBUILD)" "YES"
COMPILE_OPTIONS_ifort    +=-diag-enable sc3\
                           -diag-disable 5268\
                           -warn declarations\
                           -warn general\
                           -warn usage\
                           -warn interfaces\
                           -warn ignore_loc\
                           -warn alignments\
                           -warn unused
endif
###################################################################################################
#COMPILE SWITCHES
#-fpp:                     preprocessor
#-ftz:                     flush unterflow to zero, automatically set if O<0,1,2,3> >0
#-assume byterecl          record length is given in bytes (also set by -standard-semantics)
#-fimplicit-none:          assume "implicit-none" even if not present in source
#-diag-disable:            disables warnings, where                           
#   warning ID 5268:       the text exceeds right hand column allowed on the line (we have only comments there)
#-warn:                    enables warnings, where
#   declarations:            any undeclared names (alternative name: -implicitnone)
#   general:                 warning messages and informational messages are issued by the compiler
#   usage:                   questionable programming practices
#   interfaces:              checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks
#   ignore_loc:              %LOC is stripped from an actual argument
#   alignments:              data that is not naturally aligned
#   unused:                  declared variables that are never used
#   stderrors:               warnings about Fortran standard violations are changed to errors (STANDARD_CHECK)
#
###################################################################################################
#MORE OPTIONS FOR DEBUGGING DURING COMPILATION
#-warn:                    enables warnings, where
#   truncated_source:        Determines whether warnings occur when source exceeds the maximum column width in fixed-format files. (too many warnings because we have comments beyond character 132)
#   uncalled:                Determines whether warnings occur when a statement function is never called
#   all:
# -name as_is: case sensitive Fortran!

DEBUG_OPTIONS_ifort      :=-g\
                           -traceback\
                           -gen-interfaces\
                           -fp-stack-check\
                           -check bounds,format,output_conversion,pointers,uninit\
                           -fpe-all0\
                           -warn errors
###################################################################################################
#COMPILE SWITCHES FOR RUNTIME DEBUGGING
#-g:                       Generate symbolic debugging information in the object file
#-traceback:               Generate extra information in the object file to provide source file traceback information when a severe error occurs at run time.
#-gen-interfaces:          Generate an interface block for each routine. http://software.intel.com/en-us/blogs/2012/01/05/doctor-fortran-gets-explicit-again/  
#-fp-stack-check:          Generate extra code after every function call to ensure that the floating-point (FP) stack is in the expected state.
#-check:                   checks at runtime, where  
#   bounds:                  check if an array index is too small (<1) or too large!
#   format:                  Checking for the data type of an item being formatted for output.
#   output_conversion:       Checking for the fit of data items within a designated format descriptor field.
#   pointers:                Checking for certain disassociated or uninitialized pointers or unallocated allocatable objects.
#   uninit:                  Checking for uninitialized variables.
#-fpe-all0                capture all floating-point exceptions, sets -ftz automatically
#-warn:                    enables warnings, where
#   errors:                  warnings are changed to errors
# information on http://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors/
###################################################################################################
#MORE OPTIONS FOR RUNTIME DEBUGGING
#-heap-arrays:            should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
#-check:                   checks at runtime, where  
#   arg_temp_created:        will cause a lot of warnings because we create a bunch of temporary arrays (performance?)
#   stack:
COMPILE_OPTIONS_gfortran :=-xf95-cpp-input
ifneq "$(FASTBUILD)" "YES"
COMPILE_OPTIONS_gfortran +=-ffree-line-length-132\
                           -fimplicit-none\
                           -fmodule-private\
                           -Wall\
                           -Wextra\
                           -Wcharacter-truncation\
                           -Wunderflow\
                           -Wsuggest-attribute=pure\
                           -Wsuggest-attribute=noreturn\
                           -Wconversion-extra\
                           -Wimplicit-procedure
endif
###################################################################################################
#COMPILE SWITCHES
#-xf95-cpp-input:              preprocessor
#-ffree-line-length-132:       restrict line length to the standard 132 characters
#-fno-range-check:             disables checking if result can be represented by variable. Needs to be set to enable DAMASK_NaN
#-fimplicit-none:              assume "implicit-none" even if not present in source
#-fmodule-private:             assume "private" even if not present in source
#-Wcharacter-truncation:       warn if character expressions (strings) are truncated
#-Wunderflow:                  produce a warning when numerical constant expressions are encountered, which yield an UNDERFLOW during compilation
#-Wsuggest-attribute=pure:
#-Wsuggest-attribute=noreturn:
#-Wconversion-extra
#-Wimplicit-procedure
#-Wall:                        sets the following Fortran options:
#  -Waliasing:                   warn about possible aliasing of dummy arguments. Specifically, it warns if the same actual argument is associated with a dummy argument with "INTENT(IN)" and a dummy argument with "INTENT(OUT)" in a call with an explicit interface.
#  -Wampersand:                  checks if a character expression is continued proberly by an ampersand at the end of the line and at the beginning of the new line
#  -Warray-bounds:               checks if array reference is out of bounds at compile time. use -fcheck-bounds to also check during runtime
#  -Wconversion:                 warn about implicit conversions between different type
#  -Wsurprising:                 warn when "suspicious" code constructs are encountered. While technically legal these usually indicate that an error has been made.
#  -Wc-binding-type:
#  -Wintrinsics-std:             only standard intrisics are available, e.g. "call flush(6)" will cause an error
#  -Wno-tabs:                    do not allow tabs in source
#  -Wintrinsic-shadow:           warn if a user-defined procedure or module procedure has the same name as an intrinsic
#  -Wline-truncation:
#  -Wtarget-lifetime:
#  -Wreal-q-constant:            warn about real-literal-constants with 'q'  exponent-letter
#  -Wunused:                     a number of unused-xxx warnings
# these are general (non -Fortran options) implied by -Wall
#  -Waddress   
#  -Warray-bounds (only with -O2)  
#  -Wc++11-compat  
#  -Wchar-subscripts  
#  -Wcomment  
#  -Wformat    
#  -Wmaybe-uninitialized 
#  -Wnonnull  
#  -Wparentheses  
#  -Wpointer-sign  
#  -Wreorder   
#  -Wreturn-type  
#  -Wsequence-point   
#  -Wstrict-aliasing  
#  -Wstrict-overflow=1  
#  -Wswitch  
#  -Wtrigraphs  
#  -Wuninitialized  
#  -Wunknown-pragmas  
#  -Wunused-function  
#  -Wunused-label     
#  -Wunused-value     
#  -Wunused-variable  
#  -Wvolatile-register-var 
#-Wextra:                      sets the following Fortran options:
#  -Wunuses-parameter:
#  -Wcompare-reals:
# these are general (non -Fortran options) implied by -Wextra
#  -Wclobbered  
#  -Wempty-body  
#  -Wignored-qualifiers 
#  -Wmissing-field-initializers  
#  -Woverride-init  
#  -Wsign-compare  
#  -Wtype-limits  
#  -Wuninitialized  
#  -Wunused-but-set-parameter (only with -Wunused or -Wall)
#  -Wno-globals

###################################################################################################
#MORE OPTIONS FOR DEBUGGING DURING COMPILATION
#-Warray-temporarieswarnings:   because we have many temporary arrays (performance issue?):                     
#-Wimplicit-interface:          no interfaces for lapack routines
#-Wunsafe-loop-optimizations:   warn if the loop cannot be optimized due to nontrivial assumptions.
#-Wstrict-overflow:

DEBUG_OPTIONS_gfortran    :=-g\
                           -fbacktrace\
                           -fdump-core\
                           -fcheck=all\
                           -ffpe-trap=invalid,zero,overflow

###################################################################################################
#COMPILE SWITCHES FOR RUNTIME DEBUGGING
#-ffpe-trap=invalid,\         stop execution if floating point exception is detected (NaN is silent)
#           zero,\
#           overflow
#-fcheck=all:                   sets the following Fortran options:
#array-temps
#bounds
#do
#mem
#pointer
#recursion
###################################################################################################
#MORE OPTIONS FOR RUNTIME DEBUGGING
#-ffpe-trap=precision,\
#           denormal, \
#           underflow

ifeq "$(DEBUG)" "ON"
COMPILE_OPTIONS_$(F90) +=$(DEBUG_OPTIONS_$(F90))
endif
COMPILE_OPTIONS_$(F90) +=$(OPTIONS)
PRECISION_ifort :=-real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4
#-real-size 32:               set precision to one of those 32/64/128 (= 4/8/16 bytes) for standard real (=8 for pReal)
#-integer-size 16:            set precision to one of those 16/32/64 (= 2/4/8 bytes) for standard integer (=4 for pInt)
PRECISION_gfortran :=-fdefault-real-8 -fdefault-double-8 -DFLOAT=8 -DINT=4
#-fdefault-real-8:            set precision to 8 bytes for standard real (=8 for pReal). Will set size of double to 16 bytes as long as -fdefault-double-8 is not set
#-fdefault-double-8:          set precision to 8 bytes for double real, would be 16 bytes because -fdefault-real-8 is used
#-fdefault-integer-8:         Use it to set precision to 8 bytes for integer, don't use it for the standard case of pInt=4 (there is no -fdefault-integer-4)

###################################################################################################
COMPILE         =$(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(STANDARD_CHECK_$(F90)) $(OPTIMIZATION_$(OPTI)_$(F90)) $(INCLUDE_DIRS) $(PRECISION_$(F90)) -DSpectral
COMPILE_MAXOPTI =$(OPENMP_FLAG_$(F90)) $(COMPILE_OPTIONS_$(F90)) $(STANDARD_CHECK_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90)) $(INCLUDE_DIRS) $(PRECISION_$(F90)) -DSpectral
###################################################################################################
COMPILED_FILES = prec.o DAMASK_spectral_interface.o IO.o libs.o numerics.o debug.o math.o \
                 FEsolving.o mesh.o material.o lattice.o \
                 constitutive_dislotwin.o constitutive_j2.o constitutive_phenopowerlaw.o \
                 constitutive_titanmod.o constitutive_nonlocal.o constitutive_none.o constitutive.o crystallite.o \
                 homogenization_RGC.o homogenization_isostrain.o homogenization.o CPFEM.o \
                 DAMASK_spectral_utilities.o DAMASK_spectral_solverBasic.o 

ifdef PETSC_DIR
PETSC_FILES = DAMASK_spectral_solverAL.o \
                        DAMASK_spectral_solverBasicPETSc.o \
                        DAMASK_spectral_solverPolarisation.o
COMPILED_FILES += $(PETSC_FILES)
endif


DAMASK_spectral.exe: DAMASK_spectral_driver.o 
	$(PREFIX) $(COMPILERNAME) $(OPENMP_FLAG_$(F90)) $(OPTIMIZATION_$(MAXOPTI)_$(F90)) $(STANDARD_CHECK_$(F90)) \
  -o DAMASK_spectral.exe DAMASK_spectral_driver.o \
  $(COMPILED_FILES) $(LIB_DIRS) $(LIBRARIES) $(SUFFIX)

DAMASK_spectral_driver.o: DAMASK_spectral_driver.f90 DAMASK_spectral_solverBasic.o $(PETSC_FILES)
	$(PREFIX) $(COMPILERNAME) $(COMPILE_MAXOPTI) -c DAMASK_spectral_driver.f90 $(SUFFIX)

DAMASK_spectral_solverAL.o:    DAMASK_spectral_solverAL.f90\
                               DAMASK_spectral_utilities.o
  
DAMASK_spectral_solverBasic.o: DAMASK_spectral_solverBasic.f90\
                               DAMASK_spectral_utilities.o
  
DAMASK_spectral_solverPolarisation.o:    DAMASK_spectral_solverPolarisation.f90\
                               DAMASK_spectral_utilities.o
  
DAMASK_spectral_solverBasicPETSc.o: DAMASK_spectral_solverBasicPETSc.f90\
                               DAMASK_spectral_utilities.o

DAMASK_spectral_utilities.o:   DAMASK_spectral_utilities.f90\
                               CPFEM.o
CPFEM.o:                       CPFEM.f90\
                               homogenization.o

homogenization.o:              homogenization.f90\
                               homogenization_RGC.o \
                               homogenization_isostrain.o

homogenization_RGC.o:          homogenization_RGC.f90 \
                               crystallite.o

homogenization_isostrain.o:    homogenization_isostrain.f90 \
                               crystallite.o

crystallite.o:                 crystallite.f90 \
                               constitutive.o
                            
constitutive.o:                constitutive.f90 \
                               constitutive_nonlocal.o \
                               constitutive_titanmod.o \
                               constitutive_dislotwin.o \
                               constitutive_phenopowerlaw.o \
                               constitutive_j2.o \
                               constitutive_none.o

constitutive_nonlocal.o:       constitutive_nonlocal.f90 \
                               lattice.o

constitutive_titanmod.o:       constitutive_titanmod.f90 \
                               lattice.o

constitutive_dislotwin.o:      constitutive_dislotwin.f90 \
                               lattice.o

constitutive_phenopowerlaw.o:  constitutive_phenopowerlaw.f90 \
                               lattice.o

constitutive_j2.o:             constitutive_j2.f90 \
                               lattice.o

constitutive_none.o:           constitutive_none.f90 \
                               lattice.o

lattice.o:                     lattice.f90 \
                               material.o

material.o:                    material.f90 \
                               mesh.o

mesh.o:                        mesh.f90 \
                               FEsolving.o \
                               math.o

FEsolving.o:                   FEsolving.f90 \
                               debug.o

math.o:                        math.f90 \
                               debug.o

debug.o:                       debug.f90 \
                               numerics.o

numerics.o:                    numerics.f90 \
                               libs.o
                               
libs.o:                        libs.f90 \
                               IO.o

IO.o:                          IO.f90 \
                               DAMASK_spectral_interface.o

ifeq "$(F90)" "gfortran"
DAMASK_spectral_interface.o:   DAMASK_spectral_interface.f90 \
                               prec.o
	$(PREFIX) $(COMPILERNAME) $(COMPILE) -c -fall-intrinsics DAMASK_spectral_interface.f90 $(SUFFIX)
#-fall-intrinsics:  all intrinsic procedures (including the GNU-specific extensions) are accepted. -Wintrinsics-std will be ignored 
#                   and no user-defined procedure with the same name as any intrinsic will be called except when it is explicitly declared external
#      -->   allows the use of 'getcwd'
prec.o:                        prec.f90
	$(PREFIX) $(COMPILERNAME) $(COMPILE) -c -fno-range-check prec.f90 $(SUFFIX)
# fno-range-check:  Disable range checking on results of simplification of constant expressions during compilation
#      -->   allows the definition of DAMASK_NaN
else
DAMASK_spectral_interface.o:   DAMASK_spectral_interface.f90 \
                               prec.o
prec.o:                        prec.f90
	$(PREFIX) $(COMPILERNAME) $(COMPILE) -c prec.f90 $(SUFFIX)
endif

%.o : %.f90
	$(PREFIX) $(COMPILERNAME) $(COMPILE) -c $< $(SUFFIX)

.PHONY: tidy
tidy:
	@rm -rf *.o
	@rm -rf *.mod


.PHONY: clean
clean:
	@rm -rf *.o
	@rm -rf *.mod
	@rm -rf *.exe

