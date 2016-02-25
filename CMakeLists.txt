########################################################################################
# CMAKE for build the Material subroutine for BVP solution
########################################################################################
# OPTIONS = standard (alternative): meaning
#-------------------------------------------------------------
# F90            = ifort (gfortran): compiler type, choose Intel or GNU
# COMPILERNAME   = name of the compiler executable (if not the same as the ype), e.g. using mpich-g90 instead of ifort
# PORTABLE       = TRUE (FALSE): decision, if executable is optimized for the machine on which it was built.
# OPTIMIZATION   = DEFENSIVE (OFF,AGGRESSIVE,ULTRA): Optimization mode: O2, O0, O3 + further options for most files, O3 + further options for all files
# OPENMP         = TRUE (FALSE): OpenMP multiprocessor support
# PREFIX         = arbitrary prefix (before compilername)
# OPTION         = arbitrary option (just before file to compile)
# SUFFIX         = arbitrary suffix (after file to compile)
# STANDARD_CHECK = checking for Fortran 2008, compiler dependend
########################################################################################
cmake_minimum_required (VERSION 3.1.0)
project (DAMASK Fortran)

# The version number.
set (DAMASK_VERSION_MAJOR 1)
set (DAMASK_VERSION_MINOR ${DAMASK_V})

##
# find_package() might not work on customized system.
# it is recommended for the user to specify the location
# of PETSc and HDF5 library
if (defined PETSC_DIR)
    set (PETSC ${PETSC_DIR})
else (defined PETSC_DIR)
    find_package(PETSC)
endif (defined PETSC_DIR)

if (defined HDF5_DIR)
    set (HDF5 ${HDF5_DIR})
else (defined HDF5_DIR)
    find_package(HDF5)
endif (defined HDF5_DIR)

##
# set system include directories
include_directories(
                    ${PETSC}/lib/petsc/conf/variables
                    ${PETSC}/lib/petsc/conf/rules
                    ../lib
                    ${HDF5}/include
                    )

##
# set compile and debug flags
if (OPENMP MATCHES 'ON')
    set (OPENMP_FLAG_ifort    "-openmp -openmp-report0 -parallel")
    set (OPENMP_FLAG_gfortran "-fopenmp")
else (OPENMP MATCHES 'ON')
    message("No openmp support specified")
    set (OPENMP_FLAG_ifort    "")
    set (OPENMP_FLAG_gfortran "")
endif (OPENMP MATCHES 'ON')

if (OPTIMIZATION MATCHES "OFF")
    set (OPTIMIZATION_ifort    "-O0 -no-ip")
    set (OPTIMIZATION_gfortran "-O0")
elseif (OPTIMIZATION MATCHES "DEFENSIVE")
    set (OPTIMIZATION_ifort    "-O2")
    set (OPTIMIZATION_gfortran "-O2")
elseif (OPTIMIZATION MATCHES "AGGRESSIVE")
    set (OPTIMIZATION_ifort    "-ipo -O3 -no-prec-div -fp-model fast=2 -xHost")  #-fast = -ipo, -O3, -no-prec-div, -static, -fp-model fast=2, and -xHost"
    set (OPTIMIZATION_gfortran "-O3 -ffast-math -funroll-loops -ftree-vectorize")
elseif (OPTIMIZATION MATCHES "ULTRA")
    set (OPTIMIZATION_ifort    "-ipo -O3 -no-prec-div -fp-model fast=2 -xHost")  #-fast = -ipo, -O3, -no-prec-div, -static, -fp-model fast=2, and -xHost"
    set (OPTIMIZATION_gfortran "-O3 -ffast-math -funroll-loops -ftree-vectorize")
else (OPTIMIZATION MATCHES "OFF")
    set (OPTIMIZATION_ifort    "-O2")
    set (OPTIMIZATION_gfortran "-O2")
endif (OPTIMIZATION MATCHES "OFF")

set (STANDARD_CHECK_ifort    "-stand f08 -standard-semantics")
set (STANDARD_CHECK_gfortran "-std=f2008ts -pedantic-errors" )

###################################################################################################
# COMPILE SWITCHES
#  -shared-intel:            Link against shared Intel libraries instead of static ones
#  -fpp:                     preprocessor
#  -ftz:                     flush unterflow to zero, automatically set if O<0,1,2,3> >0
#  -assume byterecl          record length is given in bytes (also set by -standard-semantics)
#          fpe_summary       print list of floating point exceptions occured during execution
#  -fimplicit-none:          assume "implicit-none" even if not present in source
#  -diag-disable:            disables warnings, where
#     warning ID 5268:       the text exceeds right hand column allowed on the line (we have only comments there)
#  -warn:                    enables warnings, where
#     declarations:            any undeclared names (alternative name: -implicitnone)
#     general:                 warning messages and informational messages are issued by the compiler
#     usage:                   questionable programming practices
#     interfaces:              checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks
#     ignore_loc:              %LOC is stripped from an actual argument
#     alignments:              data that is not naturally aligned
#     unused:                  declared variables that are never used
#     stderrors:               warnings about Fortran standard violations are changed to errors (STANDARD_CHECK)
#
###################################################################################################
# MORE OPTIONS FOR DEBUGGING DURING COMPILATION
#  -warn:                    enables warnings, where
#     truncated_source:        Determines whether warnings occur when source exceeds the maximum column width in fixed-format files. (too many warnings because we have comments beyond character 132)
#     uncalled:                Determines whether warnings occur when a statement function is never called
#     all:
#  -name as_is: case sensitive Fortran!
###################################################################################################
set (COMPILE_OPTIONS_ifort   "-DDAMASKVERSION=${DAMASK_V}"
                             "-fpp"
                             "-ftz"
                             "-assume byterecl,fpe_summary"
                             "-diag-disable 5268"
                             "-warn declarations"
                             "-warn general"
                             "-warn usage"
                             "-warn interfaces"
                             "-warn ignore_loc"
                             "-warn alignments"
                             "-warn unused"
    )

###################################################################################################
# COMPILE SWITCHES FOR RUNTIME DEBUGGING
# -g:                       Generate symbolic debugging information in the object file
# -traceback:               Generate extra information in the object file to provide source file traceback information when a severe error occurs at run time.
# -gen-interfaces:          Generate an interface block for each routine. http://software.intel.com/en-us/blogs/2012/01/05/doctor-fortran-gets-explicit-again/
# -fp-stack-check:          Generate extra code after every function call to ensure that the floating-point (FP) stack is in the expected state.
# -ftrapuv                  Trap uninitalized variables
# -check:                   checks at runtime, where
#    bounds:                  check if an array index is too small (<1) or too large!
#    format:                  Checking for the data type of an item being formatted for output.
#    output_conversion:       Checking for the fit of data items within a designated format descriptor field.
#    pointers:                Checking for certain disassociated or uninitialized pointers or unallocated allocatable objects.
#    uninit:                  Checking for uninitialized variables.
# -fpe-all0                capture all floating-point exceptions, sets -ftz automatically
# -warn:                    enables warnings, where
#    errors:                  warnings are changed to errors
#    stderrors:               warnings about Fortran standard violations are changed to errors
#  information on http://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors/
###################################################################################################
# MORE OPTIONS FOR RUNTIME DEBUGGING
# -heap-arrays:            should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
# -check:                   checks at runtime, where
#    arg_temp_created:        will cause a lot of warnings because we create a bunch of temporary arrays (performance?)
#    stack:
###################################################################################################
set (DEBUG_OPTIONS_ifort     "-g"
                             "-traceback"
                             "-gen-interfaces"
                             "-fp-stack-check"
                             "-fp-model strict"
                             "-check bounds,format,output_conversion,pointers,uninit"
                             "-ftrapuv"
                             "-fpe-all0"
                             "-warn errors"
                             "-warn stderrors"
                             "-debug-parameters all"
    )

set (LINK_OPTIONS_ifort      "-shared-intel")

###################################################################################################
# COMPILE SWITCHES
# -shared
# -Wl,-undefined,dynamic_lookup:ensure to link against dynamic libraries
# -xf95-cpp-input:              preprocessor
# -ffree-line-length-132:       restrict line length to the standard 132 characters
# -ffpe-summary:                print summary of floating point exeptions (‘invalid’, ‘zero’, ‘overflow’, ‘underflow’, ‘inexact’ and ‘denormal’)
# -fimplicit-none:              assume "implicit-none" even if not present in source
# -fmodule-private:             assume "private" even if not present in source
# -Wcharacter-truncation:       warn if character expressions (strings) are truncated
# -Wunderflow:                  produce a warning when numerical constant expressions are encountered, which yield an UNDERFLOW during compilation
# -Wsuggest-attribute=pure:
# -Wsuggest-attribute=noreturn:
# -Wconversion-extra
# -Wimplicit-procedure
# -Wall:                        sets the following Fortran options:
#   -Waliasing:                   warn about possible aliasing of dummy arguments. Specifically, it warns if the same actual argument is associated with a dummy argument with "INTENT(IN)" and a dummy argument with "INTENT(OUT)" in a call with an explicit interface.
#   -Wampersand:                  checks if a character expression is continued proberly by an ampersand at the end of the line and at the beginning of the new line
#   -Warray-bounds:               checks if array reference is out of bounds at compile time. use -fcheck-bounds to also check during runtime
#   -Wconversion:                 warn about implicit conversions between different type
#   -Wsurprising:                 warn when "suspicious" code constructs are encountered. While technically legal these usually indicate that an error has been made.
#   -Wc-binding-type:
#   -Wintrinsics-std:             only standard intrisics are available, e.g. "call flush(6)" will cause an error
#   -Wno-tabs:                    do not allow tabs in source
#   -Wintrinsic-shadow:           warn if a user-defined procedure or module procedure has the same name as an intrinsic
#   -Wline-truncation:
#   -Wtarget-lifetime:
#   -Wreal-q-constant:            warn about real-literal-constants with 'q'  exponent-letter
#   -Wunused:                     a number of unused-xxx warnings
#  these are general (non -Fortran options) implied by -Wall
#   -Waddress
#   -Warray-bounds (only with -O2)
#   -Wc++11-compat
#   -Wchar-subscripts
#   -Wcomment
#   -Wformat
#   -Wmaybe-uninitialized
#   -Wnonnull
#   -Wparentheses
#   -Wpointer-sign
#   -Wreorder
#   -Wreturn-type
#   -Wsequence-point
#   -Wstrict-aliasing
#   -Wstrict-overflow=1
#   -Wswitch
#   -Wtrigraphs
#   -Wuninitialized
#   -Wunknown-pragmas
#   -Wunused-function
#   -Wunused-label
#   -Wunused-value
#   -Wunused-variable
#   -Wvolatile-register-var
# -Wextra:                      sets the following Fortran options:
#   -Wunuses-parameter:
#   -Wcompare-reals:
#  these are general (non -Fortran options) implied by -Wextra
#   -Wclobbered
#   -Wempty-body
#   -Wignored-qualifiers
#   -Wmissing-field-initializers
#   -Woverride-init
#   -Wsign-compare
#   -Wtype-limits
#   -Wuninitialized
#   -Wunused-but-set-parameter (only with -Wunused or -Wall)
#   -Wno-globals
# -ffpe-summary=all            only for newer gfortran
###################################################################################################
# MORE OPTIONS FOR DEBUGGING DURING COMPILATION
# -Warray-temporarieswarnings:   because we have many temporary arrays (performance issue?):
# -Wimplicit-interface:          no interfaces for lapack routines
# -Wunsafe-loop-optimizations:   warn if the loop cannot be optimized due to nontrivial assumptions.
# -Wstrict-overflow:
###################################################################################################
set (COMPILE_OPTIONS_gfortran "-DDAMASKVERSION=${DAMASKVERSION}"
                              "-xf95-cpp-input"
                              "-ffree-line-length-132"
                              "-fimplicit-none"
                              "-fmodule-private"
                              "-Wall"
                              "-Wextra"
                              "-Wcharacter-truncation"
                              "-Wunderflow"
                              "-Wsuggest-attribute=pure"
                              "-Wsuggest-attribute=noreturn"
                              "-Wconversion-extra"
                              "-Wimplicit-procedure"
                              "-Wno-unused-parameter"
    )

###################################################################################################
# COMPILE SWITCHES FOR RUNTIME DEBUGGING
# -ffpe-trap=invalid,\         stop execution if floating point exception is detected (NaN is silent)
#            zero,\
#            overflow
# -fcheck=all:                   sets the following Fortran options:
# array-temps
# bounds
# do
# mem
# pointer
# recursion
###################################################################################################
# MORE OPTIONS FOR RUNTIME DEBUGGING
# -ffpe-trap=precision,\
#            denormal, \
#            underflow
###################################################################################################
set (DEBUG_OPTIONS_gfortran "-g"
                            "-fbacktrace"
                            "-fdump-core"
                            "-fcheck=all"
                            "-ffpe-trap=invalid,zero,overflow"
    )

set (LINK_OPTIONS_gfortran "-Wl,-undefined,dynamic_lookup")

# set FLAGS
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    # for RELEASE
    set (CMAKE_Fortran_FLAGS_RELEASE     ${COMPILE_OPTIONS_ifort})
    set (CMAKE_EXE_LINKER_FLAGS_RELEASE  ${LINK_OPTIONS_ifort}
                                         ${OPTIMIZATION_ifort}   )
    # for DEBUG
    set (CMAKE_Fortran_FLAGS_DEBUG       ${COMPILE_OPTIONS_ifort}
                                         ${DEBUG_OPTIONS_ifort}  )
    set (CMAKE_EXE_LINKER_FLAGS_DEBUG    ${LINK_OPTIONS_ifort}
                                         ${DEBUG_OPTIONS_ifort}  )
    #
elseif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # for RELEASE
    set (CMAKE_Fortran_FLAGS_RELEASE     ${COMPILE_OPTIONS_gfortran})
    set (CMAKE_EXE_LINKER_FLAGS_RELEASE  ${LINK_OPTIONS_gfortran}
                                         ${OPTIMIZATION_gfortran}   )
    # for DEBUG
    set (CMAKE_Fortran_FLAGS_DEBUG       ${COMPILE_OPTIONS_gfortran}
                                         ${DEBUG_OPTIONS_gfortran}  )
    set (CMAKE_EXE_LINKER_FLAGS_DEBUG    ${LINK_OPTIONS_gfortran}
                                         ${DEBUG_OPTIONS_gfortran}  )
    #
elseif (Fortran_COMPILER_NAME MATCHES "g77")
    message (FATAL_ERROR "Fortran 77 is not supported.")
else (Fortran_COMPILER_NAME MATCHES "ifort.*")
    message (FATAL_ERROR "Require Fortran90 from GNU or Intel.")
endif (Fortran_COMPILER_NAME MATCHES "ifort.*")

##
# ADD CODE(SOURCE) DIRECTORY
add_subdirectory(code)

##
# ADD TESTING CASES
add_test (SmokeTestRun
          DAMASK_spectral.exe -g test/test1.geom -l test/test.load)

# Enable Dashboard scripting
# include (CTest)
# set (CTEST_PROJECT_NAME "DAMASK")