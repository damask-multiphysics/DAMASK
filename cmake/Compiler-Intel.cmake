###################################################################################################
# Intel Compiler
###################################################################################################
if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 18.0)
  message (FATAL_ERROR "Intel Compiler version: ${CMAKE_Fortran_COMPILER_VERSION} not supported")
endif ()

if (OPENMP)
  set (OPENMP_FLAGS "-qopenmp -parallel")
endif ()

if (OPTIMIZATION STREQUAL "OFF")
  set (OPTIMIZATION_FLAGS    "-O0 -no-ip")
elseif (OPTIMIZATION STREQUAL "DEFENSIVE")
  set (OPTIMIZATION_FLAGS    "-O2")
elseif (OPTIMIZATION STREQUAL "AGGRESSIVE")
  set (OPTIMIZATION_FLAGS    "-ipo -O3 -no-prec-div -fp-model fast=2 -xHost")
  # -fast = -ipo, -O3, -no-prec-div, -static, -fp-model fast=2, and -xHost"
endif ()

# -assume std_mod_proc_name (included in -standard-semantics) causes problems if other modules
# (PETSc, HDF5) are not compiled with this option (https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/62172)
set (STANDARD_CHECK "-stand f18 -standard-semantics -assume nostd_mod_proc_name")
set (LINKER_FLAGS   "${LINKER_FLAGS} -shared-intel")
# Link against shared Intel libraries instead of static ones

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
set (COMPILE_FLAGS "${COMPILE_FLAGS} -fpp")
# preprocessor

set (COMPILE_FLAGS "${COMPILE_FLAGS} -ftz")
# flush underflow to zero, automatically set if -O[1,2,3]

set (COMPILE_FLAGS "${COMPILE_FLAGS} -diag-disable")
# disables warnings ...
set (COMPILE_FLAGS "${COMPILE_FLAGS} 5268")
#   ... the text exceeds right hand column allowed on the line (we have only comments there)
set (COMPILE_FLAGS "${COMPILE_FLAGS},7624")
#   ... about deprecated forall (has nice syntax and most likely a performance advantage)

set (COMPILE_FLAGS "${COMPILE_FLAGS} -warn")
# enables warnings ...
set (COMPILE_FLAGS "${COMPILE_FLAGS} declarations")
#   ... any undeclared names (alternative name: -implicitnone)
set (COMPILE_FLAGS "${COMPILE_FLAGS},general")
#   ... warning messages and informational messages are issued by the compiler
set (COMPILE_FLAGS "${COMPILE_FLAGS},usage")
#   ... questionable programming practices
set (COMPILE_FLAGS "${COMPILE_FLAGS},interfaces")
#   ... checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks
set (COMPILE_FLAGS "${COMPILE_FLAGS},ignore_loc")
#   ... %LOC is stripped from an actual argument
set (COMPILE_FLAGS "${COMPILE_FLAGS},alignments")
#   ... data that is not naturally aligned
set (COMPILE_FLAGS "${COMPILE_FLAGS},unused")
#   ... declared variables that are never used

# Additional options
#  -warn:                    enables warnings, where
#     truncated_source:        Determines whether warnings occur when source exceeds the maximum column width in fixed-format files.
#                               (too many warnings because we have comments beyond character 132)
#     uncalled:                Determines whether warnings occur when a statement function is never called
#     all:
#  -name as_is: case sensitive Fortran!

#------------------------------------------------------------------------------------------------
# Runtime debugging
set (DEBUG_FLAGS "${DEBUG_FLAGS} -g")
# Generate symbolic debugging information in the object file

set (DEBUG_FLAGS "${DEBUG_FLAGS} -traceback")
# Generate extra information in the object file to provide source file traceback information when a severe error occurs at run time

set (DEBUG_FLAGS "${DEBUG_FLAGS} -gen-interfaces")
# Generate an interface block for each routine. http://software.intel.com/en-us/blogs/2012/01/05/doctor-fortran-gets-explicit-again/

set (DEBUG_FLAGS "${DEBUG_FLAGS} -fp-stack-check")
# Generate extra code after every function call to ensure that the floating-point (FP) stack is in the expected state

set (DEBUG_FLAGS "${DEBUG_FLAGS} -fp-model strict")
# Trap uninitalized variables

set (DEBUG_FLAGS "${DEBUG_FLAGS} -check" )
# Checks at runtime ...
set (DEBUG_FLAGS "${DEBUG_FLAGS} bounds")
#   ... if an array index is too small (<1) or too large!
set (DEBUG_FLAGS "${DEBUG_FLAGS},format")
#   ... for the data type of an item being formatted for output.
set (DEBUG_FLAGS "${DEBUG_FLAGS},output_conversion")
#   ... for the fit of data items within a designated format descriptor field.
set (DEBUG_FLAGS "${DEBUG_FLAGS},pointers")
#   ... for certain disassociated or uninitialized pointers or unallocated allocatable objects.
set (DEBUG_FLAGS "${DEBUG_FLAGS},uninit")
#   ... for uninitialized variables.
set (DEBUG_FLAGS "${DEBUG_FLAGS} -ftrapuv")
#   ... initializes stack local variables to an unusual value to aid error detection
set (DEBUG_FLAGS "${DEBUG_FLAGS} -fpe-all=0")
#   ... capture all floating-point exceptions, sets -ftz automatically

# disable due to compiler bug https://community.intel.com/t5/Intel-Fortran-Compiler/false-positive-stand-f18-and-IEEE-SELECTED-REAL-KIND/m-p/1227336
#set (DEBUG_FLAGS "${DEBUG_FLAGS} -warn")
# enables warnings ...
#set (DEBUG_FLAGS "${DEBUG_FLAGS} errors")
#   ... warnings are changed to errors
#set (DEBUG_FLAGS "${DEBUG_FLAGS},stderrors")
#   ... warnings about Fortran standard violations are changed to errors

set (DEBUG_FLAGS "${DEBUG_FLAGS} -debug-parameters all")
# generate debug information for parameters

# Additional options
# -heap-arrays:            Should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
# -check:                  Checks at runtime, where
#    arg_temp_created:       will cause a lot of warnings because we create a bunch of temporary arrays (performance?)
#    stack:

#------------------------------------------------------------------------------------------------
#  precision settings
set (PRECISION_FLAGS "${PRECISION_FLAGS} -real-size 64")
# set precision for standard real to 32 | 64 | 128 (= 4 | 8 | 16 bytes, type pReal is always 8 bytes)
