###################################################################################################
# IntelLLVM Compiler
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 19)

if(OPENMP)
  set(OPENMP_FLAGS "-fiopenmp")
endif()

if(OPTIMIZATION STREQUAL "OFF" OR OPTIMIZATION STREQUAL "DEBUG")
  set(OPTIMIZATION_FLAGS "-O0")
elseif(OPTIMIZATION STREQUAL "DEFENSIVE")
  set(OPTIMIZATION_FLAGS "-O2")
elseif(OPTIMIZATION STREQUAL "AGGRESSIVE")
  set(OPTIMIZATION_FLAGS "-O3 -fp-model strict -xHost -align array64byte") # -ipo/-flto give linker error
endif()

# set(STANDARD_CHECK "-stand f23 -standard-semantics -assume nostd_mod_proc_name")
# -assume std_mod_proc_name (included in -standard-semantics) causes problems if other modules
#  (PETSc, HDF5) are not compiled with this option
#  https://community.intel.com/t5/Intel-Fortran-Compiler/building-and-linking-against-a-Fortran-dynamic-library-on-Linux/td-p/1178514
# -standard-semantics causes a bizzare increase in runtime on matesting (Intel(R) Xeon(R) CPU X5670)
#  so use only -fpscomp logicals.
#  https://fortran-lang.discourse.group/t/performance-drop-when-using-intel-oneapi-with-standard-sematics-option/9944/5
set(STANDARD_CHECK "-stand f23 -fpscomp logicals")

# Link against shared Intel libraries instead of static ones:
set(LINKER_FLAGS   "${LINKER_FLAGS} -shared-intel")
# enforce use of ifx for MPI wrapper:
set(LINKER_FLAGS   "${LINKER_FLAGS} -fc=ifx")

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------
# disable flush underflow to zero, will be set if -O[1,2,3]:
set(COMPILE_FLAGS "${COMPILE_FLAGS} -no-ftz")

# disable warnings ...
set(COMPILE_FLAGS "${COMPILE_FLAGS} -diag-disable")
#   ... the text exceeds right hand column allowed on the line (enforced by pre-receive hook)
set(COMPILE_FLAGS "${COMPILE_FLAGS} 5268")
#   ... about deprecated forall (has nice syntax and most likely a performance advantage)
set(COMPILE_FLAGS "${COMPILE_FLAGS},7624")

# enable warnings ...
set(COMPILE_FLAGS "${COMPILE_FLAGS} -warn")
#   ... any undeclared names (alternative name: -implicitnone)
set(COMPILE_FLAGS "${COMPILE_FLAGS} declarations")
#   ... warning messages and informational messages are issued by the compiler
set(COMPILE_FLAGS "${COMPILE_FLAGS},general")
#   ... questionable programming practices
set(COMPILE_FLAGS "${COMPILE_FLAGS},usage")
#   ... checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks
set(COMPILE_FLAGS "${COMPILE_FLAGS},interfaces")
#   ... %LOC is stripped from an actual argument
set(COMPILE_FLAGS "${COMPILE_FLAGS},ignore_loc")
#   ... data that is not naturally aligned
set(COMPILE_FLAGS "${COMPILE_FLAGS},alignments")
#   ... declared variables that are never used
set(COMPILE_FLAGS "${COMPILE_FLAGS},unused")

# Additional options
#  -warn:                    enables warnings, where
#     truncated_source:        Determines whether warnings occur when source exceeds the maximum column width in fixed-format files.
#                               (too many warnings because we have comments beyond character 132)
#     uncalled:                Determines whether warnings occur when a statement function is never called
#     all:
#  -name as_is: case sensitive Fortran!

#------------------------------------------------------------------------------------------------
# Runtime debugging
#------------------------------------------------------------------------------------------------
# Generate symbolic debugging information in the object file
set(DEBUG_FLAGS "${DEBUG_FLAGS} -g")

# Generate extra information in the object file to provide source file traceback information when a severe error occurs at run time
set(DEBUG_FLAGS "${DEBUG_FLAGS} -traceback")

# Generate an interface block for each routine. http://software.intel.com/en-us/blogs/2012/01/05/doctor-fortran-gets-explicit-again/
set(DEBUG_FLAGS "${DEBUG_FLAGS} -gen-interfaces")

# Generate extra code after every function call to ensure that the floating-point (FP) stack is in the expected state
# set(DEBUG_FLAGS "${DEBUG_FLAGS} -fp-stack-check") not available on ifx 2025.0.4

# Trap uninitalized variables
set(DEBUG_FLAGS "${DEBUG_FLAGS} -fp-model strict")

# Check at runtime ...
set(DEBUG_FLAGS "${DEBUG_FLAGS} -check" )
#   ... if an array index is too small (<1) or too large!
set(DEBUG_FLAGS "${DEBUG_FLAGS} bounds")
#   ... for the data type of an item being formatted for output.
set(DEBUG_FLAGS "${DEBUG_FLAGS},format")
#   ... for the fit of data items within a designated format descriptor field.
set(DEBUG_FLAGS "${DEBUG_FLAGS},output_conversion")
#   ... for certain disassociated or uninitialized pointers or unallocated allocatable objects.
set(DEBUG_FLAGS "${DEBUG_FLAGS},pointers")
#   ... for uninitialized variables.
set(DEBUG_FLAGS "${DEBUG_FLAGS},nouninit") # https://fortran-lang.discourse.group/t/issue-with-stdlib-and-intel-oneapi-fortran-compiler-ifx-2024-0/7049/4
#   ... capture all floating-point exceptions, need to overwrite -no-ftz
set(DEBUG_FLAGS "${DEBUG_FLAGS} -fpe-all=0 -ftz")

# Initialize logical to false, integer to -huge, float+complex to signaling NaN
set(DEBUG_FLAGS "${DEBUG_FLAGS} -init=arrays,zero,minus_huge,snan")

# disable due to compiler bug https://community.intel.com/t5/Intel-Fortran-Compiler/false-positive-stand-f18-and-IEEE-SELECTED-REAL-KIND/m-p/1227336
# enables warnings ...
#set(DEBUG_FLAGS "${DEBUG_FLAGS} -warn")
#   ... warnings are changed to errors
#set(DEBUG_FLAGS "${DEBUG_FLAGS} errors")
#   ... warnings about Fortran standard violations are changed to errors
#set(DEBUG_FLAGS "${DEBUG_FLAGS},stderrors")

# generate debug information for parameters
set(DEBUG_FLAGS "${DEBUG_FLAGS} -debug-parameters all")

# generate complete debugging information
# Additional options
# -heap-arrays:            Should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
# -check:                  Checks at runtime, where
#    arg_temp_created:       will cause a lot of warnings because we create a bunch of temporary arrays (performance?)
#    stack:
set(DEBUG_FLAGS "${DEBUG_FLAGS} -debug all")
