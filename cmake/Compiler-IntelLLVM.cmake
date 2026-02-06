# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# IntelLLVM Compiler
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 19)

set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEBUG      "${_OPTIMIZATION_OFF}")
set(_OPTIMIZATION_DEFENSIVE  "-O2")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -fp-model strict -xHost -align array64byte") # -ipo/-flto break link

if(DEFINED _OPTIMIZATION_${OPTIMIZATION})
  set(OPTIMIZATION_FLAGS "${_OPTIMIZATION_${OPTIMIZATION}}")
else()
  message(FATAL_ERROR "Unknown OPTIMIZATION level: ${OPTIMIZATION}")
endif()

if(OPENMP)
  set(OPENMP_FLAGS "-fiopenmp")
endif()

set(STANDARD_CHECK "-stand f23 -fpscomp logicals -assume noold_unit_star")
# -standard-semantics causes a bizarre increase in runtime on matesting (Intel(R) Xeon(R) CPU
# X5670), o only enforce a logical representation that is compatible with C and enable to close the
# output unit.
# https://fortran-lang.discourse.group/t/performance-drop-when-using-intel-oneapi-with-standard-sematics-option/9944/5
# https://community.intel.com/t5/Intel-Fortran-Compiler/Redirecting-File-STDOUT-and-STDERR-in-Intel-Fortran-Without/m-p/1550590

string(APPEND LINKER_FLAGS
  " -shared-intel"   # link against shared Intel runtime libraries
  " -fc=ifx"         # force MPI wrapper to use ifx
)

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

string(APPEND COMPILE_FLAGS
  " -no-ftz"                     # disable flush-underflow-to-zero (often enabled by -O[1-3])

  " -diag-disable"               # disable selected diagnostics...
  " 5268,7624"                   # ... line too long, deprecated FORALL

  " -warn"                       # enable warnings...
  " declarations"                # ... undeclared names (aka -implicitnone)
  ",general"                     # ... compiler warnings + informational messages
  ",usage"                       # ... questionable programming practices
  ",interfaces"                  # ... check calls vs external interface blocks
  ",ignore_loc"                  # ... warn when %LOC is stripped from actual argument
  ",alignments"                  # ... data not naturally aligned
  ",unused"                      # ... declared but unused variables
)

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

string(APPEND DEBUG_FLAGS
  " -g"                              # generate symbolic debug information
  " -traceback"                      # source file traceback on severe runtime errors
  " -gen-interfaces"                 # generate interface blocks for each routine

  " -fp-model strict"                # trap use of uninitialized FP values

  " -check"                          # enable runtime checks...
  " bounds"                          # ... array bounds
  ",format"                          # ... data type vs format descriptor
  ",output_conversion"               # ... data fit in format field
  ",pointers"                        # ... disassociated/uninitialized pointers
  ",nouninit"                        # ... uninitialized variables (ifx workaround)

  " -fpe-all=0 -ftz"                 # capture all FP exceptions, override -no-ftz

  " -init=arrays,zero,minus_huge,snan" # initialize arrays and scalars to sentinel values

  " -debug-parameters all"           # generate debug info for parameters
  " -debug all"                      # generate complete debugging information
)

# Generate extra code after every function call to ensure that the floating-point (FP) stack is in the expected state
# set(DEBUG_FLAGS "${DEBUG_FLAGS} -fp-stack-check") not available on ifx 2025.0.4
# disable due to compiler bug https://community.intel.com/t5/Intel-Fortran-Compiler/false-positive-stand-f18-and-IEEE-SELECTED-REAL-KIND/m-p/1227336
# enables warnings ...
#set(DEBUG_FLAGS "${DEBUG_FLAGS} -warn")
#   ... warnings are changed to errors
#set(DEBUG_FLAGS "${DEBUG_FLAGS} errors")
#   ... warnings about Fortran standard violations are changed to errors
#set(DEBUG_FLAGS "${DEBUG_FLAGS},stderrors")
# generate complete debugging information
# Additional options
# -heap-arrays:            Should not be done for OpenMP, but set "ulimit -s unlimited" on shell. Probably it helps also to unlimit other limits
# -check:                  Checks at runtime, where
#    arg_temp_created:       will cause a lot of warnings because we create a bunch of temporary arrays (performance?)
#    stack:
