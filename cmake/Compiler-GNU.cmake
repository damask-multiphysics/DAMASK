# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# GNU Compiler
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 11.1)

set(_OPTIMIZATION_DEBUG      "-Og")
set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEFENSIVE  "-O2 -mtune=native")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -march=native -funroll-loops -ftree-vectorize -flto")

if(DEFINED _OPTIMIZATION_${OPTIMIZATION})
  set(OPTIMIZATION_FLAGS "${_OPTIMIZATION_${OPTIMIZATION}}")
else()
  message(FATAL_ERROR "Unknown OPTIMIZATION level: ${OPTIMIZATION}")
endif()

if(OPENMP)
  set(OPENMP_FLAGS "-fopenmp")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 14)
  set(STANDARD_CHECK "-std=f2018 -pedantic-errors")
else()
  set(STANDARD_CHECK "-std=f2023 -pedantic-errors")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 12)
  add_compile_definitions(OLD_STYLE_C_TO_FORTRAN_STRING)
endif()

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

string(APPEND COMPILE_FLAGS
  " -fPIE"                          # position independent code
  " -ffree-line-length-none"        # PETSc macros exceed line-length limits
  " -fimplicit-none"                # assume implicit none if not present
  " -Wall"                          # enable common warnings
  " -Wextra"                        # enable extra warnings
  " -Wcharacter-truncation"         # warn on string truncation
  " -Wunderflow"                    # warn on compile-time underflow
  " -Wsuggest-attribute=pure"       # suggest PURE where applicable
  " -Wsuggest-attribute=noreturn"   # suggest NORETURN where applicable
  " -Wconversion-extra"             # extra conversion warnings
  " -Wimplicit-procedure"           # warn on implicit procedure calls
  " -Wunused-parameter"             # warn on unused parameters
  " -Wimplicit-interface"           # warn on missing explicit interfaces
  " -Wno-maybe-uninitialized"       # suppress false positives
  " -Wno-c-binding-type"            # suppress MPI_f08 warnings
  " -ffpe-summary=all"              # report FP exceptions summary
  " -fno-unsafe-math-optimizations" # required for IEEE semantics
  " -frounding-math"                # honor rounding mode
  " -fsignaling-nans"               # enable signaling NaNs
)

# set the following Fortran options:
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
# and set the general (non-Fortran options) options:
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

# set the following Fortran options:
#   -Wunuses-parameter:
#   -Wcompare-reals:
# and set the general (non-Fortran options) options:
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

# Additional options
# -Wunsafe-loop-optimizations:   warn if the loop cannot be optimized due to nontrivial assumptions

#------------------------------------------------------------------------------------------------
# Runtime debugging
#------------------------------------------------------------------------------------------------

string(APPEND DEBUG_FLAGS
  " -ffpe-trap=invalid,zero,overflow"  # stop on FP exceptions (NaN remains silent)
  # " -ffpe-trap=precision,denormal,underflow"  # optional, more aggressive traps

  " -g"                               # generate debug symbols
  " -Og"                              # optimize for debugging experience

  " -fbacktrace"                      # runtime backtrace on error
  " -fdump-core"                      # generate core dump
  " -fcheck=all"                      # runtime checks (bounds, pointers, etc.)

  " -fstack-protector-all"            # guard variables on all stack frames

  " -finit-real=snan"                 # initialize REALs to signaling NaN
  " -finit-integer=-2147483648"       # initialize INTEGERs to sentinel value

  " -fsanitize=undefined"             # detect undefined behavior
  # " -fsanitize=address,leak,thread" # optional, heavier sanitizers
)
