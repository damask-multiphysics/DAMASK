# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# GNU C (gcc)
###################################################################################################
set(_OPTIMIZATION_DEBUG      "-Og")
set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEFENSIVE  "-O2 -mtune=native")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -march=native -funroll-loops -ftree-vectorize -flto")
set(OPTIMIZATION_FLAGS_C "${_OPTIMIZATION_${OPTIMIZATION}}")

if(OPENMP)
  set(OPENMP_FLAGS_C "-fopenmp")
endif()

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

string(APPEND C_COMPILE_FLAGS
  " -fPIE"                          # position independent code
  " -Wall"                          # enable common warnings
  " -Wextra"                        # enable extra warnings
)

#------------------------------------------------------------------------------------------------
# Runtime debugging
#------------------------------------------------------------------------------------------------

string(APPEND C_DEBUG_FLAGS
  " -g"                               # generate debug symbols
  " -fstack-protector-all"            # guard variables on all stack frames
  " -fno-omit-frame-pointer"          # keep frame pointer for accurate backtraces
  " -fstack-clash-protection"         # guard against stack-clash (deep/allocated stack)

  " -ftrivial-auto-var-init=pattern"  # initialize automatic variables to a sentinel pattern

  " -Wshadow"                         # warn on shadowed variables
  " -Wnull-dereference"               # warn on possible null pointer dereferences
  " -Wdouble-promotion"               # warn on implicit promotion of float to double
  " -Wcast-align"                     # warn on casts that increase alignment requirements
  " -Wpointer-arith"                  # warn on pointer arithmetic on unsized types
  " -Wundef"                          # warn on undefined macros in #if directives
  " -Wformat=2"                       # stricter format-string checks

  " -fsanitize=undefined"             # detect undefined behavior
)
