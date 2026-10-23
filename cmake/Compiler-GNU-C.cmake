# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# GNU C (gcc)
###################################################################################################
set(_OPTIMIZATION_DEBUG      "-Og")
set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEFENSIVE  "-O2 -mtune=native")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -march=native -funroll-loops -ftree-vectorize -flto")
set(C_OPTIMIZATION_FLAGS "${_OPTIMIZATION_${OPTIMIZATION}}")

if(OPENMP)
  set(C_OPENMP_FLAGS "-fopenmp")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/_Compiler-C-family.cmake")

string(APPEND C_DEBUG_FLAGS
  " -fsanitize=undefined"             # detect undefined behavior
)
