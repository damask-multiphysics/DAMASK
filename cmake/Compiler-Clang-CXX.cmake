# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# LLVM C++ (clang++, AppleClang)
###################################################################################################
set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEBUG      "${_OPTIMIZATION_OFF}")
set(_OPTIMIZATION_DEFENSIVE  "-O2 -mtune=native")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -march=native -funroll-loops -ftree-vectorize -flto")
set(CXX_OPTIMIZATION_FLAGS "${_OPTIMIZATION_${OPTIMIZATION}}")

# AppleClang is not built with OpenMP support, so -fopenmp is not available there.
if(OPENMP AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(CXX_OPENMP_FLAGS "-fopenmp")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/_Compiler-C-family.cmake")
