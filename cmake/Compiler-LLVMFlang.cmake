# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# LLVM Compiler
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 22)

if(OPENMP)
  set(OPENMP_FLAGS "-fopenmp")
endif()

if(OPTIMIZATION STREQUAL "DEBUG")
  set(OPTIMIZATION_FLAGS "-O0")
elseif(OPTIMIZATION STREQUAL "OFF")
  set(OPTIMIZATION_FLAGS "-O0")
elseif(OPTIMIZATION STREQUAL "DEFENSIVE")
  set(OPTIMIZATION_FLAGS "-O2 -mtune=native")
elseif(OPTIMIZATION STREQUAL "AGGRESSIVE")
  set(OPTIMIZATION_FLAGS "-O3 -march=native -funroll-loops -ftree-vectorize -flto")
endif()

set(STANDARD_CHECK "-std=f2018 -pedantic")

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

# position independent code:
set(COMPILE_FLAGS "${COMPILE_FLAGS} -fPIE")
