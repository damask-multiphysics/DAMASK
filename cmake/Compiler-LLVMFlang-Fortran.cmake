# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# LLVM Fortran (flang)
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 23)

set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEBUG      "${_OPTIMIZATION_OFF}")
set(_OPTIMIZATION_DEFENSIVE  "-O2 -mtune=native")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -march=native -funroll-loops -ftree-vectorize -flto")
set(OPTIMIZATION_FLAGS_FORTRAN "${_OPTIMIZATION_${OPTIMIZATION}}")

#if(OPENMP)
#  set(OPENMP_FLAGS_FORTRAN "-fopenmp")
#endif()

set(FORTRAN_STANDARD_CHECK "-std=f2018 -pedantic")

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

string(APPEND FORTRAN_COMPILE_FLAGS
  " -fPIE"                          # position independent code
)
