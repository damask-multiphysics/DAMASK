# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# LLVM Compiler
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 22)

set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEBUG      "${_OPTIMIZATION_OFF}")
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

set(STANDARD_CHECK "-std=f2018 -pedantic")

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

string(APPEND COMPILE_FLAGS
  " -fPIE"                          # position independent code
)
