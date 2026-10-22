# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# IntelLLVM C++ (icpx)
###################################################################################################
set(_OPTIMIZATION_OFF        "-O0")
set(_OPTIMIZATION_DEBUG      "${_OPTIMIZATION_OFF}")
set(_OPTIMIZATION_DEFENSIVE  "-O2")
set(_OPTIMIZATION_AGGRESSIVE "-O3 -fp-model strict -xHost") # -ipo/-flto break link
set(CXX_OPTIMIZATION_FLAGS "${_OPTIMIZATION_${OPTIMIZATION}}")

if(OPENMP)
  set(CXX_OPENMP_FLAGS "-fiopenmp")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/_Compiler-C-family.cmake")
