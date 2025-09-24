###################################################################################################
# LLVM Compiler
###################################################################################################
set(Fortran_COMPILER_VERSION_MIN 20)

if(OPENMP)
  set(OPENMP_FLAGS "-fopenmp")
endif()

set(STANDARD_CHECK "-std=f2018 -pedantic" )
