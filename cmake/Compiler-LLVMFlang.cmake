###################################################################################################
# LLVM Compiler
###################################################################################################
if (OPENMP)
  set (OPENMP_FLAGS "-fopenmp")
endif ()

set (STANDARD_CHECK "-std=f2018 -pedantic" )

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
set (COMPILE_FLAGS "${COMPILE_FLAGS} -cpp") # preprocessor, needed for CMake < 3.18
