###################################################################################################
# PGI Compiler
###################################################################################################
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")

  if (OPTIMIZATION STREQUAL "OFF")
    set (OPTIMIZATION_FLAGS "-O0"       )
  elseif (OPTIMIZATION STREQUAL "DEFENSIVE")
    set (OPTIMIZATION_FLAGS "-O2")
  elseif (OPTIMIZATION STREQUAL "AGGRESSIVE")
    set (OPTIMIZATION_FLAGS "-O3")
  endif ()


#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
  set (COMPILE_FLAGS "${COMPILE_FLAGS} -Mpreprocess")
  # preprocessor

  set (STANDARD_CHECK "-Mallocatable=03")

#------------------------------------------------------------------------------------------------
# Runtime debugging
  set (DEBUG_FLAGS "${DEBUG_FLAGS} -g")
  # Includes debugging information in the object module; sets the optimization level to zero unless a -⁠O option is present on the command line
