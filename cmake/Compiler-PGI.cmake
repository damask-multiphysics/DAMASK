###################################################################################################
# PGI Compiler
###################################################################################################

if (OPENMP)
  set (OPENMP_FLAGS "-mp")
else ()
  set (OPENMP_FLAGS "-nomp")
endif ()


if (OPTIMIZATION STREQUAL "OFF")
  set (OPTIMIZATION_FLAGS "-O0"       )
elseif (OPTIMIZATION STREQUAL "DEFENSIVE")
  set (OPTIMIZATION_FLAGS "-O2 -fast")
elseif (OPTIMIZATION STREQUAL "AGGRESSIVE")
  set (OPTIMIZATION_FLAGS "-O4 -fast -Mvect=sse")
endif ()

set (STANDARD_CHECK "-Mallocatable=03 -Mstandard")

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
set (COMPILE_FLAGS "${COMPILE_FLAGS} -Mpreprocess")
# preprocessor

set (COMPILE_FLAGS "${COMPILE_FLAGS} -Minfo=all")
# instructs the compiler to produce information on standard error

set (COMPILE_FLAGS "${COMPILE_FLAGS} -Minform=warn")
# instructs the compiler to display error messages at the specified and higher levels

set (COMPILE_FLAGS "${COMPILE_FLAGS} -Mdclchk")
# instructs the compiler to require that all program variables be declared

#------------------------------------------------------------------------------------------------O
# Runtime debugging
set (DEBUG_FLAGS "${DEBUG_FLAGS} -g")
# Includes debugging information in the object module; sets the optimization level to zero unless a -⁠O option is present on the command line
set (DEBUG_FLAGS "${DEBUG_FLAGS} -C")
# Generates code to check array bounds
set (DEBUG_FLAGS "${DEBUG_FLAGS} -Mchkptr")
# Check for NULL pointers (pgf95, pgfortran only)
set (DEBUG_FLAGS "${DEBUG_FLAGS} -Mchkstk")
# Check the stack for available space upon entry to and before the start of a parallel region. Useful when many private variables are declared
set (DEBUG_FLAGS "${DEBUG_FLAGS} -Mbounds")
# Specifies whether array bounds checking is enabled or disabled

#------------------------------------------------------------------------------------------------
#  precision settings
set (PRECISION_FLAGS "${PRECISION_FLAGS} -r8")
# Determines whether the compiler promotes REAL variables and constants to DOUBLE PRECISION
