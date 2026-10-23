# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# Shared options for supported C and C++ compilers
###################################################################################################

#------------------------------------------------------------------------------------------------
# Fine tuning compilation options
#------------------------------------------------------------------------------------------------

string(APPEND ${_language}_COMPILE_FLAGS
  " -fPIE"                          # position independent code
  " -Wall"                          # enable common warnings
  " -Wextra"                        # enable extra warnings
)

#------------------------------------------------------------------------------------------------
# Runtime debugging
#------------------------------------------------------------------------------------------------

string(APPEND ${_language}_DEBUG_FLAGS
  " -g"                               # generate debug symbols
  " -fstack-protector-all"            # guard variables on all stack frames
  " -fno-omit-frame-pointer"          # keep frame pointer for accurate backtraces
  " -fstack-clash-protection"         # guard against stack-clash (deep/allocated stack)

  " -ftrivial-auto-var-init=pattern"  # initialize automatic variables to a sentinel pattern

  " -Wshadow"                         # warn on shadowed variables
  " -Wnull-dereference"               # warn on possible null pointer dereferences
  " -Wdouble-promotion"               # warn on implicit promotion of float to double
  " -Wcast-align"                     # warn on casts that increase alignment requirements
  " -Wpointer-arith"                  # warn on pointer arithmetic on unsized types
  " -Wundef"                          # warn on undefined macros in #if directives
  " -Wformat=2"                       # stricter format-string checks
)

if(_language STREQUAL "CXX")
  include("${CMAKE_CURRENT_LIST_DIR}/_Compiler-CXX.cmake")
endif()
