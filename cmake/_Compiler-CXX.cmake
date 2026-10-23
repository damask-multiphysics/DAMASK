# SPDX-License-Identifier: AGPL-3.0-or-later
###################################################################################################
# CXX-specific options
###################################################################################################

string(APPEND CXX_DEBUG_FLAGS
  " -Wnon-virtual-dtor"              # warn on non-virtual destructors in polymorphic classes
  " -D_GLIBCXX_ASSERTIONS"           # enable standard library bounds/iterator checks
)
