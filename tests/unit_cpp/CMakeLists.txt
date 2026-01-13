# SPDX-License-Identifier: AGPL-3.0-or-later
find_package(GTest CONFIG REQUIRED)
include(GoogleTest)

add_executable(test_CLI test_CLI.cpp)
set_target_properties(test_CLI PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(test_CLI PRIVATE GTest::gtest_main damask_c_test ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
gtest_discover_tests(test_CLI TEST_PREFIX "cli." DISCOVERY_MODE PRE_TEST)
