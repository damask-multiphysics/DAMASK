file(GLOB damask-sources CONFIGURE_DEPENDS *.f90 *.c *.cpp)

set(executable-name "DAMASK_${DAMASK_SOLVER_LOWER}")
file(GLOB solver-sources CONFIGURE_DEPENDS ${DAMASK_SOLVER_LOWER}/*.f90)

if(NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
  add_executable(${executable-name} ${damask-sources} ${solver-sources})
  install(TARGETS ${executable-name} RUNTIME DESTINATION bin)
else()
  add_library(${executable-name} OBJECT ${damask-sources} ${solver-sources})
  exec_program(mktemp OUTPUT_VARIABLE nothing)
  exec_program(mktemp ARGS -d OUTPUT_VARIABLE black_hole)
  install(PROGRAMS ${nothing} DESTINATION ${black_hole})
endif()

set_target_properties("DAMASK_${DAMASK_SOLVER_LOWER}" PROPERTIES LINKER_LANGUAGE Fortran)

string(REPLACE ";" "\n" sources "${damask-sources};${solver-sources}")
message(${CMAKE_BINARY_DIR})
file(WRITE ${CMAKE_BINARY_DIR}/sources.txt ${sources})
