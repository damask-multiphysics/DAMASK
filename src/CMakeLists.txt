file(GLOB damask_c CONFIGURE_DEPENDS
  *.c *.cpp ${DAMASK_SOLVER_LOWER}/*.c ${DAMASK_SOLVER_LOWER}/*.cpp)
file(GLOB damask_fortran CONFIGURE_DEPENDS *.f90 ${DAMASK_SOLVER_LOWER}/*.f90)

set(executable-name "DAMASK_${DAMASK_SOLVER_LOWER}")

if(NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
  add_library(damask_c STATIC ${damask_c})
  if(Boost_FOUND)
    target_link_libraries(damask_c PUBLIC Boost::program_options)
    target_compile_definitions(damask_c PRIVATE BOOST)
  endif()

  add_library(damask_fortran OBJECT ${damask_fortran})
  add_executable(${executable-name} $<TARGET_OBJECTS:damask_fortran>)
  set_target_properties(${executable-name} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${executable-name} PRIVATE damask_c)
  install(TARGETS ${executable-name} RUNTIME DESTINATION bin)
else()
  add_library(${executable-name} OBJECT ${damask_fortran} ${damask_c})
  exec_program(mktemp OUTPUT_VARIABLE nothing)
  exec_program(mktemp ARGS -d OUTPUT_VARIABLE black_hole)
  install(PROGRAMS ${nothing} DESTINATION ${black_hole})
endif()

string(REPLACE ";" "\n" sources "${damask_c};${damask_fortran}")
message(${CMAKE_BINARY_DIR})
file(WRITE ${CMAKE_BINARY_DIR}/sources.txt ${sources})
