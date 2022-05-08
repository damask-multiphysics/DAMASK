# special flags for some files
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # long lines for interaction matrix
    SET_SOURCE_FILES_PROPERTIES("lattice.f90" PROPERTIES COMPILE_FLAGS "-ffree-line-length-240")
endif()

file(GLOB damask-sources CONFIGURE_DEPENDS *.f90 *.c)

if(PROJECT_NAME STREQUAL "damask-grid")
    set(executable-name "DAMASK_grid")
    file(GLOB solver-sources CONFIGURE_DEPENDS grid/*.f90)
elseif(PROJECT_NAME STREQUAL "damask-mesh")
    set(executable-name "DAMASK_mesh")
    file(GLOB solver-sources CONFIGURE_DEPENDS mesh/*.f90)
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
    add_executable(${executable-name} ${damask-sources} ${solver-sources})
    install(TARGETS ${executable-name} RUNTIME DESTINATION bin)
else()
    add_library(${executable-name} OBJECT ${damask-sources} ${solver-sources})
    exec_program(mktemp OUTPUT_VARIABLE nothing)
    exec_program(mktemp ARGS -d OUTPUT_VARIABLE black_hole)
    install(PROGRAMS ${nothing} DESTINATION ${black_hole})
endif()

string(REPLACE ";" "\n" sources "${damask-sources};${solver-sources}")
message(${CMAKE_BINARY_DIR})
file(WRITE ${CMAKE_BINARY_DIR}/sources.txt ${sources})
