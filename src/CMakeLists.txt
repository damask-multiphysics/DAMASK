# special flags for some files
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # long lines for interaction matrix
    SET_SOURCE_FILES_PROPERTIES("lattice.f90" PROPERTIES COMPILE_FLAGS "-ffree-line-length-240")
endif()

file(GLOB damask-sources *.f90 *.c)

# probably we should have a subfolder for MSC.Marc
list(FILTER damask-sources EXCLUDE REGEX ".*CPFEM.f90")
list(FILTER damask-sources EXCLUDE REGEX ".*DAMASK_Marc.*.f90")
list(FILTER damask-sources EXCLUDE REGEX ".*commercialFEM_fileList.*.f90")


if (PROJECT_NAME STREQUAL "damask-grid")

    file(GLOB grid-sources grid/*.f90)

    if (NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
        add_executable(DAMASK_grid ${damask-sources} ${grid-sources})
        install (TARGETS DAMASK_grid RUNTIME DESTINATION bin)
    else ()
        add_library(DAMASK_grid OBJECT ${damask-sources} ${grid-sources})
        exec_program (mktemp OUTPUT_VARIABLE nothing)
        exec_program (mktemp ARGS -d OUTPUT_VARIABLE black_hole)
        install (PROGRAMS ${nothing} DESTINATION ${black_hole})
    endif ()

elseif (PROJECT_NAME STREQUAL "damask-mesh")
      
    file(GLOB mesh-sources mesh/*.f90)

    add_executable(DAMASK_mesh ${damask-sources} ${mesh-sources})
    install (TARGETS DAMASK_mesh RUNTIME DESTINATION bin)
      
endif ()
