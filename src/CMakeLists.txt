# special flags for some files
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # long lines for interaction matrix
    SET_SOURCE_FILES_PROPERTIES("lattice.f90" PROPERTIES COMPILE_FLAGS "-ffree-line-length-240")
endif()

file(GLOB_RECURSE sources *.f90 *.c)

# probably we should have subfolders for abaqus and MSC.Marc
list(FILTER sources EXCLUDE REGEX ".*CPFEM\\.f90")
list(FILTER sources EXCLUDE REGEX ".*DAMASK_marc.*\\.f90")
list(FILTER sources EXCLUDE REGEX ".*mesh_marc.*\\.f90")
list(FILTER sources EXCLUDE REGEX ".*mesh_abaqus.*\\.f90")
list(FILTER sources EXCLUDE REGEX ".*commercialFEM_fileList.*\\.f90")


if (PROJECT_NAME STREQUAL "damask-grid")

      # probably we should have subfolders for FEM and spectral
      list(FILTER sources EXCLUDE REGEX ".*DAMASK_FEM.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*FEM_utilities.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*FEM_zoo.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*mesh_FEM.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*FEM_mech.*\\.f90")

      if(NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
        add_executable(DAMASK_spectral ${sources})
        install (TARGETS DAMASK_spectral RUNTIME)
      else()
        add_library(DAMASK_spectral OBJECT ${sources})
        exec_program (mktemp OUTPUT_VARIABLE nothing)
        exec_program (mktemp ARGS -d OUTPUT_VARIABLE black_hole)
        install (PROGRAMS ${nothing} DESTINATION ${black_hole})
      endif()

elseif (PROJECT_NAME STREQUAL "damask-mesh")
      
      # probably we should have subfolders for FEM and spectral
      list(FILTER sources EXCLUDE REGEX ".*DAMASK_grid.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*grid_mech_FEM.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*grid_mech_spectral_basic.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*grid_mech_spectral_polarisation.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*grid_damage_spectral.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*grid_thermal_spectral.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*spectral_utilities.*\\.f90")
      list(FILTER sources EXCLUDE REGEX ".*mesh_grid.*\\.f90")

      add_executable(DAMASK_FEM ${sources})
      install (TARGETS DAMASK_FEM RUNTIME)
      
endif()
