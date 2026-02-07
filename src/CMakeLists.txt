# SPDX-License-Identifier: AGPL-3.0-or-later
function(fs_case_sensitive base_path out_var)
  set(_test_dir "${base_path}/neveruse")
  file(MAKE_DIRECTORY "${_test_dir}")
  file(WRITE "${_test_dir}/never.use" "test\n")
  file(WRITE "${_test_dir}/NEVER.USE" "test\n")

  # Count distinct directory entries
  file(GLOB _case_test_files "${_test_dir}/*")
  list(LENGTH _case_test_files _case_test_count)

  # If both names coexist as separate entries -> case-sensitive
  set(_is_case_sensitive FALSE)
  if(_case_test_count GREATER 1)
    set(_is_case_sensitive TRUE)
  endif()

  file(REMOVE_RECURSE "${_test_dir}")
  set(${out_var} ${_is_case_sensitive} PARENT_SCOPE)
endfunction()

foreach(VARIANT GRID MESH TEST)
  if(${VARIANT})
    string(TOLOWER "${VARIANT}" VARIANT_LOWER)
    set(executable_name "DAMASK_${VARIANT_LOWER}")
    if(VARIANT_LOWER STREQUAL "test")
      set(variant_path "${PROJECT_SOURCE_DIR}/tests/unit_fortran")
    else()
      set(variant_path ${VARIANT_LOWER})
    endif()

    file(GLOB damask_fortran CONFIGURE_DEPENDS *.f90 ${variant_path}/*.f90)
    file(GLOB damask_c CONFIGURE_DEPENDS *.c *.cpp ${variant_path}/*.c ${variant_path}/*.cpp)

    if(NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
      add_library(damask_fortran_${VARIANT_LOWER} OBJECT ${damask_fortran})
      target_compile_definitions(damask_fortran_${VARIANT_LOWER} PUBLIC ${VARIANT})
      set_target_properties(damask_fortran_${VARIANT_LOWER}
        PROPERTIES Fortran_MODULE_DIRECTORY
          ${CMAKE_CURRENT_BINARY_DIR}/${VARIANT}_${CMAKE_BUILD_TYPE})
      add_executable(${executable_name} $<TARGET_OBJECTS:damask_fortran_${VARIANT_LOWER}>)
      set_target_properties(${executable_name} PROPERTIES LINKER_LANGUAGE Fortran)

      add_library(damask_c_${VARIANT_LOWER} STATIC ${damask_c})
      target_compile_definitions(damask_c_${VARIANT_LOWER} PUBLIC ${VARIANT})
      if(Boost_FOUND)
        target_link_libraries(damask_c_${VARIANT_LOWER} PUBLIC Boost::program_options)
        target_compile_definitions(damask_c_${VARIANT_LOWER} PUBLIC BOOST)
        add_compile_definitions(BOOST)
      endif()
      target_link_libraries(${executable_name} PRIVATE damask_c_${VARIANT_LOWER})
      install(TARGETS ${executable_name} RUNTIME DESTINATION bin)
      fs_case_sensitive(${CMAKE_INSTALL_PREFIX}/bin fs_is_case_sensitive)
      if(fs_is_case_sensitive)
        INSTALL(CODE "execute_process( \
          COMMAND ${CMAKE_COMMAND} -E create_symlink \
          ${executable_name} \
          ${CMAKE_INSTALL_PREFIX}/bin/damask_${VARIANT_LOWER} \
          )")
      endif()
    else()
      add_library(${executable_name} OBJECT ${damask_fortran} ${damask_c})
      exec_program(mktemp OUTPUT_VARIABLE nothing)
      exec_program(mktemp ARGS -d OUTPUT_VARIABLE black_hole)
      install(PROGRAMS ${nothing} DESTINATION ${black_hole})
    endif()
  endif()
  string(REPLACE ";" "\n" sources_${VARIANT_LOWER} "${damask_fortran};${damask_c};${damask_c_${VARIANT_LOWER}}")
  file(WRITE ${CMAKE_BINARY_DIR}/sources_${VARIANT_LOWER}.txt ${sources_${VARIANT_LOWER}})
endforeach()
