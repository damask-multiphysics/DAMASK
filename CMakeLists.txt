cmake_minimum_required(VERSION 3.18.0)
set(CMAKE_CXX_STANDARD 23)
include(FindPkgConfig REQUIRED)

if(DEFINED ENV{PETSC_DIR})
  message("PETSC_DIR:\n$ENV{PETSC_DIR}\n")
else()
  message(FATAL_ERROR "PETSc location (PETSC_DIR) is not defined")
endif()

# Dummy project to determine compiler names and version
project(Prerequisites LANGUAGES)
set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_check_modules(PETSC_MIN REQUIRED PETSc>=3.15.0 QUIET) #CMake does not support version range
pkg_check_modules(PETSC REQUIRED PETSc<3.24.0)

pkg_get_variable(CMAKE_Fortran_COMPILER PETSc fcompiler)
pkg_get_variable(CMAKE_C_COMPILER PETSc ccompiler)
pkg_get_variable(CMAKE_CXX_COMPILER PETSc cxxcompiler)

# Get version
file(STRINGS ${PROJECT_SOURCE_DIR}/VERSION DAMASK_VERSION_STRING)
string(REPLACE "-" ";" DAMASK_VERSION_LIST "${DAMASK_VERSION_STRING}")
list(LENGTH DAMASK_VERSION_LIST CHUNKS)
if(CHUNKS GREATER 2)
  list(GET DAMASK_VERSION_LIST 0 DAMASK_VERSION_MAIN)
  list(GET DAMASK_VERSION_LIST 1 DAMASK_VERSION_TWEAK)
  list(GET DAMASK_VERSION_LIST 2 PROJECT_VERSION_HASH)
  set(DAMASK_VERSION_CMAKE "${DAMASK_VERSION_MAIN}.${DAMASK_VERSION_TWEAK}")
else()
  set(DAMASK_VERSION_CMAKE "${DAMASK_VERSION_STRING}")
endif()
# Solver determines name of project
string(TOUPPER "${DAMASK_SOLVER}" DAMASK_SOLVER_UPPER)
string(TOLOWER "${DAMASK_SOLVER}" DAMASK_SOLVER_LOWER)
if("${DAMASK_SOLVER_UPPER}" MATCHES "^(GRID|MESH|TEST)$")
  project("damask-${DAMASK_SOLVER_LOWER}"
          DESCRIPTION "${DAMASK_SOLVER_LOWER} solver for DAMASK"
          HOMEPAGE_URL https://damask-multiphysics.org
          LANGUAGES Fortran C CXX VERSION "${DAMASK_VERSION_CMAKE}")
else()
  message(FATAL_ERROR "Invalid solver: DAMASK_SOLVER=${DAMASK_SOLVER}")
endif()
add_compile_definitions("-D${DAMASK_SOLVER_UPPER}")

set(CMAKE_Fortran_PREPROCESS "ON")

# EXPERIMENTAL: This might help to detect external packages in future if PETSc is not aware of them
set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/externalpackages:$ENV{PKG_CONFIG_PATH}")
pkg_check_modules(HDF5 hdf5)
pkg_check_modules(fYAML libfyaml)
pkg_check_modules(LAPACK lapack)
find_package(Boost COMPONENTS program_options)
if(DAMASK_SOLVER_UPPER STREQUAL "GRID")
  pkg_check_modules(zlib zlib)
  pkg_check_modules(FFTW3 fftw3)
endif()

message("\nBuilding ${CMAKE_PROJECT_NAME} ${DAMASK_VERSION_STRING}\n")

add_compile_definitions(-DPETSC)
add_compile_definitions(-DDAMASK_VERSION="${DAMASK_VERSION_STRING}")
add_compile_definitions(-DDAMASK_VERSION_MAJOR="${PROJECT_VERSION_MAJOR}")
add_compile_definitions(-DDAMASK_VERSION_MINOR="${PROJECT_VERSION_MINOR}")
add_compile_definitions(-DDAMASK_VERSION_PATCH="${PROJECT_VERSION_PATCH}")
if(DEFINED PROJECT_VERSION_HASH)
  add_compile_definitions(-DDAMASK_VERSION_HASH="${PROJECT_VERSION_HASH}")
endif()

add_compile_definitions(-DCMAKE_SYSTEM_NAME="${CMAKE_SYSTEM_NAME}")
add_compile_definitions(-DCMAKE_SYSTEM_PROCESSOR="${CMAKE_SYSTEM_PROCESSOR}")

if(PETSC_VERSION VERSION_EQUAL 3.17.0)
  add_compile_definitions("-DCHKERRQ=PetscCall")
endif()
if(PETSC_VERSION VERSION_LESS 3.23.0)
  add_compile_definitions(
    $<$<COMPILE_LANGUAGE:Fortran>:VecGetArray=VecGetArrayF90>
    $<$<COMPILE_LANGUAGE:Fortran>:VecRestoreArray=VecRestoreArrayF90>

    $<$<COMPILE_LANGUAGE:Fortran>:VecGetArrayRead=VecGetArrayReadF90>
    $<$<COMPILE_LANGUAGE:Fortran>:VecRestoreArrayRead=VecRestoreArrayReadF90>

    $<$<COMPILE_LANGUAGE:Fortran>:DMDAVecGetArray=DMDAVecGetArrayF90>
    $<$<COMPILE_LANGUAGE:Fortran>:DMDAVecRestoreArray=DMDAVecRestoreArrayF90>

    $<$<COMPILE_LANGUAGE:Fortran>:DMDAVecGetArrayRead=DMDAVecGetArrayReadF90>
    $<$<COMPILE_LANGUAGE:Fortran>:DMDAVecRestoreArrayRead=DMDAVecRestoreArrayReadF90>

    $<$<COMPILE_LANGUAGE:Fortran>:ISGetIndices=ISGetIndicesF90>
    $<$<COMPILE_LANGUAGE:Fortran>:ISRestoreIndices=ISRestoreIndicesF90>

    $<$<COMPILE_LANGUAGE:Fortran>:SNES_CONVERGED_USER=7>
    $<$<COMPILE_LANGUAGE:Fortran>:SNES_CONVERGED_ITERATING=0>
    $<$<COMPILE_LANGUAGE:Fortran>:SNES_DIVERGED_USER=-12>

    $<$<COMPILE_LANGUAGE:Fortran>:PETSCDTSIMPLEXQUAD_DEFAULT=-1>
    $<$<COMPILE_LANGUAGE:Fortran>:INSERT_ALL_VALUES=5>
  )
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)

# Predefined sets for OPTIMIZATION/OPENMP based on BUILD_TYPE
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
  set(DEBUG_FLAGS "${DEBUG_FLAGS} -DDEBUG")
  set(PARALLEL "OFF")
  set(OPTI "DEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  set(PARALLEL "ON")
  set(OPTI "DEFENSIVE")
elseif(CMAKE_BUILD_TYPE STREQUAL "DEBUGRELEASE")
  set(DEBUG_FLAGS "${DEBUG_FLAGS} -DDEBUG")
  set(PARALLEL "ON")
  set(OPTI "DEFENSIVE")
elseif(CMAKE_BUILD_TYPE STREQUAL "PERFORMANCE")
  set(PARALLEL "ON")
  set(OPTI "AGGRESSIVE")
else()
  message(FATAL_ERROR "Invalid build type: CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
endif()

# $OPTIMIZATION takes precedence over $BUILD_TYPE defaults
if(OPTIMIZATION STREQUAL "" OR NOT DEFINED OPTIMIZATION)
  set(OPTIMIZATION "${OPTI}")
else()
  set(OPTIMIZATION "${OPTIMIZATION}")
endif()

# $OPENMP takes precedence over $BUILD_TYPE defaults
if(OPENMP STREQUAL "" OR NOT DEFINED OPENMP)
  set(OPENMP "${PARALLEL}")
else()
  set(OPENMP "${OPENMP}")
endif()

# syntax check only(mainly for pre-receive hook)
if(CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
  set(BUILDCMD_POST "${BUILDCMD_POST} -fsyntax-only")
endif()


if((CMAKE_Fortran_COMPILER_ID STREQUAL "GNU") OR
   (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM") OR
   (CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang"))
  include("${PROJECT_SOURCE_DIR}/cmake/Compiler-${CMAKE_Fortran_COMPILER_ID}.cmake")
else()
  message(FATAL_ERROR "Compiler '${CMAKE_Fortran_COMPILER_ID}' not supported")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS Fortran_COMPILER_VERSION_MIN)
  message(FATAL_ERROR "Version '${CMAKE_Fortran_COMPILER_VERSION}' of '${CMAKE_Fortran_COMPILER_ID}' is not supported")
endif()


file(STRINGS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/petsc/conf/petscvariables" PETSC_EXTERNAL_LIB REGEX "PETSC_EXTERNAL_LIB_BASIC = .*$?")
string(REPLACE "PETSC_EXTERNAL_LIB_BASIC = " "" PETSC_EXTERNAL_LIB "${PETSC_EXTERNAL_LIB}")
message("PETSC_EXTERNAL_LIB:\n${PETSC_EXTERNAL_LIB}\n")

file(STRINGS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/petsc/conf/petscvariables" PETSC_INCLUDES REGEX "PETSC_FC_INCLUDES = .*$?")
string(REPLACE "PETSC_FC_INCLUDES = " "" PETSC_INCLUDES "${PETSC_INCLUDES}")
message("PETSC_INCLUDES:\n${PETSC_INCLUDES}\n")

set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} "${BUILDCMD_PRE} ${OPENMP_FLAGS} ${STANDARD_CHECK} ${OPTIMIZATION_FLAGS} ${COMPILE_FLAGS}")
set(CMAKE_Fortran_LINK_EXECUTABLE           "${BUILDCMD_PRE} ${CMAKE_Fortran_COMPILER} ${OPENMP_FLAGS} ${OPTIMIZATION_FLAGS} ${LINKER_FLAGS}")

set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${PETSC_INCLUDES}")   # zlib installed via PETSc

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${DEBUG_FLAGS}")
  set(CMAKE_Fortran_LINK_EXECUTABLE           "${CMAKE_Fortran_LINK_EXECUTABLE} ${DEBUG_FLAGS}")
endif()

set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${PETSC_INCLUDES} ${BUILDCMD_POST}")
set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}     "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}     ${PETSC_INCLUDES} ${BUILDCMD_POST}")

set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} <LINK_FLAGS> <FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES> -L${PETSC_LIBRARY_DIRS} -lpetsc ${PETSC_EXTERNAL_LIB}")
if(DAMASK_SOLVER_UPPER STREQUAL "GRID")
  set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} -lz")
endif()

if(fYAML_FOUND STREQUAL "1")
  set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} -L${fYAML_LIBRARY_DIRS}")
  foreach(fYAML_LIBRARY ${fYAML_LIBRARIES})
    set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} -L${fYAML_LIBRARY_DIRS} -l${fYAML_LIBRARY}")
  endforeach()
  add_compile_definitions(-DFYAML)
  pkg_get_variable(fYAML_INCLUDE_DIR libfyaml includedir) # fYAML_INCLUDE_DIRS and fYAML_CFLAGS are not working
  set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} -I${fYAML_INCLUDE_DIR}")
endif()

# Link with Fortran compiler requires to add C++ standard library
set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} -lstdc++ ${BUILDCMD_POST}")

message("Fortran Compiler Flags:\n${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\n")
message("C Compiler Flags:\n${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\n")
message("C++ Compiler Flags:\n${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}\n")
message("Fortran Linker Command:\n${CMAKE_Fortran_LINK_EXECUTABLE}\n")

# location of code
add_subdirectory(src)
