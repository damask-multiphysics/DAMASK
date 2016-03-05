#!/bin/bash

cat  README
echo
echo "Building FEM solver with ${FC}"
DAMASKVERSION=$(cat VERSION)

# prepare building directory
# structure:
#   BUILD_DIR
#   |-BUILD_SPECTRAL
#   |-BUILD_FEM
#   |-BUILD_MARC
if [ ! -d build ]; then
    mkdir build
fi
cd build
if [ -d FEM ] ; then
    rm -rf FEM
fi
mkdir FEM
cd FEM

##
# CMake call
# PETSC_DIR          |  PETSC directory
# DAMASK_V           |  DAMASK current revision
# CMAKE_BUILD_TYPE   |  Default set to release (no debugging output)
# OPENMP             |  [ON/OFF]
# OPTIMIZATION       |  [OFF,DEFENSIVE,AGGRESSIVE,ULTRA]
# DAMASK_DRIVER      |  [SPECTRAL, FEM]
# DAMASK_INSTALL     |  Directory to install binary output
cmake -D PETSC_DIR=${PETSC_DIR}             \
      -D DAMASK_V=${DAMASKVERSION}          \
      -D CMAKE_BUILD_TYPE=RELEASE           \
      -D OPENMP=ON                          \
      -D OPTIMIZATION=DEFENSIVE             \
      -D DAMASK_DRIVER=FEM                  \
      -D DAMASK_INSTALL=${HOME}/bin         \
      ../..

echo
echo "Please move to the build directory using"
echo "    cd build/FEM"
echo "Using the following command to build DAMASK spectral solver"
echo "    make clean all install"
