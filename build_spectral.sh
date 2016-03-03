#!/bin/bash

cat  README
echo
echo "Building spectral solver with ${FC}"
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
if [ -d build_spectral ] ; then
    rm -rf build_spectral
fi
mkdir build_spectral
cd build_spectral

##
# CMake call
# PETSC_DIR         |  PETSC directory
# HDF5_DIR          |  HDF5 library (same compiler for DAMASK)
# DAMASK_V          |  DAMASK current revision
# CMAKE_BUILD_TYPE  |  Default set to release (no debugging output)
# OPENMP            |  [ON/OFF]
# OPTIMIZATION      |  [OFF,DEFENSIVE,AGGRESSIVE,ULTRA]
# DAMASK_DRIVER     |  [SPECTRAL, FEM]
cmake -D PETSC_DIR=${PETSC_DIR}    \
      -D HDF5_DIR=${HDF5_DIR}      \
      -D DAMASK_V=${DAMASKVERSION} \
      -D CMAKE_BUILD_TYPE=RELEASE  \
      -D OPENMP=ON                 \
      -D OPTIMIZATION=DEFENSIVE    \
      -D DAMASK_DRIVER=SPECTRAL    \
      ../..