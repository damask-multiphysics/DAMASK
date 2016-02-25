#!/bin/bash

cat  README
echo "Building spectral solver with ${FC}"
DAMASKVERSION :=$(shell cat ../VERSION)

# prepare building directory
if [ -d build_spectral ] ; then
    rm -rf build_spectral
fi
mkdir build_spectral
cd build_spectral

##
# CMake call
# PETSC_DIR         |  PETSC directory
# HDF5_DIR          |  HDF5 library (same compiler for DAMASK)
# DAMASK_ROOT       |  DAMASK source location
# DAMASK_V          |  DAMASK current revision
# CMAKE_BUILD_TYPE  |  Default set to release (no debugging output)
# OPENMP            |  "ON" will turn on OPENMP support
# OPTIMIZATION      |  [OFF,DEFENSIVE,AGGRESSIVE,ULTRA]
cmake -D PETSC_DIR=${PETSC_DIR}    \
      -D HDF5_DIR=${HDF5_DIR}      \
      -D DAMASK_ROOT=..            \
      -D DAMASK_V=${DAMASKVERSION} \
      -D CMAKE_BUILD_TYPE=RELEASE  \
      -D OPENMP=ON                 \
      -D OPTIMIZATION=DEFENSIVE    \
      ..

# instruction for compiling
echo "Please go into build_spectral directory and use make"
echo "to build DAMASK_spectal.exe"