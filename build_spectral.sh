#!/bin/bash

cat  README
echo "Building spectral solver with ${FC}"

# prepare building directory
if [ ! -d build_spectral ] ; then
    mkdir build_spectral
fi
cd build_spectral

# start building
cmake -D PETSC_DIR=${PETSC_DIR} \
      -D HDF5=${HDF5} \
      ..

# instruction for compiling
echo "Please go into build_spectral directory and use make"
echo "to build DAMASK_spectal.exe"