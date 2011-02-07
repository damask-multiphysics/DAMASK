#!/bin/bash

# This script is used to compile the python module used for geometry reconstruction.
# It uses the fortran wrapper f2py that is included in the numpy package to construct the
# module reconstruct.so out of the fortran code reconstruct.f90
# written by M. Diehl, m.diehl@mpie.de


#f2py -m reconstruct2 -h reconstruct2.pyf --overwrite-signature reconstruct_new.f90
#f2py -m reconstruct -h reconstruct.pyf reconstruct.f90
# f2py -c \
# --f90flags="-heap-arrays 500000000" \  # preventing segmentation fault for large arrays
# reconstruct.pyf \
# reconstruct.f90
f2py -c --f90flags="-heap-arrays 500000000" reconstruct.pyf reconstruct.f90