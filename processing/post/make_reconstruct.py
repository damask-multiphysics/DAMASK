#!/bin/bash

# This script is used to compile the python module used for geometry reconstruction.
# It uses the fortran wrapper f2py that is included in the numpy package to construct the
# module reconstruct.so out of the fortran code reconstruct.f90
# written by M. Diehl, m.diehl@mpie.de

f2py -c -m reconstruct reconstruct.f90
