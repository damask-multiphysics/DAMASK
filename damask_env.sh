#!/bin/bash

# it would be nice to set the variables
script_location = `dirname $0`
echo I am here: $script_location

export DAMASK_ROOT=/nethome/$USER/DAMASK
export DAMASK_BIN=/nethome/$USER/DAMASK/bin
export PATH=${PATH}:${DAMASK_BIN}
export PYTHONPATH=${PYTHONPATH}:${DAMASK_ROOT}/lib
export LD_LIBRARY_PATH=${DAMASK_ROOT}/lib/fftw/lib:lib/acml4.4.0/ifort64_mp/lib:lib/acml4.4.0/ifort64/lib:${LD_LIBRARY_PATH}
