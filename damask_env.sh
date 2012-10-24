#!/bin/bash
# sets up an environment for DAMASK
# usage:  source damask_env.sh

# It would be nice to set the variables
# depending on the location of the script
# /!\ CAREFUL WITH HARDLINKS
# but the ideas below don't work with sourcing this snippet.
# scriptloc = `dirname $0`
scriptloc=`dirname $(readlink -f $0)`
#echo I am here: $scriptloc
#export DAMASK_ROOT=$scriptloc
#export DAMASK_BIN=$scriptloc/bin


# Until the above works, we assume DAMASK
# lives at the top of the user's home dir
# strip the username from "MPIE\"
me=`echo $USER |cut -d'\\' -f 2` 
#echo user: $me
#export DAMASK_ROOT=`echo /nethome/$me/DAMASK`
#export DAMASK_BIN=`echo /nethome/$me/DAMASK/bin`
export DAMASK_ROOT=~/DAMASK
export DAMASK_BIN=~/DAMASK/bin

echo DAMASK_ROOT = $DAMASK_ROOT
echo DAMASK_BIN = $DAMASK_BIN
# set more variables depending on the above ones
export PATH=${PATH}:${DAMASK_BIN}
export PYTHONPATH=${PYTHONPATH}:${DAMASK_ROOT}/lib
export LD_LIBRARY_PATH=${DAMASK_ROOT}/lib/fftw/lib:lib/acml4.4.0/ifort64_mp/lib:lib/acml4.4.0/ifort64/lib:${LD_LIBRARY_PATH}
