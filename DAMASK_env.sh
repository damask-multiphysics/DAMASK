# sets up an environment for DAMASK on bash
# usage:  source DAMASK_env.sh
LOCATION=$(readlink -f "`dirname $BASH_SOURCE`")
export DAMASK_ROOT=${LOCATION}
export DAMASK_NUM_THREADS=2
# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
  echo Max-Planck-Institut für Eisenforschung, Düsseldorf
  echo http://damask.mpie.de
  echo
  echo Preparing environment ...
  echo "DAMASK_ROOT=$DAMASK_ROOT"
  echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
fi
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited
export DAMASK_BIN=${DAMASK_ROOT}/bin
export PATH=${PATH}:${DAMASK_BIN}
export PYTHONPATH=${PYTHONPATH}:${DAMASK_ROOT}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
