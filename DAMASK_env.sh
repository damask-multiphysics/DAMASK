# sets up an environment for DAMASK on bash
# usage:  source damask_env.sh
export DAMASK_ROOT=$HOME/DAMASK                                                                     
export DAMASK_BIN=$DAMASK_ROOT/bin
export PYTHONPATH=${PYTHONPATH}:${DAMASK_ROOT}/lib
export DAMASK_NUM_THREADS=2
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited

# dissable output in case of scp
[ -z "$PS1" ] && return
echo
echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
echo Max-Planck-Institut für Eisenforschung, Düsseldorf
echo http://damask.mpie.de
echo
echo Preparing environment ...
echo "DAMASK_ROOT=$DAMASK_ROOT"
echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"


