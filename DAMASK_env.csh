# sets up an environment for DAMASK on csh
# usage:  source DAMASK_env.csh
setenv DAMASK_ROOT $HOME/DAMASK                                                                     
setenv DAMASK_BIN $DAMASK_ROOT/bin
setenv PYTHONPATH ${PYTHONPATH}:${DAMASK_ROOT}/lib
setenv DAMASK_NUM_THREADS 2
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited

echo
echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
echo Max-Planck-Institut für Eisenforschung, Düsseldorf
echo http://damask.mpie.de
echo
echo Preparing environment ...
echo "DAMASK_ROOT=$DAMASK_ROOT"
echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"

