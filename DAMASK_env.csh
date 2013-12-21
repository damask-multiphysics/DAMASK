# sets up an environment for DAMASK on tcsh
# usage:  source DAMASK_env.csh
set MAGIG=($_)
set FILENAME=`readlink -f $called[2]`
set LOCATION = `dirname $FILENAME`
setenv DAMASK_ROOT ${LOCATION}
setenv DAMASK_NUM_THREADS 2
# disable output in case of scp
if($?prompt) then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
  echo Max-Planck-Institut für Eisenforschung, Düsseldorf
  echo http://damask.mpie.de
  echo
  echo Preparing environment ...
  echo "DAMASK_ROOT=$DAMASK_ROOT"
  echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
endif
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited
setenv DAMASK_BIN ${DAMASK_ROOT}:bin
setenv PATH ${PATH}:${DAMASK_BIN}
setenv PYTHONPATH ${PYTHONPATH}:${DAMASK_ROOT}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}
