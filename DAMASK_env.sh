# sets up an environment for DAMASK on bash
# usage:  source DAMASK_env.sh

if [ "$OSTYPE" == "linux-gnu" ] || [ "$OSTYPE" == 'linux' ]
  then DAMASK_ROOT=$(readlink -f "`dirname $BASH_SOURCE`")
else
  STAT=$(stat "`dirname $BASH_SOURCE`")
  DAMASK_ROOT=${STAT##* }
fi

if [ -f $HOME/.damask/damask.conf ]; then
   source $HOME/.damask/damask.conf
else
   source /etc/damask.conf
fi

if [ "x$DAMASK_NUM_THREADS" != "x" ] 
  then export DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS
fi

if [ "x$F90" != "x" ] 
  then export F90=$F90
fi

# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
  echo Max-Planck-Institut für Eisenforschung, Düsseldorf
  echo http://damask.mpie.de
  echo
  echo Preparing environment ...
  echo "DAMASK installation in $DAMASK_ROOT"
  echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
  echo "F90=$F90"
fi
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited
export PYTHONPATH=$PYTHONPATH:$DAMASK_ROOT/lib
