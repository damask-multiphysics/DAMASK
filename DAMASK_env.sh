# sets up an environment for DAMASK on bash
# usage:  source DAMASK_env.sh

if [ "$OSTYPE" == "linux-gnu" ] || [ "$OSTYPE" == 'linux' ]
  then DAMASK_ROOT=$(readlink -f "`dirname $BASH_SOURCE`")
else
  STAT=$(stat "`dirname $BASH_SOURCE`")
  DAMASK_ROOT=${STAT##* }
  unset STAT
fi

if [ -f $HOME/.damask/damask.conf ]; then
   source $HOME/.damask/damask.conf
else
   source /etc/damask.conf
fi

# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
  echo Max-Planck-Institut für Eisenforschung, Düsseldorf
  echo http://damask.mpie.de
  echo
  echo Using environment with ...
  echo "DAMASK installation in $DAMASK_ROOT"
  echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
  echo "F90=$F90"
  echo "FFTW_ROOT=$FFTW_ROOT"
  if [ "x$LAPACK_ROOT" != "x" ]; then
      echo "LAPACK_ROOT=$LAPACK_ROOT"
  fi
  if [ "x$ACML_ROOT" != "x" ]; then
      echo "ACML_ROOT=$ACML_ROOT"
  fi
  if [ "x$IMKL_ROOT" != "x" ]; then
      echo "IMKL_ROOT=$IMKL_ROOT"
  fi
  echo "MARC_ROOT=$MARC_ROOT"
  echo "HDF5_ROOT=$HDF5_ROOT (future use)"
fi
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited
export PYTHONPATH=$PYTHONPATH:$DAMASK_ROOT/lib
unset DAMASK_ROOT

