# sets up an environment for DAMASK on bash
# usage:  source DAMASK_env.sh

if [ "$OSTYPE" == "linux-gnu" ] || [ "$OSTYPE" == 'linux' ]; then
  DAMASK_ROOT=$(readlink -f "`dirname $BASH_SOURCE`")
else
  STAT=$(stat "`dirname $BASH_SOURCE`")
  DAMASK_ROOT=${STAT##* }
  unset STAT
fi

[[ -f $HOME/.damask/damask.conf ]] && source $HOME/.damask/damask.conf || source /etc/damask.conf


# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
  echo Max-Planck-Institut für Eisenforschung, Düsseldorf
  echo http://damask.mpie.de
  echo
  echo Using environment with ...
  echo "DAMASK           $DAMASK_ROOT"
  echo "Multithreading   DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
  echo "Compiler         F90=$F90"
  ([[ "x$IMKL_ROOT"   != "x" ]] && echo "IMKL             $IMKL_ROOT") || \
  ([[ "x$ACML_ROOT"   != "x" ]] && echo "ACML             $ACML_ROOT") || \
  ([[ "x$LAPACK_ROOT" != "x" ]] && echo "LAPACK           $LAPACK_ROOT")
  echo "MSC.Marc/Mentat  $MSC_ROOT"
  echo "FFTW             $FFTW_ROOT"
  echo "HDF5             $HDF5_ROOT (for future use)"
  echo
fi

ulimit -s unlimited 2>/dev/null
ulimit -c 0         2>/dev/null
ulimit -v unlimited 2>/dev/null
ulimit -m unlimited 2>/dev/null

export PYTHONPATH=$DAMASK_ROOT/lib:$PYTHONPATH

for var in DAMASK IMKL ACML LAPACK MSC FFTW HDF5; do
  unset "${var}_ROOT"
done

