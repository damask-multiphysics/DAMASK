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

# if DAMASK_BIN is present and not in $PATH, add it
if [[ "x$DAMASK_BIN" != "x" && ! `echo ":$PATH:" | grep $DAMASK_BIN:` ]]; then
  export PATH=$DAMASK_BIN:$PATH
fi

SOLVER=`which DAMASK_spectral`
if [ "x$SOLVER" == "x" ]; then
  export SOLVER='Not found!'
fi
PROCESSING=`which postResults`
if [ "x$PROCESSING" == "x" ]; then
  export PROCESSING='Not found!'
fi
  echo "DAMASK_spectral  `which DAMASK_spectral`"

# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
  echo Max-Planck-Institut für Eisenforschung, Düsseldorf
  echo http://damask.mpie.de
  echo
  echo Using environment with ...
  echo "DAMASK           $DAMASK_ROOT"
  ([[ "x$SOLVER"       != "x" ]] && echo "Spectral Solver  $SOLVER") 
  ([[ "x$PROCESSING"   != "x" ]] && echo "Post Processing  $PROCESSING")
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
# http://superuser.com/questions/220059/what-parameters-has-ulimit
ulimit -s unlimited   2>/dev/null # maximum stack size (kB)
ulimit -h unlimited   2>/dev/null # maximum heap size (kB)
ulimit -c 0           2>/dev/null # core  file size (512-byte blocks)
ulimit -v unlimited   2>/dev/null # maximum virtual memory size
ulimit -m unlimited   2>/dev/null # maximum physical memory size

export DAMASK_NUM_THREADS
export PYTHONPATH=$DAMASK_ROOT/lib:$PYTHONPATH

for var in DAMASK IMKL ACML LAPACK MSC FFTW HDF5; do
  unset "${var}_ROOT"
done

