# sets up an environment for DAMASK on bash
# usage:  source DAMASK_env.sh

if [ "$OSTYPE" == "linux-gnu" ] || [ "$OSTYPE" == 'linux' ]
  then LOCATION=$(readlink -f "`dirname $BASH_SOURCE`")
else
  STAT=$(stat "`dirname $BASH_SOURCE`")
  LOCATION=${STAT##* }
fi
export DAMASK_ROOT=${LOCATION}
source $DAMASK_ROOT/installation/options
if [ "x$DAMASK_NUM_THREADS" != "x" ] 
  then export DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS
fi
export FFTW_ROOT=$FFTW_ROOT
LD_NEW=$FFTW_ROOT/lib
if [ "x$LAPACK_ROOT" != "x" ] 
  then LD_NEW=$LD_NEW:$LAPACK_ROOT/lib:$LAPACK_ROOT/lib64
fi
if [ "x$ACML_ROOT" != "x" ] 
  then LD_NEW=$LD_NEW:$ACML_ROOT/ifort64_mp/lib:$ACML_ROOT/ifort64/lib:$ACML_ROOT/gfortran64_mp/lib:$ACML_ROOT/gfortran64/lib
fi
if [ "x$IMKL_ROOT" != "x" ] 
  then LD_NEW=$LD_NEW:$IMKL_ROOT/lib/intel64
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
  echo "DAMASK_ROOT=$DAMASK_ROOT"
  echo "DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
  echo "F90=$F90"
  echo "prepending to LD_LIBRARY_PATH: $LD_NEW"
fi
ulimit -s unlimited
ulimit -c 0
ulimit -v unlimited
ulimit -m unlimited
export DAMASK_BIN=$DAMASK_ROOT/bin
export PATH=$PATH:$DAMASK_BIN
export PYTHONPATH=$PYTHONPATH:$DAMASK_ROOT/lib
export LD_LIBRARY_PATH=$LD_NEW:$LD_LIBRARY_PATH
