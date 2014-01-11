# sets up an environment for DAMASK on bash
# usage:  source DAMASK_env.sh
FFTWROOT=/usr/local
LAPACKROOT=/usr
ACMLROOT=
IMKLROOT=
DAMASKFORTRAN=gfortran

if [ "$OSTYPE" == "linux-gnu" ]
  then LOCATION=$(readlink -f "`dirname $BASH_SOURCE`")
else
  STAT=$(stat "`dirname $BASH_SOURCE`")
  LOCATION=${STAT##* }
fi
export DAMASK_ROOT=${LOCATION}
export DAMASK_NUM_THREADS=2
export FFTWROOT=${FFTWROOT}
LD_NEW=$FFTWROOT/lib
if [ "x$LAPACKROOT" != "x" ] 
  then export LAPACKROOT=$LAPACKROOT
  LD_NEW=$LD_NEW:$LAPACKROOT/lib:$LAPACKROOT/lib64
fi
if [ "x$ACMLROOT" != "x" ] 
  then export ACMLROOT=$ACMLROOT
  LD_NEW=$LD_NEW:$ACMLROOT/ifort64_mp/lib:$ACMLROOT/ifort64/lib:$ACMLROOT/gfortran64_mp/lib:$ACMLROOT/gfortran64/lib
fi
if [ "x$IMKLROOT" != "x" ] 
  then export IMKLROOT=${IMKLROOT}
fi

if [ "x$DAMASKFORTRAN" != "x" ] 
  then export F90=$DAMASKFORTRAN
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
