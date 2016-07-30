# sets up an environment for DAMASK on zsh
# usage:  source DAMASK_env.zsh

DAMASK_ROOT=${0:a:h}

# defining set() allows to source the same file for tcsh and zsh, with and without space around =
set() {
    export $1$2$3
 }
source $DAMASK_ROOT/CONFIG
unset -f set

# add DAMASK_BIN if present but not yet in $PATH
MATCH=`echo ":$PATH:" | grep $DAMASK_BIN:`
if [[ ( "x$DAMASK_BIN" != "x" ) && ( "x$MATCH" = "x" ) ]]; then
  export PATH=$DAMASK_BIN:$PATH
fi

SOLVER=`which DAMASK_spectral 2>/dev/null`
PROCESSING=`which postResults 2>/dev/null`
if [ "x$DAMASK_NUM_THREADS" = "x" ]; then
  DAMASK_NUM_THREADS=1
fi

# according to http://software.intel.com/en-us/forums/topic/501500
# this seems to make sense for the stack size
if [ "`which free 2>/dev/null`" != "free not found" ]; then
  freeMem=`free -k | grep -E '(Mem|Speicher):' | awk '{print $4;}'`

  # http://superuser.com/questions/220059/what-parameters-has-ulimit             
  ulimit -s `expr $freeMem / $DAMASK_NUM_THREADS / 2` 2>/dev/null # maximum stack size (kB)
  ulimit -d `expr $freeMem                       / 2` 2>/dev/null # maximum  heap size (kB)
fi
ulimit -c 0           2>/dev/null # core  file size (512-byte blocks)
ulimit -v unlimited   2>/dev/null # maximum virtual memory size
ulimit -m unlimited   2>/dev/null # maximum physical memory size

# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit --- DAMASK
  echo Max-Planck-Institut für Eisenforschung GmbH, Düsseldorf
  echo https://damask.mpie.de
  echo
  echo Using environment with ...
  echo "DAMASK             $DAMASK_ROOT"
  echo "Spectral Solver    $SOLVER" 
  echo "Post Processing    $PROCESSING"
  echo "Multithreading     DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
  if [ "x$PETSC_DIR"   != "x" ]; then
    echo "PETSc location     $PETSC_DIR"
    [[ $(python -c "import os,sys; print(os.path.realpath(os.path.expanduser(sys.argv[1])))" "$PETSC_DIR") == $PETSC_DIR ]] \
    || echo "               ~~> "$(python -c "import os,sys; print(os.path.realpath(os.path.expanduser(sys.argv[1])))" "$PETSC_DIR")
  fi
  [[ "x$PETSC_ARCH"  == "x" ]] \
  || echo "PETSc architecture $PETSC_ARCH"
  echo "MSC.Marc/Mentat    $MSC_ROOT"
  echo
  echo -n "heap  size         "
   [[ "$(ulimit -d)" == "unlimited" ]] \
   && echo "unlimited" \
   || echo $(python -c \
          "import math; \
           size=$(( 1024*$(ulimit -d) )); \
           print '{:.4g} {}'.format(size / (1 << ((int(math.log(size,2) / 10) if size else 0) * 10)), \
           ['bytes','KiB','MiB','GiB','TiB','EiB','ZiB'][int(math.log(size,2) / 10) if size else 0])")
  echo -n "stack size         "
   [[ "$(ulimit -s)" == "unlimited" ]] \
   && echo "unlimited" \
   || echo $(python -c \
          "import math; \
           size=$(( 1024*$(ulimit -s) )); \
           print '{:.4g} {}'.format(size / (1 << ((int(math.log(size,2) / 10) if size else 0) * 10)), \
           ['bytes','KiB','MiB','GiB','TiB','EiB','ZiB'][int(math.log(size,2) / 10) if size else 0])")
fi

export DAMASK_NUM_THREADS
export PYTHONPATH=$DAMASK_ROOT/lib:$PYTHONPATH

for var in BASE STAT SOLVER PROCESSING FREE DAMASK_BIN MATCH; do
  unset "${var}"
done
for var in DAMASK MSC; do
  unset "${var}_ROOT"
done
for var in ABAQUS MARC; do
  unset "${var}_VERSION"
done
