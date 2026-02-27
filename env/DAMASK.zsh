# SPDX-License-Identifier: AGPL-3.0-or-later
# set up an environment for DAMASK on zsh
# usage:  source DAMASK.zsh

function canonicalPath {
  python3 -c "import os,sys; print(os.path.normpath(os.path.realpath(os.path.expanduser(sys.argv[1]))))" $1
}

function blink {
  echo -e "\033[2;5m$1\033[0m"
}

ENV_ROOT=$(canonicalPath "${0:a:h}")
DAMASK_ROOT=$(canonicalPath "${0:a:h}'/..")

# add BRANCH if DAMASK_ROOT is a git repository
cd $DAMASK_ROOT >/dev/null; BRANCH=$(git branch 2>/dev/null| grep -E '^\* '); cd - >/dev/null

PATH=${DAMASK_ROOT}/bin:$PATH

SOLVER_GRID=$(which damask_grid || true 2>/dev/null)
[[ "x$SOLVER_GRID" == "x" ]] && SOLVER_GRID=$(blink 'Not found!')
SOLVER_MESH=$(which damask_mesh || true 2>/dev/null)
[[ "x$SOLVER_MESH" == "x" ]] && SOLVER_MESH=$(blink 'Not found!')


# currently, there is no information that unlimited stack size causes problems
# still, http://software.intel.com/en-us/forums/topic/501500 suggest to fix it
# more info https://jblevins.org/log/segfault
#           https://stackoverflow.com/questions/79923/what-and-where-are-the-stack-and-heap
#           http://superuser.com/questions/220059/what-parameters-has-ulimit
ulimit -s unlimited 2>/dev/null # maximum stack size (kB)

[[ "x$OMP_NUM_THREADS" == "x" ]] && export OMP_NUM_THREADS=4
[[ "x$OPENBLAS_NUM_THREADS" == "x" ]] && export OPENBLAS_NUM_THREADS=1 # avoid nested threads
[[ "x$I_MPI_JOB_ABORT_SIGNAL" == "x" ]] && export I_MPI_JOB_ABORT_SIGNAL=15 # SIGTERM
[[ "x$I_MPI_JOB_SIGNAL_PROPAGATION" == "x" ]] && export I_MPI_JOB_SIGNAL_PROPAGATION=yes

# disable output in case of scp
if [ ! -z "$PS1" ]; then
  echo
  echo Düsseldorf Advanced Materials Simulation Kit --- DAMASK
  echo Max-Planck-Institut für Nachhaltige Materialien GmbH, Düsseldorf
  echo https://damask-multiphysics.org
  echo
  echo "Using environment with ..."
  echo "DAMASK             $DAMASK_ROOT $BRANCH"
  echo "Grid Solver        $SOLVER_GRID"
  echo "Mesh Solver        $SOLVER_MESH"
  if [ "x$PETSC_DIR" != "x" ]; then
    echo -n "PETSc location     "
    [ -d $PETSC_DIR ] && echo $PETSC_DIR || blink $PETSC_DIR
    [[ $(canonicalPath "$PETSC_DIR") == $PETSC_DIR ]] \
    || echo "               ~~> "$(canonicalPath "$PETSC_DIR")
  fi
  [[ "x$PETSC_ARCH" != "x" ]] && echo "PETSc architecture $PETSC_ARCH"
  echo "Multithreading     OMP_NUM_THREADS=$OMP_NUM_THREADS"
  echo "                   OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
  echo "IntelMPI           I_MPI_JOB_ABORT_SIGNAL=$I_MPI_JOB_ABORT_SIGNAL"
  echo "                   I_MPI_JOB_SIGNAL_PROPAGATION=$I_MPI_JOB_SIGNAL_PROPAGATION"
  echo -n "heap  size         "
   [[ "$(ulimit -d)" == "unlimited" ]] \
   && echo "unlimited" \
   || echo $(python3 -c \
          "import math; \
           size=$(( 1024*$(ulimit -d) )); \
           print('{:.4g} {}'.format(size / (1 << ((int(math.log(size,2) / 10) if size else 0) * 10)), \
           ['bytes','KiB','MiB','GiB','TiB','EiB','ZiB'][int(math.log(size,2) / 10) if size else 0]))")
  echo -n "stack size         "
   [[ "$(ulimit -s)" == "unlimited" ]] \
   && echo "unlimited" \
   || echo $(python3 -c \
          "import math; \
           size=$(( 1024*$(ulimit -s) )); \
           print('{:.4g} {}'.format(size / (1 << ((int(math.log(size,2) / 10) if size else 0) * 10)), \
           ['bytes','KiB','MiB','GiB','TiB','EiB','ZiB'][int(math.log(size,2) / 10) if size else 0]))")
  echo
fi

export DAMASK_ROOT
export PYTHONPATH=$DAMASK_ROOT/python:$PYTHONPATH
export FPATH=$DAMASK_ROOT/env:$FPATH

for var in SOLVER BRANCH ENV_ROOT; do
  unset "${var}"
done
