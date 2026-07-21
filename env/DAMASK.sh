# SPDX-License-Identifier: AGPL-3.0-or-later
# set up an environment for DAMASK
# usage:  source DAMASK.sh

canonical_path() {
  python3 -c "import os,sys; print(os.path.normpath(os.path.realpath(os.path.expanduser(sys.argv[1]))))" "$1"
}

print_memory_limit() {
  python3 -c \
"import math,sys
if sys.argv[1] == 'unlimited':
    print('unlimited')
else:
    size=1024*int(sys.argv[1]); \
    print('{:.4g} {}'.format(size / (1 << ((int(math.log(size,2) / 10) if size else 0) * 10)), \
      ['bytes','KiB','MiB','GiB','TiB','EiB','ZiB'][int(math.log(size,2) / 10) if size else 0]))" "$1"
}

blink() {
  printf '\033[2;5m%s\033[0m\n' "$1"
}

if [ -n "${BASH_VERSINFO:-}" ]; then # bash
  ENV_ROOT=$(canonical_path "${BASH_SOURCE[0]}/../")
elif [ -n "${ZSH_VERSION:-}" ]; then # zsh
  ENV_ROOT=$(canonical_path "${0:a:h}")
else # fallback
  ENV_ROOT=$(canonical_path "${0}")
fi
DAMASK_ROOT=$(canonical_path "$ENV_ROOT/../")

# shorthand command to change to DAMASK_ROOT directory
eval "DAMASK_root() { cd \"$DAMASK_ROOT\"; }"

# add BRANCH if DAMASK_ROOT is a git repository
BRANCH="$(git branch C "${DAMASK_ROOT}" --show-current 2>/dev/null || true)"

PATH=${DAMASK_ROOT}/bin${PATH:+:$PATH}

SOLVER_GRID=$(command -v damask_grid || true 2>/dev/null)
[ -z "$SOLVER_GRID" ] && SOLVER_GRID=$(blink 'Not found!')
SOLVER_MESH=$(command -v damask_mesh || true 2>/dev/null)
[ -z "$SOLVER_MESH" ] && SOLVER_MESH=$(blink 'Not found!')


# currently, there is no information that unlimited stack size causes problems
# still, http://software.intel.com/en-us/forums/topic/501500 suggest to fix it
# more info https://jblevins.org/log/segfault
#           https://stackoverflow.com/questions/79923/what-and-where-are-the-stack-and-heap
#           http://superuser.com/questions/220059/what-parameters-has-ulimit
ulimit -s unlimited 2>/dev/null # maximum stack size (kB)

[ -z "$OMP_NUM_THREADS" ] && export OMP_NUM_THREADS=4
[ -z "$OPENBLAS_NUM_THREADS" ] && export OPENBLAS_NUM_THREADS=1 # avoid nested threads
[ -z "$I_MPI_JOB_ABORT_SIGNAL" ] && export I_MPI_JOB_ABORT_SIGNAL=15 # SIGTERM
[ -z "$I_MPI_JOB_SIGNAL_PROPAGATION" ] && export I_MPI_JOB_SIGNAL_PROPAGATION=yes

# disable output in case of scp
if [ "$-" != "${-#*i}" ]; then
  echo
  echo "Düsseldorf Advanced Materials Simulation Kit --- DAMASK"
  echo "Max-Planck-Institut für Nachhaltige Materialien GmbH, Düsseldorf"
  echo "https://damask-multiphysics.org"
  echo
  echo "Using environment with ..."
  echo "DAMASK             $DAMASK_ROOT $BRANCH"
  echo "Grid Solver        $SOLVER_GRID"
  echo "Mesh Solver        $SOLVER_MESH"
  if [ -n "${PETSC_DIR:-}" ]; then
    printf "PETSc location     "
    [ -d "$PETSC_DIR" ] && echo "$PETSC_DIR" || blink "$PETSC_DIR"
    [ "$(canonical_path "$PETSC_DIR")" = "$PETSC_DIR" ] \
    || echo "               ~~> $(canonical_path "$PETSC_DIR")"
  fi
  [ -n "${PETSC_ARCH:-}" ] && echo "PETSc architecture $PETSC_ARCH"
  echo "Multithreading     OMP_NUM_THREADS=$OMP_NUM_THREADS"
  echo "                   OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
  echo "IntelMPI           I_MPI_JOB_ABORT_SIGNAL=$I_MPI_JOB_ABORT_SIGNAL"
  echo "                   I_MPI_JOB_SIGNAL_PROPAGATION=$I_MPI_JOB_SIGNAL_PROPAGATION"
  printf "heap  size         "; print_memory_limit "$(ulimit -d)"
  printf "stack size         "; print_memory_limit "$(ulimit -s)"
  echo
fi

export DAMASK_ROOT
export PYTHONPATH=${DAMASK_ROOT}/python${PYTHONPATH:+:$PYTHONPATH}

if [ -n "${BASH_VERSINFO:-}" ]; then
  # shellcheck source=/dev/null
  source "$ENV_ROOT/damask_grid"
  # shellcheck source=/dev/null
  source "$ENV_ROOT/damask_mesh"
elif [ -n "${ZSH_VERSION:-}" ]; then
  export FPATH=${DAMASK_ROOT}/env${FPATH:+:$FPATH}
fi

for var in SOLVER_GRID SOLVER_MESH BRANCH ENV_ROOT; do
  unset "${var}"
done
