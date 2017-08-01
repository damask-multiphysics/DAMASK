# sets up an environment for DAMASK on tcsh
# usage:  source DAMASK_env.csh

set CALLED=($_)
set DIRNAME=`dirname $CALLED[2]`

# transition compatibility (renamed $DAMASK_ROOT/DAMASK_env.csh to $DAMASK_ROOT/env/DAMASK.csh)
set FILENAME=`basename $CALLED[2]`
if ($FILENAME == "DAMASK.csh") then
  set DIRNAME=$DIRNAME"/../"
endif

set DAMASK_ROOT=`python -c "import os,sys; print(os.path.realpath(os.path.expanduser(sys.argv[1])))" $DIRNAME`


source $DAMASK_ROOT/CONFIG

# if DAMASK_BIN is present
if ( $?DAMASK_BIN) then
  set path = ($DAMASK_BIN $path)
endif

set SOLVER=`which DAMASK_spectral`                                                          
set PROCESSING=`which postResults`                                                          
if ( "x$DAMASK_NUM_THREADS" == "x" ) then                                                                  
  set DAMASK_NUM_THREADS=1
endif

# currently, there is no information that unlimited causes problems
# still,  http://software.intel.com/en-us/forums/topic/501500 suggest to fix it
# http://superuser.com/questions/220059/what-parameters-has-ulimit
limit datasize  unlimited  # maximum  heap size (kB)
limit stacksize unlimited  # maximum stack size (kB)
endif
if ( `limit | grep memoryuse` != "" ) then
  limit memoryuse  unlimited  # maximum physical memory size
endif
if ( `limit | grep vmemoryuse` != "" ) then
  limit vmemoryuse unlimited  # maximum virtual memory size
endif

# disable output in case of scp
if ( $?prompt ) then
  echo ''
  echo Düsseldorf Advanced Materials Simulation Kit --- DAMASK
  echo Max-Planck-Institut für Eisenforschung GmbH, Düsseldorf
  echo https://damask.mpie.de
  echo
  echo Using environment with ...
  echo "DAMASK             $DAMASK_ROOT"
  echo "Spectral Solver    $SOLVER" 
  echo "Post Processing    $PROCESSING"
  echo "Multithreading     DAMASK_NUM_THREADS=$DAMASK_NUM_THREADS"
  if ( $?PETSC_DIR) then
    echo "PETSc location     $PETSC_DIR"
  endif
  if ( $?MSC_ROOT) then
    echo "MSC.Marc/Mentat    $MSC_ROOT"
  endif
  echo
  echo `limit datasize`
  echo `limit stacksize`
endif

setenv DAMASK_NUM_THREADS $DAMASK_NUM_THREADS
if ( ! $?PYTHONPATH ) then
  setenv PYTHONPATH $DAMASK_ROOT/lib
else
  setenv PYTHONPATH $DAMASK_ROOT/lib:$PYTHONPATH
endif
