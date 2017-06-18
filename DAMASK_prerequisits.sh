#!/usr/bin/env bash

OUTFILE="bla.txt"
echo date +"%m-%d-%y" >OUTFILE

# redirect STDOUT and STDERR to logfile
# https://stackoverflow.com/questions/11229385/redirect-all-output-in-a-bash-script-when-using-set-x^
exec > $OUTFILE 2>&1

# directory, file is not a symlink by definition
# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DAMASK_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo ==============================================================================================
echo DAMASK settings
echo ==============================================================================================
echo
echo DAMASK_ROOT: $DAMASK_ROOT
echo
echo Settings in CONFIG:
cat  CONFIG
echo
echo ==============================================================================================
echo Python
echo ==============================================================================================

DEFAULT_PYTHON=python2.7
for executable in python python2 python3 python2.7; do
  if [[ "$(which $executable)x" != "x" ]]; then
     echo $executable version: $($executable --version 2>&1)
  else
     echo $executable does not exist
  fi
done
echo Location of $DEFAULT_PYTHON: $(ls -la $(which $DEFAULT_PYTHON))
echo
for module in numpy scipy;do
  echo ----------------------------------------------------------------------------------------------
  echo $module
  $DEFAULT_PYTHON -c "import $module; print('Version: {}'.format($module.__version__));print('Location: {}'.format($module.__file__))"
done
echo ----------------------------------------------------------------------------------------------
echo vtk
$DEFAULT_PYTHON -c "import vtk; print('Location: {}'.format(vtk.__file__))"
