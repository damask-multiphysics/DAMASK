#!/usr/bin/env bash

OUTFILE="system_report.txt"
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
  echo ----------------------------------------------------------------------------------------------
  $DEFAULT_PYTHON -c "import $module; \
                      print('Version: {}'.format($module.__version__)); \
                      print('Location: {}'.format($module.__file__))"
done
echo ----------------------------------------------------------------------------------------------
echo vtk
echo ----------------------------------------------------------------------------------------------
$DEFAULT_PYTHON -c "import vtk; \
                    print('Version: {}'.format(vtk.vtkVersion.GetVTKVersion())); \
                    print('Location: {}'.format(vtk.__file__))"
echo ----------------------------------------------------------------------------------------------
echo h5py
echo ----------------------------------------------------------------------------------------------
$DEFAULT_PYTHON -c "import h5py; \
                    print('Version: {}'.format(h5py.version.version)); \
                    print('Location: {}'.format(h5py.__file__))"
echo
echo ==============================================================================================
echo GCC
echo ==============================================================================================
for executable in gcc g++ gfortran ;do
  if [[ "$(which $executable)x" != "x" ]]; then
    echo $(which $executable) version: $($executable --version 2>&1)
  else
     echo $executable does not exist
  fi
done
echo
echo ==============================================================================================
echo Intel Compiler Suite
echo ==============================================================================================
for executable in icc icpc ifort ;do
  if [[ "$(which $executable)x" != "x" ]]; then
    echo $(which $executable) version: $($executable --version 2>&1)
  else
     echo $executable does not exist
  fi
done
echo
echo ==============================================================================================
echo MPI Wrappers
echo ==============================================================================================
for executable in mpicc mpiCC mpicxx mpicxx mpifort mpif90 mpif77; do
  if [[ "$(which $executable)x" != "x" ]]; then
    echo $(which $executable) version: $($executable --show 2>&1)
  else
     echo $executable does not exist
  fi
done
