#!/usr/bin/env bash

#==================================================================================================
# Execute this script (type './DAMASK_prerequisites.sh') 
# and send system_report.txt to damask@mpie.de for support
#==================================================================================================

OUTFILE="system_report.txt"
echo ===========================================
echo +  Generating $OUTFILE                    
echo +  Send to damask@mpie.de for support
echo +  view with \'cat $OUTFILE\'
echo ===========================================

function getDetails {
if which $1 &> /dev/null; then
  echo ----------------------------------------------------------------------------------------------
  echo $1:
  echo ----------------------------------------------------------------------------------------------
  echo + location:
  which $1
  echo + $1 $2:
  $1 $2
  echo -e '\n'
else
  echo $ does not exist
fi
}

# redirect STDOUT and STDERR to logfile
# https://stackoverflow.com/questions/11229385/redirect-all-output-in-a-bash-script-when-using-set-x^
exec > $OUTFILE 2>&1
 
# directory, file is not a symlink by definition
# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DAMASK_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
echo System report for \'$(hostname)\' created on $(date '+%Y-%m-%d %H:%M:%S') by \'$(whoami)\'
echo XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
echo
echo ==============================================================================================
echo DAMASK settings
echo ==============================================================================================
echo ----------------------------------------------------------------------------------------------
echo DAMASK_ROOT:
echo ----------------------------------------------------------------------------------------------
echo $DAMASK_ROOT
echo
echo ----------------------------------------------------------------------------------------------
echo Version:
echo ----------------------------------------------------------------------------------------------
cat  VERSION
echo
echo ----------------------------------------------------------------------------------------------
echo Settings in CONFIG:
echo ----------------------------------------------------------------------------------------------
cat  CONFIG
echo
echo ==============================================================================================
echo System
echo ==============================================================================================
uname -a
echo
echo PATH: $PATH
echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
echo PYTHONPATH: $PYTHONPATH
echo SHELL: $SHELL
echo PETSC_ARCH: $PETSC_ARCH
echo PETSC_DIR: $PETSC_DIR
echo 
echo ==============================================================================================
echo Python
echo ==============================================================================================

DEFAULT_PYTHON=python2.7
for executable in python python2 python3 python2.7; do
  getDetails $executable '--version'
done
echo ----------------------------------------------------------------------------------------------
echo Details on $DEFAULT_PYTHON:
echo ----------------------------------------------------------------------------------------------
echo $(ls -la $(which $DEFAULT_PYTHON))
for module in numpy scipy;do
  echo -e '\n----------------------------------------------------------------------------------------------'
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
  getDetails $executable '--version'
done
echo
echo ==============================================================================================
echo Intel Compiler Suite
echo ==============================================================================================
for executable in icc icpc ifort ;do
  getDetails $executable '--version'
done
echo
echo ==============================================================================================
echo MPI Wrappers
echo ==============================================================================================
for executable in mpicc mpiCC mpic++ mpicpc mpicxx mpifort mpif90 mpif77; do
  getDetails $executable '-show'
done
echo
echo ==============================================================================================
echo MPI Launchers
echo ==============================================================================================
for executable in mpirun mpiexec; do
  getDetails $executable '--version'
done
echo
echo ==============================================================================================
echo Abaqus
echo ==============================================================================================
cd installation/mods_Abaqus                                                                         # to have the right environment file
for executable in abaqus abq2016 abq2017; do
  getDetails $executable 'information=all'
done
cd ../..

