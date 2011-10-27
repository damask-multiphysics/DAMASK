#!/bin/bash
# Sets environment variable DAMASK_ROOT to allow python to locate various scripts.
# run with --bashrc to set DAMASK_ROOT permanently.
# should maybe got to /etc/profile.d/ for permanent "installation"
#http://stackoverflow.com/questions/630372/determine-the-path-of-the-executing-bash-script
MY_PATH="`dirname \"$0\"`"              # relative
MY_PATH="`( cd \"$MY_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$MY_PATH" ] ; then
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  echo "set_damask_root.sh failed"
  exit 1  # fail
fi
#echo "$MY_PATH"
DAMASK_ROOT=$MY_PATH
export DAMASK_ROOT
echo "DAMASK_ROOT: $DAMASK_ROOT"

PYTHONPATH=$DAMASK_ROOT:$PYTHONPATH
export PYTHONPATH
#echo "PYTHONPATH: $PYTHONPATH"

if [ "$1" = "--bashrc" ] ; then
  echo -e "\nDAMASK_ROOT=$DAMASK_ROOT; export DAMASK_ROOT" >> $HOME/.bashrc
  echo -e "\nPYTHONPATH=$DAMASK_ROOT:$PYTHONPATH; export PYTHONPATH" >> $HOME/.bashrc
  else
  echo "Run with option --bashrc to set DAMASK_ROOT and adjust PYTHONPATH permanently through .bashrc"
fi
echo "Now starting a new child bash that inherits the appropriate DAMASK_ROOT and PYTHONPATH variables. Have even more fun ..."
echo "To run the automated tests call ./testing/run_tests.py"
bash
