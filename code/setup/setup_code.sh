#!/usr/bin/env sh

${DAMASK_ROOT}/code/setup/modify_Files.py
#remove -c or --compile from args an in case call compile spectral solver
compile=0
myArgs=''
for args
do
  if [ "$args" = "-c" -o "$args" = "--compile" ]
    then compile=1
  else
    myArgs="$myArgs $args"
  fi
done
echo $myArgs
if [ "$compile" = 1 ]
 then echo 'Compiling spectral solver' 
 ${DAMASK_ROOT}/code/setup/compile_SpectralSolver.py $myArgs
fi
${DAMASK_ROOT}/code/setup/symlink_Code.py
