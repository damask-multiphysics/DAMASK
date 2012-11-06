# sets up an environment for DAMASK
# usage:  source damask_env.sh

currentDir=${PWD}
if [ "${BASH_SOURCE[0]:0:1}" = "/" ]
then cd -P `dirname "${BASH_SOURCE[0]}"`
else cd -P `dirname "${PWD}/${BASH_SOURCE[0]}"`
fi
theRoot=${PWD}
cd $currentDir

export DAMASK_ROOT=$theRoot
export DAMASK_BIN=$theRoot/bin
[[ "${PATH}" == *"${DAMASK_BIN}"* ]]            || export PATH=${PATH}:${DAMASK_BIN}
[[ "${PYTHONPATH}" == *"${DAMASK_ROOT}/lib"* ]] || export PYTHONPATH=${PYTHONPATH}:${DAMASK_ROOT}/lib
#export LD_LIBRARY_PATH=${DAMASK_ROOT}/lib/fftw/lib:lib/acml4.4.0/ifort64_mp/lib:lib/acml4.4.0/ifort64/lib:${LD_LIBRARY_PATH}
