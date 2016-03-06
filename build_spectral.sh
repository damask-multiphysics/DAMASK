#!/bin/bash

cat  README
echo

if [ "$OSTYPE" == "linux-gnu" ] || [ "$OSTYPE" == 'linux' ]; then
  DAMASK_ROOT=$(readlink -f "`dirname $BASH_SOURCE`")
else
  [[ "${BASH_SOURCE::1}" == "/" ]] && BASE="" || BASE="`pwd`/"
  STAT=$(stat "`dirname $BASE$BASH_SOURCE`")
  DAMASK_ROOT=${STAT##* }
fi

DAMASKVERSION=$(cat VERSION)
BUILDROOT=$DAMASK_ROOT/build
BUILDDIR=spectral

# prepare building directory
# structure:
#   BUILD_DIR
#   |-BUILD_SPECTRAL
#   |-BUILD_FEM
#   |-BUILD_MARC
if [ ! -d $BUILDROOT ]; then
    mkdir $BUILDROOT
fi
cd $BUILDROOT
if [ -d $BUILDDIR ] ; then
    rm -rf $BUILDDIR
fi
mkdir $BUILDDIR
cd $BUILDDIR

##
# CMake call
# PETSC_DIR                |  PETSC directory
# DAMASK_V                 |  DAMASK current revision
# CMAKE_BUILD_TYPE         |  Default set to release (no debugging output)
# OPENMP                   |  [ON/OFF]
# OPTIMIZATION             |  [OFF,DEFENSIVE,AGGRESSIVE,ULTRA]
# DAMASK_DRIVER            |  [SPECTRAL, FEM]
# DAMASK_INSTALL           |  Directory to install binary output
cmake -D PETSC_DIR=${PETSC_DIR}           \
      -D DAMASK_V=${DAMASKVERSION}        \
      -D CMAKE_BUILD_TYPE=RELEASE         \
      -D OPENMP=ON                        \
      -D OPTIMIZATION=DEFENSIVE           \
      -D DAMASK_DRIVER=SPECTRAL           \
      -D DAMASK_INSTALL=${HOME}/bin       \
      ../..

echo
echo "Please move to the build directory using"
echo "    cd build/spectral"
echo "Using the following command to build DAMASK spectral solver"
echo "    make clean all install"
