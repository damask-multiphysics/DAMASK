SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
DAMASK_ROOT = $(shell python -c "import os,sys; print(os.path.normpath(os.path.realpath(os.path.expanduser('$(pwd)'))))")
.PHONY: all
all: grid mesh python processing

.PHONY: grid
grid: build/grid
	@(cd build/grid;make -j4 all install;)
.PHONY: spectral
spectral: grid

.PHONY: mesh
mesh: build/mesh
	@(cd build/mesh; make -j4 all install;)
.PHONY: FEM
FEM: mesh

.PHONY: build/grid
build/grid:
	@mkdir -p build/grid
	@(cd build/grid; cmake -Wno-dev -DDAMASK_SOLVER=GRID -DCMAKE_INSTALL_PREFIX=${DAMASK_ROOT} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP} ../../;)

.PHONY: build/mesh
build/mesh:
	@mkdir -p build/mesh
	@(cd build/mesh; cmake -Wno-dev -DDAMASK_SOLVER=FEM -DCMAKE_INSTALL_PREFIX=${DAMASK_ROOT}  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP} ../../;)

.PHONY: clean
clean:
	@rm -rf build

.PHONY: processing
processing:
	@./installation/symlink_Processing.py ${MAKEFLAGS}

.PHONY: python
python:
	@pip install --user ./python

