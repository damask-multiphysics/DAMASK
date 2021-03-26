SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: grid mesh

.PHONY: grid
grid:
	@cmake -B build/grid -DDAMASK_SOLVER=GRID -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/grid --parallel ${DAMASK_NUM_THRADS}
	@cmake --install build/grid

.PHONY: mesh
mesh:
	@cmake -B build/mesh -DDAMASK_SOLVER=MESH -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/mesh --parallel ${DAMASK_NUM_THRADS}
	@cmake --install build/mesh

.PHONY: clean
clean:
	@rm -rf build
