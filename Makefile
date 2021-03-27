SHELL = /bin/sh

###################################################################################################
# One-command-build invoking CMake (meant for developers, should not be part of the distribution)
###################################################################################################

.PHONY: all
all: grid mesh

.PHONY: grid
grid:
	@cmake -B build/grid -DDAMASK_SOLVER=GRID -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/grid --parallel
	@cmake --install build/grid

.PHONY: mesh
mesh:
	@cmake -B build/mesh -DDAMASK_SOLVER=MESH -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/mesh --parallel
	@cmake --install build/mesh

.PHONY: clean
clean:
	@rm -rf build
