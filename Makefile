SHELL = /bin/sh

###################################################################################################
# One-command-build invoking CMake (meant for developers, should not be part of the distribution)
###################################################################################################

.PHONY: all
all: grid mesh

.PHONY: grid
grid:
	@cmake -B build/grid -DDAMASK_SOLVER=grid -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/grid --parallel --target install

.PHONY: mesh
mesh:
	@cmake -B build/mesh -DDAMASK_SOLVER=mesh -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/mesh --parallel --target install

.PHONY: test
test:
	@cmake -B build/test -DDAMASK_SOLVER=test -DCMAKE_INSTALL_PREFIX=${PWD} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE=${BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}
	@cmake --build build/test --parallel --target install

.PHONY: clean
clean:
	@rm -rf build
