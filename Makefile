SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: spectral FEM marc processing
spectral: build/spectral
	@(cd build/spectral;make --no-print-directory -ws all install VERBOSE=1;)

FEM: build/FEM
	@(cd build/FEM; make --no-print-directory -ws all install;)


OPTIONS="-Wno-dev -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILDCMD_POST=${BUILDCMD_POST} -DBUILDCMD_PRE={BUILDCMD_PRE} -DOPTIMIZATION=${OPTIMIZATION} -DOPENMP=${OPENMP}"

build/spectral: build
	@mkdir -p build/spectral
	@(cd build/spectral; cmake -DDAMASK_SOLVER=SPECTRAL ${OPTIONS} ../..;)

build/FEM: build
	@mkdir -p build/FEM
	@(cd build/FEM; cmake -DDAMASK_SOLVER=FEM ${OPTIONS} ../..;)


build:
	@mkdir -p build

.PHONY: marc
marc:
	@./installation/symLink_Code.sh
	@./installation/mods_MarcMentat/apply_DAMASK_modifications.sh ${MAKEFLAGS}

.PHONY: clean
clean:
	@rm -rf build

.PHONY: processing
install:
	@./installation/symlink_Processing.py ${MAKEFLAGS}

