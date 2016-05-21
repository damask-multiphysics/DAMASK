SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: spectral FEM

spectral: build/spectral
	@(cd build/spectral; make --no-print-directory -ws all install VERBOSE=1;)

build/spectral: build
	@mkdir build/spectral
	@(cd build/spectral; cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DDAMASK_SOLVER=SPECTRAL ../..;)

build: bin
	@mkdir build

bin:
	@mkdir bin

FEM: build/FEM
	@(cd build/FEM; make --no-print-directory -ws all install;)

build/FEM: build
	@mkdir build
	@(cd build/FEM; cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DDAMASK_SOLVER=FEM ../..;)

.PHONY: clean
clean:
	rm -rvf build     # standard build directory
	rm -rvf testing   # for testing build (testing script in PRIVATE)
	rm -rvf bin       # default binary store location
