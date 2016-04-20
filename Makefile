SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: spectral FEM

spectral: build/spectral
	@(cd build/spectral; make -s all install;)

build/spectral: build
	@mkdir build/spectral
	@(cd build/spectral; cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DDAMASK_DRIVER=SPECTRAL ../..;)

build: bin
	@mkdir build

bin:
	@mkdir bin

FEM: build/FEM
	@(cd build/FEM; make -s all install;)

build/FEM: build
	@mkdir build
	@(cd build/FEM; cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DDAMASK_DRIVER=FEM ../..;)

.PHONY: clean
clean:
	rm -rvf build
