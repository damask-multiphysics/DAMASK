SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: spectral FEM marc
spectral: build/spectral
	@(cd build/spectral; make --no-print-directory -ws all install VERBOSE=1;)

build/spectral: build
	@mkdir build/spectral
	@(cd build/spectral; cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DDAMASK_SOLVER=SPECTRAL ../..;)

FEM: build/FEM
	@(cd build/FEM; make --no-print-directory -ws all install;)

build/FEM: build
	@mkdir build
	@(cd build/FEM; cmake -Wno-dev -DCMAKE_BUILD_TYPE=RELEASE -DDAMASK_SOLVER=FEM ../..;)

build: bin
	@mkdir build

bin:
	@mkdir bin


.PHONY: marc
marc:
	@./installation/mods_MarcMentat/apply_DAMASK_modifications.sh ${MAKEFLAGS}


.PHONY: tidy
tidy:
	@$(MAKE) tidy -C code >/dev/null

.PHONY: clean
clean:
	@$(MAKE) cleanDAMASK -C code >/dev/null

.PHONY: install
install:
	@./installation/symlink_Code.py ${MAKEFLAGS}
	@./installation/symlink_Processing.py ${MAKEFLAGS}

