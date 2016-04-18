SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: spectral FEM


spectral: build
	@if [ ! -d build/spectral ]; then mkdir build/spectral; fi
	(cd build/spectral; cmake ../../ ;)

build:
	mkdir build

.PHONY: FEM
FEM:
	$(MAKE) DAMASK_FEM.exe -C src

.PHONY: marc
marc:
	@./installation/mods_MarcMentat/apply_DAMASK_modifications.sh ${MAKEFLAGS}

.PHONY: processing
processing:
	@if hash cython 2>/dev/null; then \
		cd ./lib/damask; \
	    CC=gcc python setup_corientation.py build_ext --inplace; \
		rm -rv build; \
		rm *.c; \
	fi
	@./installation/compile_CoreModule.py ${MAKEFLAGS}

.PHONY: tidy
tidy:
	@$(MAKE) tidy -C src >/dev/null

.PHONY: clean


clean:
	rm -rvf build
