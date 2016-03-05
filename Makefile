SHELL = /bin/sh
########################################################################################
# Makefile for the installation of DAMASK
########################################################################################
.PHONY: all
all: spectral marc processing

.PHONY: spectral
spectral:
	$(MAKE) DAMASK_spectral.exe -C src

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
	@$(MAKE) cleanDAMASK -C src >/dev/null

.PHONY: install
install:
	@./installation/symlink_Code.py ${MAKEFLAGS}
	@./installation/symlink_Processing.py ${MAKEFLAGS}

