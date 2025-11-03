SHELL = /bin/sh
.RECIPEPREFIX = >

###################################################################################################
# One-command-build invoking CMake (meant for developers, should not be part of the distribution)
###################################################################################################


.PHONY: grid mesh test build clean

SELECTED_TARGETS := $(shell printf '%s' '$(MAKECMDGOALS)' | tr '[:lower:]' '[:upper:]')

grid mesh test: build

build:
> @cmake -S . -B build \
>   -DGRID=OFF -DMESH=OFF -DTEST=OFF \
>   $(foreach VARIANT,$(if $(SELECTED_TARGETS),$(SELECTED_TARGETS),GRID MESH),-D$(VARIANT)=ON) \
>   -DCMAKE_INSTALL_PREFIX=${PWD} \
>   -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
>   -DBUILDCMD_POST=$(BUILDCMD_POST) \
>   -DBUILDCMD_PRE=$(BUILDCMD_PRE) \
>   -DOPTIMIZATION=$(OPTIMIZATION) \
>   -DOPENMP=$(OPENMP)
> @cmake --build build --parallel --target install

clean:
> @rm -rf build

.DEFAULT:
> @echo "Error: unknown target '$@'. Valid targets: [grid, mesh, test]"; exit 2
