####################################################################################################
# One-command-build invoking CMake (meant for developers, should not be part of the distribution)
####################################################################################################
SHELL = /bin/sh
.RECIPEPREFIX = >
ROOTDIR=$(shell pwd)

.PHONY: grid mesh test build clean iwyu clang-tidy clang-format runtest

SELECTED_TARGETS := $(shell printf '%s' '$(MAKECMDGOALS)' | tr '[:lower:]' '[:upper:]')

grid mesh test: build

build:
> @cmake -S . -B $(ROOTDIR)/build \
>   -DGRID=OFF -DMESH=OFF -DTEST=OFF \
>   $(foreach VARIANT,$(if $(SELECTED_TARGETS),$(SELECTED_TARGETS),GRID MESH),-D$(VARIANT)=ON) \
>   -DCMAKE_INSTALL_PREFIX=${PWD} \
>   -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
>   -DBUILDCMD_POST=$(BUILDCMD_POST) \
>   -DBUILDCMD_PRE=$(BUILDCMD_PRE) \
>   -DOPTIMIZATION=$(OPTIMIZATION) \
>   -DOPENMP=$(OPENMP)
> @cmake --build $(ROOTDIR)/build --parallel --target install

clean:
> @rm -rf $(ROOTDIR)/build $(ROOTDIR)/build-clang-tidy

iwyu:
> @cmake -B build-iwyu -DTEST=ON \
>   -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE="include-what-you-use;-Xiwyu;--error;-Xiwyu;--mapping_file=$(/usr/share/include-what-you-use/boost-all.imp);-Xiwyu;--verbose=3"
>   -DCMAKE_INSTALL_PREFIX=$(shell mktemp -d)
> @cmake --build build-iwyu --target install --verbose

clang-tidy:
> @cmake -B build-clang-tidy -DTEST=ON \
>   -DCMAKE_CXX_CLANG_TIDY="clang-tidy;--config-file=$(ROOTDIR)/.clang-tidy"
> @cmake --build build-clang-tidy --target install --verbose

clang-format:
> @clang-format --style=file:.clang-format -i --Werror \
>   $(shell find src tests -name '*.c' -o -name '*.cpp' -o -name '*.h')

runtest:
> @$(ROOTDIR)/bin/damask_test
> @ctest --test-dir build/tests/unit_cpp/ --output-on-failure

.DEFAULT:
> @echo "Error: unknown target '$@'"; exit 2
