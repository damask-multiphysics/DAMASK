FROM ubuntu:latest

RUN <<SysReq
    apt-get update
    apt-get install -y wget git software-properties-common make cmake pkg-config
    add-apt-repository ppa:ubuntu-toolchain-r/test -y
    apt-get update
SysReq

ARG GCC_V=12
RUN apt-get install -y gcc-${GCC_V} gfortran-${GCC_V} g++-${GCC_V}
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-${GCC_V} 100 \
          --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-${GCC_V} \
          --slave /usr/bin/g++      g++      /usr/bin/g++-${GCC_V} \
          --slave /usr/bin/gcov     gcov     /usr/bin/gcov-${GCC_V}

ARG PETSC_VERSION=3.18.4
RUN <<GetPETSc
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz
    tar -xf petsc-${PETSC_VERSION}.tar.gz -C .
    rm petsc-${PETSC_VERSION}.tar.gz
GetPETSc
ENV PETSC_DIR=${PWD}/petsc-${PETSC_VERSION}
ENV PETSC_ARCH=gcc${GCC_V}
WORKDIR /petsc-${PETSC_VERSION}
RUN ./configure --with-fc=gfortran --with-cc=gcc --with-cxx=g++ --download-openmpi --download-fftw --download-hdf5 --download-hdf5-fortran-bindings=1 --download-zlib --with-mpi-f90module-visibility=1 --download-fblaslapack=1
RUN make all

RUN git clone https://github.com/eisenforschung/DAMASK.git /DAMASK
WORKDIR /DAMASK

# DAMASK_grid
RUN <<DAMASK_grid
    cmake -B build/grid -DDAMASK_SOLVER=grid -DCMAKE_INSTALL_PREFIX=${PWD}
    cmake --build build/grid --parallel
    cmake --install build/grid
    cp ./bin/DAMASK_grid /usr/bin/
DAMASK_grid
# DAMASK_mesh
RUN <<DAMASK_mesh
    cmake -B build/mesh -DDAMASK_SOLVER=mesh -DCMAKE_INSTALL_PREFIX=${PWD}
    cmake --build build/mesh --parallel
    cmake --install build/mesh
    cp ./bin/DAMASK_mesh /usr/bin/
DAMASK_mesh


WORKDIR /wd

## Tests
RUN DAMASK_grid -l tensionX.yaml -g 20grains16x16x16.vti -m material.yaml -w /DAMASK/examples/grid
RUN DAMASK_mesh -h

ENTRYPOINT ["DAMASK_grid"]
CMD ["bash"]