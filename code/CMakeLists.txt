# The dependency detection for Fortran is not working in CMake
# Have to define them explicitly
add_library(DAMASK_prec "prec.f90")

###################################################################################################

###################################################################################################
if (SPECTRAL)
      add_library(DAMASK_interface "spectral_interface.f90")
      target_link_libraries(DAMASK_interface DAMASK_prec)
endif(SPECTRAL)

add_library(DAMASK_IO "IO.f90")
target_link_libraries(DAMASK_IO DAMASK_interface)

add_library(DAMASK_LIBS "libs.f90")
target_link_libraries(DAMASK_LIBS DAMASK_IO)

add_library(DAMASK_NUMERICS "numerics.f90")
target_link_libraries(DAMASK_NUMERICS DAMASK_LIBS)

add_library(DAMASK_DEBUG "debug.f90")
target_link_libraries(DAMASK_DEBUG DAMASK_NUMERICS)

add_library(DAMASK_BASICS ALIAS DAMASK_DEBUG)

# group sources for base modules
# the FEM modules would require special attention
# will take care of it later.
set (SRC    "CPFEM"
            "CPFEM2"
            "core_quit"
            "commercialFEM_fileList"
            "compilation_info"
            "constitutive"
            "crystallite"
            "damask_hdf5"
            "lattice"
            "material"
            "math"
            "mesh"
            "quit__genmod"
    )

# compiler base modules
foreach (p ${SRC})
      add_library (${p} "${p}.f90")
      target_link_libraries(${p} DAMASK_BASICS)
endforeach (p)

# set libraries/modules for linking
foreach (p ${SRC})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)

# compile each sub moudel
add_subdirectory(damage)
add_subdirectory(homogenization)
add_subdirectory(hydrogenflux)
add_subdirectory(kinematics)
add_subdirectory(plastic)
add_subdirectory(porosity)
add_subdirectory(source)
add_subdirectory(spectral)
add_subdirectory(thermal)
add_subdirectory(vacancyflux)

# compile spectral solver
add_executable(DAMASKSpectral.exe DAMASK_spectral.f90)
target_link_libraries (DAMASKSpectral.exe ${DAMASK_LIB})