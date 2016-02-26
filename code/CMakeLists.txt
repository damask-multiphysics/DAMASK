
include_directories(${PROJECT_SRC_DIR}/code)
include_directories(${PROJECT_SRC_DIR}/code/damage)
include_directories(${PROJECT_SRC_DIR}/code/homogenization)
include_directories(${PROJECT_SRC_DIR}/code/hydrogenflux)
include_directories(${PROJECT_SRC_DIR}/code/kinematics)
include_directories(${PROJECT_SRC_DIR}/code/plastic)
include_directories(${PROJECT_SRC_DIR}/code/porosity)
include_directories(${PROJECT_SRC_DIR}/code/source)
include_directories(${PROJECT_SRC_DIR}/code/spectral)
include_directories(${PROJECT_SRC_DIR}/code/thermal)
include_directories(${PROJECT_SRC_DIR}/code/vacancyflux)

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
            "debug"
            "FEsolving"
            "IO"
            "lattice"
            "libs"
            "material"
            "math"
            "mesh"
            "numerics"
            "prec"
            "quit__genmod"
    )

# compiler base modules
foreach (p ${SRC})
      add_library (${p} "${p}.f90")
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