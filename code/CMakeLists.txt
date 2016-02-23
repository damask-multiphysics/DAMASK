# group sources for base modules
set (SRC    "CPFEM"
            "CPFEM2"
            "core_quit"
            "commercialFEM_fileList"
            "compilation_info"
            "constitutive"
            "crystallite"
            "damask_hdf5.f90"
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
      add_library (${p} MODULE "${p}.f90")
endforeach (p)

# compile each sub moudel
add_subdirectory(damage)
add_subdirectory(homogenization)
add_subdirectory(hydrogenflux)
add_subdirectory(kinematics)
add_subdirectory(plastic)
add_subdirectory(porosity)
add_subdirectory(sources)
add_subdirectory(spectral)
add_subdirectory(thermal)
add_subdirectory(vacancyflux)

# compile spectral solver 
add_executable(DAMASKSpectral.exe DAMASK_spectral.f90)
