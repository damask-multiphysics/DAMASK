# The dependency detection in CMake is not functioning for Fortran
# !!! EXPLICIT DEPENDENCY DECLARATION !!!

add_library(DAMASK_PREC "prec.f90")

if (SPECTRAL)
      add_library(DAMASK_INTERFACE "spectral_interface.f90")
      target_link_libraries(DAMASK_INTERFACE DAMASK_PREC)
elseif(SPECTRAL)
      message(FATAL_ERROR "NOT IMPLEMENTED YET")
endif(SPECTRAL)

add_library(DAMASK_IO "IO.f90")
target_link_libraries(DAMASK_IO DAMASK_INTERFACE)

add_library(DAMASK_LIBS "libs.f90")
target_link_libraries(DAMASK_LIBS DAMASK_IO)

add_library(DAMASK_NUMERICS "numerics.f90")
target_link_libraries(DAMASK_NUMERICS DAMASK_LIBS)

add_library(DAMASK_DEBUG "debug.f90")
target_link_libraries(DAMASK_DEBUG DAMASK_NUMERICS)

add_library(DAMASK_FEsolving "FEsolving.f90")
target_link_libraries(DAMASK_FEsolving DAMASK_DEBUG)

add_library(DAMASK_MATH "math.f90")
target_link_libraries(DAMASK_MATH DAMASK_FEsolving)

# SPECTRAL solver and FEM solver use different mesh
# source files
if (SPECTRAL)
      add_library(DAMASK_MESH "mesh.f90")
endif(SPECTRAL)
if (FEM)
      add_library(DAMASK_MESH "${FEM_DIR}/meshFEM.f90")
endif(FEM)
target_link_libraries(DAMASK_MESH DAMASK_MATH)

add_library(DAMASK_MATERIAL "material.f90")
target_link_libraries(DAMASK_MATERIAL DAMASK_MESH)

add_library(DAMASK_LATTICE "lattice.f90")
target_link_libraries(DAMASK_LATTICE DAMASK_MATERIAL)
add_library(DAMASK_DRIVERS ALIAS DAMASK_LATTICE)

# For each modular section
add_library (DAMASK_PLASTIC "plastic_dislotwin.f90"
                            "plastic_disloUCLA.f90"
                            "plastic_isotropic.f90"
                            "plastic_j2.f90"
                            "plastic_phenopowerlaw.f90"
                            "plastic_titanmod.f90"
                            "plastic_nonlocal.f90"
                            "plastic_none.f90"
                            "plastic_phenoplus.f90")
target_link_libraries(DAMASK_PLASTIC DAMASK_DRIVERS)

add_library (DAMASK_KINEMATICS "kinematics_cleavage_opening.f90"
                               "kinematics_slipplane_opening.f90"
                               "kinematics_thermal_expansion.f90"
                               "kinematics_vacancy_strain.f90"
                               "kinematics_hydrogen_strain.f90")
target_link_libraries(DAMASK_KINEMATICS DAMASK_DRIVERS)

add_library (DAMASK_SOURCE   "source_thermal_dissipation.f90"
                             "source_thermal_externalheat.f90"
                             "source_damage_isoBrittle.f90"
                             "source_damage_isoDuctile.f90"
                             "source_damage_anisoBrittle.f90"
                             "source_damage_anisoDuctile.f90"
                             "source_vacancy_phenoplasticity.f90"
                             "source_vacancy_irradiation.f90"
                             "source_vacancy_thermalfluc.f90")
target_link_libraries(DAMASK_SOURCE DAMASK_DRIVERS)

add_library(DAMASK_CONSTITUTIVE "constitutive.f90")
target_link_libraries(DAMASK_CONSTITUTIVE DAMASK_PLASTIC DAMASK_KINEMATICS DAMASK_SOURCE)

# compile spectral solver
add_executable(DAMASKSpectral.exe DAMASK_spectral.f90)
target_link_libraries (DAMASKSpectral.exe DAMASK_CONSTITUTIVE)