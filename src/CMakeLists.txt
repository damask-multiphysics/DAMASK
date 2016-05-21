# SOME FILES REQUIRE SPECIAL FLAGS
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    # fno-range-check:  Disable range checking on results of simplification of constant expressions during compilation
    #      -->   allows the definition of DAMASK_NaN
    #-fall-intrinsics:  all intrinsic procedures (including the GNU-specific extensions) are accepted. -Wintrinsics-std will be ignored
    #                   and no user-defined procedure with the same name as any intrinsic will be called except when it is explicitly declared external
    #      -->   allows the use of 'isnan'
    #-fno-fast-math:
    #      -->   otherwise, when setting -ffast-math, isnan always evaluates to false (I would call it a bug)
    SET_SOURCE_FILES_PROPERTIES( "prec.f90" PROPERTIES
                                    COMPILE_FLAGS "-fno-range-check -fall-intrinsics -fno-fast-math")
    # long lines for interaction matrix
    SET_SOURCE_FILES_PROPERTIES( "lattice.f90" PROPERTIES
                                COMPILE_FLAGS "-ffree-line-length-240")
endif()

# The dependency detection in CMake is not functioning for Fortran,
# hence we declare the dependencies from top to bottom in the following

add_library(C_ROUTINES "C_routines.c")

add_library(SYSTEM_ROUTINES "system_routines.f90")
target_link_libraries(SYSTEM_ROUTINES C_ROUTINES)

add_library(PREC "prec.f90")
target_link_libraries(PREC SYSTEM_ROUTINES)

if ("${PROJECT_NAME}" STREQUAL "DAMASK_spectral")
      add_library(DAMASK_INTERFACE "spectral_interface.f90")
elseif ("${PROJECT_NAME}" STREQUAL "DAMASK_FEM")
      add_library(DAMASK_INTERFACE "DAMASK_FEM_interface.f90")
endif()
target_link_libraries(DAMASK_INTERFACE PREC)

add_library(IO "IO.f90")
target_link_libraries(IO DAMASK_INTERFACE)

add_library(EXTERNALLIBS "libs.f90")
target_link_libraries(EXTERNALLIBS IO)

add_library(NUMERICS "numerics.f90")
target_link_libraries(NUMERICS EXTERNALLIBS)

add_library(DEBUG "debug.f90")
target_link_libraries(DEBUG NUMERICS)

add_library(FEsolving "FEsolving.f90")
target_link_libraries(FEsolving DEBUG)

add_library(DAMASK_MATH "math.f90")
target_link_libraries(DAMASK_MATH FEsolving)

# SPECTRAL solver and FEM solver use different mesh
# source files
if ("${PROJECT_NAME}" STREQUAL "DAMASK_spectral")
      add_library(MESH "mesh.f90")
      target_link_libraries(MESH DAMASK_MATH)
elseif ("${PROJECT_NAME}" STREQUAL "DAMASK_FEM")
      add_library(FEZoo "FEZoo.f90")
      target_link_libraries(FEZoo DAMASK_MATH)
      add_library(MESH  "meshFEM.f90")
      target_link_libraries(MESH  FEZoo)
endif()

add_library(MATERIAL "material.f90")
target_link_libraries(MATERIAL MESH)

add_library(LATTICE "lattice.f90")
target_link_libraries(LATTICE MATERIAL)
add_library(DAMASK_HELPERS ALIAS LATTICE)

# For each modular section
add_library (PLASTIC
  "plastic_dislotwin.f90"
  "plastic_disloUCLA.f90"
  "plastic_isotropic.f90"
  "plastic_phenopowerlaw.f90"
  "plastic_titanmod.f90"
  "plastic_nonlocal.f90"
  "plastic_none.f90"
  "plastic_phenoplus.f90")
target_link_libraries(PLASTIC DAMASK_HELPERS)

add_library (KINEMATICS
  "kinematics_cleavage_opening.f90"
  "kinematics_slipplane_opening.f90"
  "kinematics_thermal_expansion.f90"
  "kinematics_vacancy_strain.f90"
  "kinematics_hydrogen_strain.f90")
target_link_libraries(KINEMATICS DAMASK_HELPERS)

add_library (SOURCE
  "source_thermal_dissipation.f90"
  "source_thermal_externalheat.f90"
  "source_damage_isoBrittle.f90"
  "source_damage_isoDuctile.f90"
  "source_damage_anisoBrittle.f90"
  "source_damage_anisoDuctile.f90"
  "source_vacancy_phenoplasticity.f90"
  "source_vacancy_irradiation.f90"
  "source_vacancy_thermalfluc.f90")
target_link_libraries(SOURCE DAMASK_HELPERS)

add_library(CONSTITUTIVE "constitutive.f90")
target_link_libraries(CONSTITUTIVE PLASTIC)
target_link_libraries(CONSTITUTIVE KINEMATICS)
target_link_libraries(CONSTITUTIVE SOURCE)

add_library(CRYSTALLITE "crystallite.f90")
target_link_libraries(CRYSTALLITE CONSTITUTIVE)

add_library(HOMOGENIZATION
  "homogenization_RGC.f90"
  "homogenization_isostrain.f90"
  "homogenization_none.f90")
target_link_libraries(HOMOGENIZATION CRYSTALLITE)

add_library(HYDROGENFLUX
  "hydrogenflux_isoconc.f90"
  "hydrogenflux_cahnhilliard.f90")
target_link_libraries(HYDROGENFLUX CRYSTALLITE)

add_library(POROSITY
  "porosity_none.f90"
  "porosity_phasefield.f90")
target_link_libraries(POROSITY CRYSTALLITE)

add_library(VACANCYFLUX
  "vacancyflux_isoconc.f90"
  "vacancyflux_isochempot.f90"
  "vacancyflux_cahnhilliard.f90")
target_link_libraries(VACANCYFLUX CRYSTALLITE)

add_library(DAMAGE
  "damage_none.f90"
  "damage_local.f90"
  "damage_nonlocal.f90")
target_link_libraries(DAMAGE CRYSTALLITE)

add_library(THERMAL
  "thermal_isothermal.f90"
  "thermal_adiabatic.f90"
  "thermal_conduction.f90")
target_link_libraries(THERMAL CRYSTALLITE)

add_library(DAMASK_ENGINE "homogenization.f90")
target_link_libraries(DAMASK_ENGINE THERMAL       )
target_link_libraries(DAMASK_ENGINE DAMAGE        )
target_link_libraries(DAMASK_ENGINE VACANCYFLUX   )
target_link_libraries(DAMASK_ENGINE POROSITY      )
target_link_libraries(DAMASK_ENGINE HYDROGENFLUX  )
target_link_libraries(DAMASK_ENGINE HOMOGENIZATION)


if ("${PROJECT_NAME}" STREQUAL "DAMASK_spectral")
      add_library(DAMASK_CPFE "CPFEM2.f90")
      target_link_libraries(DAMASK_CPFE DAMASK_ENGINE)

      add_library(SPECTRAL_UTILITIES spectral_utilities.f90)
      target_link_libraries(SPECTRAL_UTILITIES DAMASK_CPFE)

      add_library(SPECTRAL_SOLVER "spectral_thermal.f90"
                                  "spectral_damage.f90"
                                  "spectral_mech_AL.f90"
                                  "spectral_mech_Polarisation.f90"
                                  "spectral_mech_Basic.f90")
      target_link_libraries(SPECTRAL_SOLVER SPECTRAL_UTILITIES)
      add_executable(DAMASK_spectral "DAMASK_spectral.f90")
      target_link_libraries(DAMASK_spectral SPECTRAL_SOLVER)
elseif ("${PROJECT_NAME}" STREQUAL "DAMASK_FEM")
      add_library(DAMASK_CPFE "CPFEM.f90")
      target_link_libraries(DAMASK_CPFE DAMASK_ENGINE)

      add_library(DAMASK_FEM_UTILITY "FEM_utilities.f90")
      target_link_libraries(DAMASK_FEM_UTILITY DAMASK_CPFE)

      add_library(DAMASK_FEM_BASE "FEM_hydrogenflux.f90"
                                  "FEM_porosity.f90"
                                  "FEM_vacancyflux.f90"
                                  "FEM_damage.f90"
                                  "FEM_thermal.f90"
                                  "FEM_mech.f90")
      target_link_libraries(DAMASK_FEM_BASE DAMASK_FEM_UTILITY)

      add_library(DAMASK_FEM_DRIVER "DAMASK_FEM_driver.f90")
      target_link_libraries(DAMASK_FEM_DRIVER DAMASK_FEM_BASE)

      add_executable(DAMASK_FEM "DAMASK_FEM_driver.f90")
      target_link_libraries(DAMASK_FEM DAMASK_FEM_EXE)
endif()
