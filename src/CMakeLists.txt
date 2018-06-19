# special flags for some files
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET_SOURCE_FILES_PROPERTIES( "lattice.f90" PROPERTIES
                                COMPILE_FLAGS "-ffree-line-length-240")
    # long lines for interaction matrix
endif()

# The dependency detection in CMake is not functioning for Fortran,
# hence we declare the dependencies from top to bottom in the following
add_library(C_ROUTINES OBJECT "C_routines.c")
set(OBJECTFILES $<TARGET_OBJECTS:C_ROUTINES>)

add_library(SYSTEM_ROUTINES OBJECT "system_routines.f90")
add_dependencies(SYSTEM_ROUTINES C_ROUTINES)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:SYSTEM_ROUTINES>)

add_library(PREC OBJECT "prec.f90")
list(APPEND OBJECTFILES $<TARGET_OBJECTS:PREC>)

if (PROJECT_NAME STREQUAL "DAMASK_spectral")
      add_library(DAMASK_INTERFACE OBJECT "spectral_interface.f90")
elseif (PROJECT_NAME STREQUAL "DAMASK_FEM")
      add_library(DAMASK_INTERFACE OBJECT "FEM_interface.f90")
else ()
    message (FATAL_ERROR "Build target (PROJECT_NAME) is not defined")
endif()
add_dependencies(DAMASK_INTERFACE PREC SYSTEM_ROUTINES)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DAMASK_INTERFACE>)

add_library(IO OBJECT "IO.f90")
add_dependencies(IO DAMASK_INTERFACE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:IO>)

add_library(NUMERICS OBJECT "numerics.f90")
add_dependencies(NUMERICS IO)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:NUMERICS>)

add_library(DEBUG OBJECT "debug.f90")
add_dependencies(DEBUG NUMERICS)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DEBUG>)

add_library(CONFIG OBJECT "config.f90")
add_dependencies(CONFIG DEBUG)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:CONFIG>)

add_library(FEsolving OBJECT "FEsolving.f90")
add_dependencies(FEsolving DEBUG)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:FEsolving>)

add_library(DAMASK_MATH OBJECT "math.f90")
add_dependencies(DAMASK_MATH FEsolving)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DAMASK_MATH>)

# SPECTRAL solver and FEM solver use different mesh files
if (PROJECT_NAME STREQUAL "DAMASK_spectral")
      add_library(MESH OBJECT "mesh.f90")
      add_dependencies(MESH DAMASK_MATH)
      list(APPEND OBJECTFILES $<TARGET_OBJECTS:MESH>)
elseif (PROJECT_NAME STREQUAL "DAMASK_FEM")
      add_library(FEZoo OBJECT "FEZoo.f90")
      add_dependencies(FEZoo DAMASK_MATH)
      list(APPEND OBJECTFILES $<TARGET_OBJECTS:FEZoo>)
      add_library(MESH  OBJECT "meshFEM.f90")
      add_dependencies(MESH FEZoo)
      list(APPEND OBJECTFILES $<TARGET_OBJECTS:MESH>)
endif()

add_library(MATERIAL OBJECT "material.f90")
add_dependencies(MATERIAL MESH CONFIG)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:MATERIAL>)

add_library(DAMASK_HELPERS OBJECT "lattice.f90")
add_dependencies(DAMASK_HELPERS MATERIAL)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DAMASK_HELPERS>)

# For each modular section
add_library (PLASTIC OBJECT
  "plastic_dislotwin.f90"
  "plastic_disloUCLA.f90"
  "plastic_isotropic.f90"
  "plastic_phenopowerlaw.f90"
  "plastic_kinematichardening.f90"
  "plastic_nonlocal.f90"
  "plastic_none.f90")
add_dependencies(PLASTIC DAMASK_HELPERS)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:PLASTIC>)

add_library (KINEMATICS OBJECT
  "kinematics_cleavage_opening.f90"
  "kinematics_slipplane_opening.f90"
  "kinematics_thermal_expansion.f90"
  "kinematics_vacancy_strain.f90"
  "kinematics_hydrogen_strain.f90")
add_dependencies(KINEMATICS DAMASK_HELPERS)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:KINEMATICS>)

add_library (SOURCE OBJECT
  "source_thermal_dissipation.f90"
  "source_thermal_externalheat.f90"
  "source_damage_isoBrittle.f90"
  "source_damage_isoDuctile.f90"
  "source_damage_anisoBrittle.f90"
  "source_damage_anisoDuctile.f90"
  "source_vacancy_phenoplasticity.f90"
  "source_vacancy_irradiation.f90"
  "source_vacancy_thermalfluc.f90")
add_dependencies(SOURCE DAMASK_HELPERS)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:SOURCE>)

add_library(CONSTITUTIVE OBJECT "constitutive.f90")
add_dependencies(CONSTITUTIVE PLASTIC KINEMATICS SOURCE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:CONSTITUTIVE>)

add_library(CRYSTALLITE OBJECT "crystallite.f90")
add_dependencies(CRYSTALLITE CONSTITUTIVE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:CRYSTALLITE>)

add_library(HOMOGENIZATION OBJECT
  "homogenization_RGC.f90"
  "homogenization_isostrain.f90"
  "homogenization_none.f90")
add_dependencies(HOMOGENIZATION CRYSTALLITE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:HOMOGENIZATION>)

add_library(HYDROGENFLUX OBJECT
  "hydrogenflux_isoconc.f90"
  "hydrogenflux_cahnhilliard.f90")
add_dependencies(HYDROGENFLUX CRYSTALLITE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:HYDROGENFLUX>)

add_library(POROSITY OBJECT
  "porosity_none.f90"
  "porosity_phasefield.f90")
add_dependencies(POROSITY CRYSTALLITE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:POROSITY>)

add_library(VACANCYFLUX OBJECT
  "vacancyflux_isoconc.f90"
  "vacancyflux_isochempot.f90"
  "vacancyflux_cahnhilliard.f90")
add_dependencies(VACANCYFLUX CRYSTALLITE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:VACANCYFLUX>)

add_library(DAMAGE OBJECT
  "damage_none.f90"
  "damage_local.f90"
  "damage_nonlocal.f90")
add_dependencies(DAMAGE CRYSTALLITE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DAMAGE>)

add_library(THERMAL OBJECT
  "thermal_isothermal.f90"
  "thermal_adiabatic.f90"
  "thermal_conduction.f90")
add_dependencies(THERMAL CRYSTALLITE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:THERMAL>)

add_library(DAMASK_ENGINE OBJECT "homogenization.f90")
add_dependencies(DAMASK_ENGINE THERMAL DAMAGE VACANCYFLUX POROSITY HYDROGENFLUX HOMOGENIZATION)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DAMASK_ENGINE>)

add_library(DAMASK_CPFE OBJECT "CPFEM2.f90")
add_dependencies(DAMASK_CPFE DAMASK_ENGINE)
list(APPEND OBJECTFILES $<TARGET_OBJECTS:DAMASK_CPFE>)

if (PROJECT_NAME STREQUAL "DAMASK_spectral")
      add_library(SPECTRAL_UTILITIES OBJECT "spectral_utilities.f90")
      add_dependencies(SPECTRAL_UTILITIES DAMASK_CPFE)
      list(APPEND OBJECTFILES $<TARGET_OBJECTS:SPECTRAL_UTILITIES>)

      add_library(SPECTRAL_SOLVER OBJECT
      "spectral_thermal.f90"
      "spectral_damage.f90"
      "spectral_mech_Polarisation.f90"
      "spectral_mech_Basic.f90")
      add_dependencies(SPECTRAL_SOLVER SPECTRAL_UTILITIES)
      list(APPEND OBJECTFILES $<TARGET_OBJECTS:SPECTRAL_SOLVER>)
      if(NOT CMAKE_BUILD_TYPE STREQUAL "SYNTAXONLY")
        add_executable(DAMASK_spectral "DAMASK_spectral.f90" ${OBJECTFILES})
      else()
        add_library(DAMASK_spectral OBJECT "DAMASK_spectral.f90")
      endif()
      add_dependencies(DAMASK_spectral SPECTRAL_SOLVER)
elseif (PROJECT_NAME STREQUAL "DAMASK_FEM")
      add_library(FEM_UTILITIES OBJECT "FEM_utilities.f90")
      add_dependencies(FEM_UTILITIES DAMASK_CPFE)

      add_library(FEM_SOLVER OBJECT
      "FEM_hydrogenflux.f90"
      "FEM_porosity.f90"
      "FEM_vacancyflux.f90"
      "FEM_damage.f90"
      "FEM_thermal.f90"
      "FEM_mech.f90")
      add_dependencies(FEM_SOLVER FEM_UTILITIES)

      add_executable(DAMASK_FEM "DAMASK_FEM_driver.f90")
      add_dependencies(DAMASK_FEM FEM_SOLVER)
endif()
