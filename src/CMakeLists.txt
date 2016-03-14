# The dependency detection in CMake is not functioning for Fortran
# !!! EXPLICIT DEPENDENCY DECLARATION !!!

# SOME FILES REQUIRE SPECIAL FLAGS
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    # -diag-disable 7410 should disable warning about directory statement
    # in inquire function, but does not work. hence the other 2 statements
    SET_SOURCE_FILES_PROPERTIES( "spectral_interface.f90" PROPERTIES
                                    COMPILE_FLAGS "-diag-remark 7410 -stand none -warn nostderrors" )
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    # fno-range-check:  Disable range checking on results of simplification of constant expressions during compilation
    #      -->   allows the definition of DAMASK_NaN
    #-fall-intrinsics:  all intrinsic procedures (including the GNU-specific extensions) are accepted. -Wintrinsics-std will be ignored
    #                   and no user-defined procedure with the same name as any intrinsic will be called except when it is explicitly declared external
    #      -->   allows the use of 'isnan'
    #-fno-fast-math:
    #      -->   otherwise, when setting -ffast-math, isnan always evaluates to false (I would call it a bug)
    SET_SOURCE_FILES_PROPERTIES( "prec.f90" PROPERTIES
                                    COMPILE_FLAGS "-fno-range-check -fall-intrinsics -fno-fast-math")
    #-fall-intrinsics:  all intrinsic procedures (including the GNU-specific extensions) are accepted. -Wintrinsics-std will be ignored
    #                   and no user-defined procedure with the same name as any intrinsic will be called except when it is explicitly declared external
    #      -->   allows the use of 'getcwd'
    SET_SOURCE_FILES_PROPERTIES( "spectral_interface.f90" PROPERTIES
                                    COMPILE_FLAGS "-fall-intrinsics")
    # long lines for interaction matrix
    SET_SOURCE_FILES_PROPERTIES( "lattice.f90" PROPERTIES
                                COMPILE_FLAGS "-ffree-line-length-240")
endif()

add_library(DAMASK_PREC "prec.f90")

if (SPECTRAL)
      add_library(DAMASK_INTERFACE "spectral_interface.f90")
elseif(FEM)
      add_library(DAMASK_INTERFACE "DAMASK_FEM_interface.f90")
endif(SPECTRAL)
target_link_libraries(DAMASK_INTERFACE DAMASK_PREC)

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
      target_link_libraries(DAMASK_MESH DAMASK_MATH)
endif(SPECTRAL)
if (FEM)
      add_library(DAMASK_FEZoo "FEZoo.f90")
      target_link_libraries(DAMASK_FEZoo DAMASK_MATH)
      add_library(DAMASK_MESH  "meshFEM.f90")
      target_link_libraries(DAMASK_MESH  DAMASK_FEZoo)
endif(FEM)


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
target_link_libraries(DAMASK_CONSTITUTIVE DAMASK_PLASTIC   )
target_link_libraries(DAMASK_CONSTITUTIVE DAMASK_KINEMATICS)
target_link_libraries(DAMASK_CONSTITUTIVE DAMASK_SOURCE    )

add_library(DAMASK_CRYSTALLITE "crystallite.f90")
target_link_libraries(DAMASK_CRYSTALLITE DAMASK_CONSTITUTIVE)

add_library(DAMASK_HOMOGENIZATION "homogenization_RGC.f90"
                                  "homogenization_isostrain.f90"
                                  "homogenization_none.f90")
target_link_libraries(DAMASK_HOMOGENIZATION DAMASK_CRYSTALLITE)

add_library(DAMASK_HYDROGENFLUX "hydrogenflux_isoconc.f90"
                                "hydrogenflux_cahnhilliard.f90")
target_link_libraries(DAMASK_HYDROGENFLUX DAMASK_CRYSTALLITE)

add_library(DAMASK_POROSITY "porosity_none.f90"
                            "porosity_phasefield.f90")
target_link_libraries(DAMASK_POROSITY DAMASK_CRYSTALLITE)

add_library(DAMASK_VACANCYFLUX "vacancyflux_isoconc.f90"
                               "vacancyflux_isochempot.f90"
                               "vacancyflux_cahnhilliard.f90")
target_link_libraries(DAMASK_VACANCYFLUX DAMASK_CRYSTALLITE)

add_library(DAMASK_DAMAGE "damage_none.f90"
                          "damage_local.f90"
                          "damage_nonlocal.f90")
target_link_libraries(DAMASK_DAMAGE DAMASK_CRYSTALLITE)

add_library(DAMASK_THERMAL "thermal_isothermal.f90"
                           "thermal_adiabatic.f90"
                           "thermal_conduction.f90")
target_link_libraries(DAMASK_THERMAL DAMASK_CRYSTALLITE)

add_library(DAMASK_ENGINE "homogenization.f90")
target_link_libraries(DAMASK_ENGINE DAMASK_THERMAL       )
target_link_libraries(DAMASK_ENGINE DAMASK_DAMAGE        )
target_link_libraries(DAMASK_ENGINE DAMASK_VACANCYFLUX   )
target_link_libraries(DAMASK_ENGINE DAMASK_POROSITY      )
target_link_libraries(DAMASK_ENGINE DAMASK_HYDROGENFLUX  )
target_link_libraries(DAMASK_ENGINE DAMASK_HOMOGENIZATION)


if (FEM)
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

      add_executable(DAMASK_FEM.exe "DAMASK_FEM_driver.f90")
      target_link_libraries(DAMASK_FEM.exe DAMASK_FEM_DRIVER)
endif(FEM)

if (SPECTRAL)
      add_library(DAMASK_CPFE "CPFEM2.f90")
      target_link_libraries(DAMASK_CPFE DAMASK_ENGINE)

      add_library(DAMASK_SPECTRAL_UTILITY spectral_utilities.f90)
      target_link_libraries(DAMASK_SPECTRAL_UTILITY DAMASK_CPFE)

      add_library(DAMASK_SPECTRAL_BASE "spectral_thermal.f90"
                                       "spectral_damage.f90")
      target_link_libraries(DAMASK_SPECTRAL_BASE DAMASK_SPECTRAL_UTILITY)

      add_library(DAMASK_SPECTRAL_MECH "spectral_mech_AL.f90"
                                       "spectral_mech_Polarisation.f90"
                                       "spectral_mech_Basic.f90")
      target_link_libraries(DAMASK_SPECTRAL_MECH DAMASK_SPECTRAL_UTILITY)

      add_library(DAMASK_EXE "DAMASK_spectral.f90")
      target_link_libraries(DAMASK_EXE DAMASK_CPFE         )
      target_link_libraries(DAMASK_EXE DAMASK_SPECTRAL_BASE)
      target_link_libraries(DAMASK_EXE DAMASK_SPECTRAL_MECH)

      add_executable(DAMASKSpectral.exe "DAMASK_spectral.f90")
      target_link_libraries(DAMASKSpectral.exe DAMASK_EXE)
endif(SPECTRAL)