# The dependency detection in CMake is not functioning for Fortran
# !!! EXPLICIT DEPENDENCY DECLARATION !!!
###################################################################################################
# fno-range-check:  Disable range checking on results of simplification of constant expressions during compilation
#      -->   allows the definition of DAMASK_NaN
#-fall-intrinsics:  all intrinsic procedures (including the GNU-specific extensions) are accepted. -Wintrinsics-std will be ignored
#                   and no user-defined procedure with the same name as any intrinsic will be called except when it is explicitly declared external
#      -->   allows the use of 'isnan'
#-fno-fast-math:
#      -->   otherwise, when setting -ffast-math, isnan always evaluates to false (I would call it a bug)
###################################################################################################
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
      set (SPEICAL_FLAGS "")
elseif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
      set (SPEICAL_FLAGS "-fno-range-check -fall-intrinsics -fno-fast-math")
else (Fortran_COMPILER_NAME MATCHES "ifort.*")
      message (FATAL_ERROR "Unknown Fortrna compiler.")
endif (Fortran_COMPILER_NAME MATCHES "ifort.*")

set (CMAKE_Fortran_FLAGS_RELEASE "${DAMASK_Fortran_FLAGS_RELEASE} ${SPEICAL_FLAGS}")

add_library(DAMASK_prec "prec.f90")

set (CMAKE_Fortran_FLAGS_RELEASE "${DAMASK_Fortran_FLAGS_RELEASE}")

###################################################################################################
# -diag-disable 7410 should disable warning about directory statement in inquire function,
#  but does not work. hence the other 2 statements
###################################################################################################
if (SPECTRAL)
      if (Fortran_COMPILER_NAME MATCHES "ifort.*")
            set (SPEICAL_FLAGS "-diag-remark 7410 -stand none -warn nostderrors")
      elseif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
            set (SPEICAL_FLAGS "-fall-intrinsics")
      else (Fortran_COMPILER_NAME MATCHES "ifort.*")
            message (FATAL_ERROR "Unknown Fortrna compiler.")
      endif (Fortran_COMPILER_NAME MATCHES "ifort.*")

      set (CMAKE_Fortran_FLAGS_RELEASE "${DAMASK_Fortran_FLAGS_RELEASE} ${SPEICAL_FLAGS}")

      add_library(DAMASK_interface "spectral_interface.f90")
      target_link_libraries(DAMASK_interface DAMASK_prec)

      set (CMAKE_Fortran_FLAGS_RELEASE "${DAMASK_Fortran_FLAGS_RELEASE}")
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