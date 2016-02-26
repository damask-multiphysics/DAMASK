# group source for sepctral solver driver
set (SPECTRAL "spectral_damage"
              "spectral_mech_AL"
              "spectral_mech_Basic"
              "spectral_mech_Polarisation"
              "spectral_thermal"
              "spectral_utilities"
    )

# compile spectral solver driver module
foreach (p ${SPECTRAL})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${SPECTRAL})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)