# group sources
set (HYDROGENFLUX   "hydrogenflux_isoconc"
                    "hydrogenflux_cahnhilliard"
    )

# compile hydrogenflux modules
foreach (p ${HYDROGENFLUX})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${HYDROGENFLUX})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)