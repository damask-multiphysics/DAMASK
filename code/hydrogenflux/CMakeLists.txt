# group sources
set (HYDROGENFLUX   "hydrogenflux_isoconc"
                    "hydrogenflux_cahnhilliard"
    )

# compile hydrogenflux modules
foreach (p ${HYDROGENFLUX})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${HYDROGENFLUX})
    set (AUX_LIB ${AUX_LIB} ${p})
endforeach (p)