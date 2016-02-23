# group sources
set (HYDROGENFLUX   "hydrogenflux_isoconc"
                    "hydrogenflux_cahnhilliard"
    )

# compile hydrogenflux modules
foreach (p ${HYDROGENFLUX})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)
