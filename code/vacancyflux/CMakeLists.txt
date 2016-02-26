# group source file
set (VACANCYFLUX    "vacancyflux_isoconc"
                    "vacancyflux_isochempot"
                    "vacancyflux_cahnhilliard"
    )

# compiler as module
foreach (p ${VACANCYFLUX})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${VACANCYFLUX})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)