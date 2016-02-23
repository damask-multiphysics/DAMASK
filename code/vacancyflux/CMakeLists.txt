# group source file
set (VACANCYFLUX    "vacancyflux_isoconc"
                    "vacancyflux_isochempot"
                    "vacancyflux_cahnhilliard"
    )

# compiler as module
foreach (p ${VACANCYFLUX})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)
