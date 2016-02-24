# group sources
set (POROSITY "porosity_none"
              "porosity_phasefield"
    )

# compile porosity modules
foreach (p ${POROSITY})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${POROSITY})
    set (AUX_LIB ${AUX_LIB} ${p})
endforeach (p)