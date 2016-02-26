# group sources
set (POROSITY "porosity_none"
              "porosity_phasefield"
    )

# compile porosity modules
foreach (p ${POROSITY})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${POROSITY})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)