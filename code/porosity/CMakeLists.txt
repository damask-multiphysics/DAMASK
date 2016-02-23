# group sources
set (POROSITY "porosity_none"
              "porosity_phasefield"
    )

# compile porosity modules
foreach (p ${POROSITY})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)
