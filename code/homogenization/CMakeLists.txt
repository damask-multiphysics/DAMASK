# group sources
set (HOMOGENIZATION   "homogenization"
                      "homogenization_RGC"
                      "homogenization_isostrain"
                      "homogenization_none"
    )

# compile modules
foreach (p ${HOMOGENIZATION})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${HOMOGENIZATION})
    set (AUX_LIB ${AUX_LIB} ${p})
endforeach (p)