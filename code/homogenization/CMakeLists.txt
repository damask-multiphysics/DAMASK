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
