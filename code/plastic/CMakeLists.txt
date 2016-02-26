# group sources
set (PLASTIC  "plastic_dislotwin"
              "plastic_disloUCLA"
              "plastic_isotropic"
              "plastic_j2"
              "plastic_phenopowerlaw"
              "plastic_titanmod"
              "plastic_nonlocal"
              "plastic_none"
              "plastic_phenoplus"
    )

# compile module
foreach (p ${PLASTIC})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${PLASTIC})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)