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

# compile module and cumulatively link the
# compiled libraries
add_library (DAMASK_PLASTIC "plastic_dislotwin.f90"
                            "plastic_disloUCLA.f90"
                            "plastic_isotropic.f90"
                            "plastic_j2.f90"
                            "plastic_phenopowerlaw.f90"
                            "plastic_titanmod.f90"
                            "plastic_nonlocal.f90"
                            "plastic_none.f90"
                            "plastic_phenoplus.f90")
target_link_libraries(DAMASK_PLASTIC DAMASK_DRIVERS)
# foreach (p ${PLASTIC})
#     add_library (${p} "${p}.f90")
#     target_link_libraries(${p} DAMASK_DRIVERS)
#     add_library (DAMASK_DRIVERS ALIAS ${p})
# endforeach (p)