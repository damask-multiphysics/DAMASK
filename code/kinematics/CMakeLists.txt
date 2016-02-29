# group sources
set (KINEMATICS   "kinematics_cleavage_opening"
                  "kinematics_slipplane_opening"
                  "kinematics_thermal_expansion"
                  "kinematics_vacancy_strain"
                  "kinematics_hydrogen_strain"
    )

# compile module and cumulatively link the
# compiled libraries
foreach (p ${KINEMATICS})
    add_library (${p} "${p}.f90")
    target_link_libraries(${p} DAMASK_DRIVERS)
    add_library (DAMASK_DRIVERS ALIAS ${p})
endforeach (p)