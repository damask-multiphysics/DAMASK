# group sources
set (KINEMATICS   "kinematics_cleavage_opening"
                  "kinematics_slipplane_opening"
                  "kinematics_thermal_expansion"
                  "kinematics_vacancy_strain"
                  "kinematics_hydrogen_strain"
    )

# compile kinamtic modules
foreach (p ${KINEMATICS})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)
