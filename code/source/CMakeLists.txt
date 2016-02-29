# group source
set (SOURCE "source_thermal_dissipation"
            "source_thermal_externalheat"
            "source_damage_isoBrittle"
            "source_damage_isoDuctile"
            "source_damage_anisoBrittle"
            "source_damage_anisoDuctile"
            "source_vacancy_phenoplasticity"
            "source_vacancy_irradiation"
            "source_vacancy_thermalfluc"
    )

# compile module and cumulatively link the
# compiled libraries
foreach (p ${KINEMATICS})
    add_library (${p} "${p}.f90")
    target_link_libraries(${p} DAMASK_DRIVERS)
    add_library (DAMASK_DRIVERS ALIAS ${p})
endforeach (p)