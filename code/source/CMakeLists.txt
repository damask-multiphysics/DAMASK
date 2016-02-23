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

# compile modules for source
foreach (p ${SOURCE})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)
