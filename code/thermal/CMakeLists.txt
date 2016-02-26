# group source related themal module
set (THERMAL  "thermal_isothermal"
              "thermal_adiabatic"
              "thermal_conduction"
    )

# compiler theraml module
foreach (p ${THERMAL})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries/modules for linking
foreach (p ${THERMAL})
    set (DAMASK_LIB ${DAMASK_LIB} ${p})
endforeach (p)