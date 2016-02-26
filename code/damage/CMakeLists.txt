# group sources
set (DAMAGE   "damage_none"
              "damage_local"
              "damage_nonlocal"
    )

# compile damage module
foreach (p ${DAMAGE})
    add_library (${p} "${p}.f90")
endforeach (p)

# set libraries for linking
foreach (p ${DAMAGE})
    set (DAMASK_LIB ${DAMASK_LIB} ${P})
endforeach (p)