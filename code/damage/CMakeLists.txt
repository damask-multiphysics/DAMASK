# group sources
set (DAMAGE   "damage_none"
              "damage_local"
              "damage_nonlocal"
    )

# compile damage module
foreach (p ${DAMAGE})
    add_library (${p} MODULE "${p}.f90")
endforeach (p)
