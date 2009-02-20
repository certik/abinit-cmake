set(SRC
    contract_dp_ge_val.F90
    contract_int_ge_val.F90
    contract_int_le_val.F90
    contract_int_list.F90
)

include_directories(../defs)
add_library(11contract ${SRC})
add_dependencies(11contract defs)
