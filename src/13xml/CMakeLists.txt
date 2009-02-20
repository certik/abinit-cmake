set(SRC
    append_cml.F90
    append_cml2.F90
    findmarkup.F90
    getattribute.F90
    importcml.F90
    prt_cml.F90
    prt_cml2.F90
)

include_directories(../defs)
add_library(13xml ${SRC})
add_dependencies(13xml defs)
