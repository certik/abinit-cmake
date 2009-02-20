set(SRC
    pstate.F90
)

include_directories(../defs)
add_library(21paral_md ${SRC})
add_dependencies(21paral_md defs)
