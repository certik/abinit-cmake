set(SRC
    instrng.F90
    inread.F90
    inreplsp.F90
    intagm.F90
    inarray.F90
    incomprs.F90
)

include_directories(../defs)
add_library(12parser ${SRC})
add_dependencies(12parser defs)
