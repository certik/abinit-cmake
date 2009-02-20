set(SRC
    contistr01.F90
    contistr03.F90
    contistr12.F90
    contstr21.F90
    contstr23.F90
    contstr25.F90
    contstr25a.F90
    contstr26.F90
)

include_directories(../defs)
add_library(12nlstrain ${SRC})
add_dependencies(12nlstrain defs)
