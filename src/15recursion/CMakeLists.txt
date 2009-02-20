set(SRC
    getngrec.F90
    recursion.F90
    entropyrec.F90
    free_energyrec.F90
    fermisolverec.F90
    green_kernel.F90
    vtorhorec.F90
)

include_directories(../defs)
add_library(15recursion ${SRC})
add_dependencies(15recursion defs)
