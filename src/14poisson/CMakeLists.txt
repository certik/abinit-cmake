set(SRC
    psolver_kernel.F90
    psolver_hartree.F90
    psolver_rhohxc.F90
)

include_directories(../defs)
add_library(14poisson ${SRC})
add_dependencies(14poisson defs)
