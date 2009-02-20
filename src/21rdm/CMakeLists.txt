set(SRC
    crho.F90
    cvxclda.F90
    fftwfn.F90
    lattice.F90
    occred.F90
    old_setmesh.F90
    rdm.F90
    setup_G_rotation_old.F90
)

include_directories(../defs)
add_library(21rdm ${SRC})
add_dependencies(21rdm defs)
