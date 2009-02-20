set(SRC
    calc_b_matrix.F90
    calc_prim_int.F90
    deloc2xcart.F90
    fred2fdeloc.F90
    hessinit.F90
    hessupdt.F90
    make_angles.F90
    make_bonds.F90
    make_dihedrals.F90
    make_prim_internals.F90
    xcart2deloc.F90
    xcart2deloc_fixb.F90
    xfpack.F90
)

include_directories(../defs)
add_library(16geomoptim ${SRC})
add_dependencies(16geomoptim defs)
