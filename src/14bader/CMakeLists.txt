set(SRC
    addout.F90
    adini.F90
    aim_follow.F90
    consist.F90
    cpdrv.F90
    critic.F90
    critics.F90
    defad.F90
    drvaim.F90
    evspln.F90
    graph.F90
    initaim.F90
    inpar.F90
    inspln.F90
    integrho.F90
    integvol.F90
    onestep.F90
    plint.F90
    rsurf.F90
    surf.F90
    vgh_rho.F90
    vnorm.F90
)

include_directories(../defs)
add_library(14bader ${SRC})
add_dependencies(14bader defs)
