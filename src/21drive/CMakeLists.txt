set(SRC
    afterscfloop.F90
    brdmin.F90
    delocint.F90
    diisrelax.F90
    driver.F90
    elpolariz.F90
    gstate.F90
    isotemp.F90
    loop3dte.F90
    moldyn.F90
    move.F90
    nonlinear.F90
    outscfcv.F90
    scfcv.F90
    screening.F90
    sigma.F90
    testfi.F90
    timana.F90
    wannier.F90
)

include_directories(../defs)
add_library(21drive ${SRC})
add_dependencies(21drive defs)
