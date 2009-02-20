set(SRC
    bound.F90
    getcut.F90
    getkgrid.F90
    getkpgnorm.F90
    getmpw.F90
    getng.F90
    getph.F90
    getwtk.F90
    get_irredg.F90
    get_full_gsphere.F90
    get_full_kgrid.F90
    get_tetra.F90
    laplacian.F90
    initylmg.F90
    irrzg.F90
    kpgio.F90
    kpgsph.F90
    merge_kgirr.F90
    mkkin.F90
    pmat2cart.F90
    setsym.F90
    smpbz.F90
    symkchk.F90
    symkpt.F90
    symg.F90
    symq3.F90
    testkgrid.F90
)

include_directories(../defs)
add_library(13recipspace ${SRC})
add_dependencies(13recipspace defs)
