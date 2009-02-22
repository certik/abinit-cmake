set(SRC
    hirsh.F90
    kptindex.F90
    lineint.F90
    localorb_S.F90
    normalize.F90
    overlap_wf.F90
    planeint.F90
    pointint.F90
    recip.F90
    reduce.F90
    rrho.F90
    rtau.F90
    vdot.F90
    volumeint.F90
    wffile.F90
    wfread.F90
)

include_directories(../defs)
add_library(19cut3d ${SRC})
add_dependencies(19cut3d defs)
