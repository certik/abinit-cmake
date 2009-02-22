set(SRC
    bldlwf.F90
    chkilwf.F90
    invars7w.F90
    overlap_ph.F90
    readeig.F90
    rwwan.F90
    secinit.F90
    shellin.F90
    wanvec.F90
    zmnbld.F90
)

include_directories(../defs)
add_library(17lwf ${SRC})
add_dependencies(17lwf defs)
