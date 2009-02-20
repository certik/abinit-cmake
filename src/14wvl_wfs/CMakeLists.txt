set(SRC
    wvl_setngfft.F90
    wvl_free_type.F90
    wvl_init_type_wfs.F90
    wvl_init_type_proj.F90
    wvl_tail_corrections.F90
    wvl_nl_gradient.F90
    wvl_rwwf.F90
)

include_directories(../defs)
add_library(14wvl_wfs ${SRC})
add_dependencies(14wvl_wfs defs)
