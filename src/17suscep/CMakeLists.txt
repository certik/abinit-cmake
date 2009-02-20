set(SRC
    inwffil3.F90
    k_rpa.F90
    kxc_pgg.F90
    kxc_alda.F90
    kxc_eok.F90
    klocal.F90
    dyson_de.F90
    dyson_gl.F90
    dyson_ls.F90
    dyson_sc.F90
    acfd_dyson.F90
    acfd_intexact.F90
    getlambda.F90
    getfreqsus.F90
    get_g_tiny.F90
    geteexc_cc.F90
    geteexc_uc.F90
    get_susd_null.F90
    prtsusd.F90
    suscep.F90
    suscep_dyn.F90
    suscep_stat.F90
    susk.F90
    susk_dyn.F90
    susk_dyn_pgg.F90
    suskmm.F90
    suskmm_dyn.F90
    xcacfd.F90
)

include_directories(../defs)
add_library(17suscep ${SRC})
add_dependencies(17suscep defs)
