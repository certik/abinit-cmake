set(SRC
    calc_lifetime.F90
    calc_xc_ep.F90
    drivexc.F90
    hartre.F90
    invcb.F90
    lifetime_bn.F90
    lifetime_psn.F90
    lifetime_rpa.F90
    mkcore.F90
    mkdenpos.F90
    mkvxc3.F90
    mkvxcgga3.F90
    phase.F90
    rhohxc.F90
    xc_kernel.F90
    xcden.F90
    xce_ap.F90
    xcepsn_tcdft.F90
    xchcth.F90
    xchelu.F90
    xclb.F90
    xcmult.F90
    xcp_ap.F90
    xcpot.F90
    xcpbe.F90
    xcppsn_tcdft.F90
    xcpzca.F90
    xcspol.F90
    xctetr.F90
    xcwign.F90
    xcxalp.F90
)

include_directories(../defs)
add_library(13xc ${SRC})
add_dependencies(13xc defs)
