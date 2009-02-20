set(SRC
    dens_in_sph.F90
    dos_hdr_write.F90
    getnel.F90
    get_dos_1band.F90
    get_dos_1band_m.F90
    get_fsurf_1band.F90
    get_tetra_weight.F90
    init_bess_spl.F90
    init_ylm_spl.F90
    newocc.F90
    occeig.F90
    pareigocc.F90
    partial_dos_fractions.F90
    prtxcfermsurf.F90
    printbxsf.F90
    recip_ylm.F90
    simpson_int.F90
    sphericaldens.F90
    tetrahedron.F90
)

include_directories(../defs)
add_library(14occeig ${SRC})
add_dependencies(14occeig defs)
