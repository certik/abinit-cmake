set(SRC
    berryphase_new.F90
    cgwf.F90
    inwffil.F90
    iofn1.F90
    iofn2.F90
    loper3.F90
    lobpcgccwf.F90
    lobpcgccIIwf.F90
    lobpcgwf.F90
    lobpcgIIwf.F90
    mv_3dte.F90
    outkss.F90
    outwf.F90
    prep_nonlop.F90
    respfn.F90
    subdiago.F90
    tddft.F90
    vtorho.F90
    vtowfk.F90
    wfsinp.F90
    wvl_wfsinp_disk.F90
    wvl_wfsinp_scratch.F90
    wvl_wfsinp_reformat.F90
    wvl_vtorho.F90
)

include_directories(../defs)
add_library(18seqpar ${SRC})
add_dependencies(18seqpar defs)
