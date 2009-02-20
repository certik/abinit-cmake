set(SRC
    bestwfs.F90
    envlop.F90
    fxphas.F90
    getghc.F90
    listkk.F90
    pw_orthon.F90
    precon.F90
    precon2.F90
    prep_fourwf.F90
    prep_getghc.F90
    projbd.F90
    sdirot.F90
    wfconv.F90
    prep_kg_sym_do.F90
    prep_wavef_sym_do.F90
    prep_wavef_sym_undo.F90
    prep_index_wavef_bandpp.F90
    zprecon3.F90
)

include_directories(../defs)
add_library(14wfs ${SRC})
add_dependencies(14wfs defs)
