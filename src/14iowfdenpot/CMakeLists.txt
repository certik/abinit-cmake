set(SRC
    calcdensph.F90
    fappnd.F90
    hdr_check.F90
    hdr_clean.F90
    hdr_init.F90
    hdr_update.F90
    initwf.F90
    ioarr.F90
    mk_hdr_check_fmt.F90
    out1dm.F90
    outwant.F90
    randac.F90
    rdnpw.F90
    rdkss.F90
    read_wfrspa.F90
    testlda.F90
    WffReadEigK.F90
    WffReadSkipK.F90
)

include_directories(../defs)
add_library(14iowfdenpot ${SRC})
add_dependencies(14iowfdenpot defs)
