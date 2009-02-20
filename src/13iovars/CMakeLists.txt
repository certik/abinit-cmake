set(SRC
    chkdpr.F90
    chkinp.F90
    chkint.F90
    chkneu.F90
    chkvars.F90
    getdim_nloc.F90
    indefo.F90
    ingeo.F90
    ingeobld.F90
    inkpts.F90
    invacuum.F90
    invars0.F90
    invars1.F90
    invars1m.F90
    invars2.F90
    invars2m.F90
    memana.F90
    memorf.F90
    memory.F90
    outqmc.F90
    outvars.F90
    outvar1.F90
    out_geometry_xml.F90
    out_resultsgs_xml.F90
    parsefile.F90
    prtocc.F90
    prttagm.F90
    wvl_memory.F90
)

include_directories(../defs)
add_library(13iovars ${SRC})
add_dependencies(13iovars defs)
