set(SRC
    chkpawovlp.F90
    clp.F90
    dbeta.F90
    gaunt.F90
    gipaw_aug_fields.F90
    gipaw_j_dia_aug.F90
    indgrid.F90
    initang.F90
    initrhoij.F90
    initylmr.F90
    int_ang.F90
    make_cs_dia.F90
    make_efg_paw.F90
    make_fc_paw.F90
    mkeuler.F90
    nhatgrid.F90
    optics_paw.F90
    partial_dos_fractions_paw.F90
    pawaccrhoij.F90
    pawalloc.F90
    pawdenpot.F90
    pawdij.F90
    pawgrnl.F90
    pawgylm.F90
    pawgylmg.F90
    pawinit.F90
    pawlsylm.F90
    pawmknhat.F90
    pawmknhat3.F90
    pawmkrhoij.F90
    pawmkrhoij3.F90
    pawnabla_init.F90
    pawpolev.F90
    pawprt.F90
    pawpupot.F90
    pawpuxinit.F90
    pawshpfun.F90
    pawsushat.F90
    pawuenergy.F90
    pawxc.F90
    pawxcm.F90
    pawxcsph.F90
    pawxenergy.F90
    pawxpot.F90
    permutations.F90
    phim.F90
    realgaunt.F90
    setnoccmmp.F90
    setsymrhoij.F90
    simple_j_dia.F90
    smatrix_paw.F90
    smatrix_pawinit.F90
    spline_paw_fncs.F90
    symdij.F90
    symrhoij.F90
    symrhoij3.F90
    transgrid.F90
)

include_directories(../defs)
add_library(13paw ${SRC})
add_dependencies(13paw defs)
