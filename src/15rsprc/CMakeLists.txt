set(SRC
    mod_prc_memory.F90
    lavnl.F90
    ladielmt.F90
    moddiel_csrb.F90
    prctfw3.F90
    prcrskerker2.F90
    prcrskerker1.F90
    prctfvw2.F90
    prctfvw1.F90
)

include_directories(../defs)
add_library(15rsprc ${SRC})
add_dependencies(15rsprc defs)
