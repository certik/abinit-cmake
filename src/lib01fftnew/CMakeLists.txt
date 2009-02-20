set(SRC
    accrho.F90
    addrho.F90
    applypot.F90
    back.F90
    back_wf.F90
    ctrig.F90
    fftstp.F90
    fill_cent.F90
    fill.F90
    forw.F90
    forw_wf.F90
    indirect_parallel_Fourier.F90
    mpiswitch.F90
    mpiswitch_cent.F90
    multpot.F90
    scramble.F90
    switch.F90
    switch_cent.F90
    switchreal.F90
    switchreal_cent.F90
    unfill.F90
    unfill_cent.F90
    unmpiswitch.F90
    unmpiswitch_cent.F90
    unscramble.F90
    unswitch.F90
    unswitch_cent.F90
    unswitchreal.F90
    unswitchreal_cent.F90
)

include_directories(../defs)
add_library(lib01fftnew ${SRC})
add_dependencies(lib01fftnew defs)
