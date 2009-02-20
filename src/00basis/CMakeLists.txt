set(SRC
    leave_myproc.F90
    initpapichoice.F90
    papi_init.F90
    timab.F90
    time_accu.F90
    timein.F90
    wrtout_myproc.F90
    )

add_executable(test ${SRC})
