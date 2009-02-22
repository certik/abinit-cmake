set(ABINIT_LIBS
    21drive
    21paral_md
    21rdm
    18seqpar
    17suscep
    16response
    16geomoptim
    15gw
    15common
    15recursion
    15rsprc
    14occeig
    14iowfdenpot
    14wvl_wfs
    14wfs
    14poisson
    13ionetcdf
    13iovars
    13io_mpi
    13paw
    13recipspace
    13xc
    13xml
    13nonlocal
    12nlstrain
    12ffts
    13psp
    12geometry
    12parser
    12spacepar
    11contract
    11util
    01manage_mpi
    00basis
    lib01cg
    lib01hidempi
    lib01fftnew
    lib00numeric
    defs
    )
include_directories(../defs)
find_package(LAPACK REQUIRED)

add_executable(abinis abinit.F90)
add_dependencies(abinis ${ABINIT_LIBS})
target_link_libraries(abinis ${ABINIT_LIBS} ${LAPACK_LIBRARIES})

add_executable(mrgddb mrgddb.F90)
target_link_libraries(mrgddb 17ddb 16response 13io_mpi 11util 01manage_mpi
    00basis lib01hidempi defs)

add_executable(mrggkk mrggkk.F90)
target_link_libraries(mrggkk 14iowfdenpot 13io_mpi 11util 01manage_mpi 00basis
    lib01hidempi defs)

add_executable(mrgscr mrgscr.F90)
target_link_libraries(mrgscr 15gw 14iowfdenpot 13io_mpi 12geometry 11util
    01manage_mpi 00basis lib01hidempi lib00numeric defs)

add_executable(anaddb anaddb.F90)
target_link_libraries(anaddb  17ddb 16response 15common 14occeig 14wfs
    14iowfdenpot 13recipspace 13nonlocal 13io_mpi 12geometry 12parser 11util
    01manage_mpi lib00numeric 00basis lib01hidempi lib00numeric defs
    ${LAPACK_LIBRARIES}
    )

add_executable(aim aim.F90)
target_link_libraries(aim 14bader 14iowfdenpot 13io_mpi 12parser 12geometry
    11util 01manage_mpi lib00numeric 00basis lib01hidempi defs)

add_executable(cut3d cut3d.F90)
target_link_libraries(cut3d 19cut3d 15common 14occeig 14iowfdenpot 13recipspace
    13nonlocal 13ionetcdf 13io_mpi 12geometry 12ffts 12spacepar 12parser lib01fftnew
    lib01hidempi 11contract 11util 01manage_mpi lib00numeric 00basis defs
    ${LAPACK_LIBRARIES}
    )

add_executable(macroave macroave.F90)
target_link_libraries(macroave 14iowfdenpot 13io_mpi 11util 01manage_mpi
    lib00macroav 00basis defs)
