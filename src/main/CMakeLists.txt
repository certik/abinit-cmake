include_directories(../defs)
add_executable(abinis abinit.F90)
add_dependencies(abinis 00basis 01manage_mpi)
target_link_libraries(abinis
    defs 00basis 01manage_mpi
    #11util 11contract 12spacepar 12parser 12geometry
    #13psp 12ffts 12nlstrain 13nonlocal 13xml 13xc 13recipspace 13paw 13io_mpis
    #13iovars 13ionetcdf 14poissons 14wfs 14wvl_wfs 14iowfdenpot 15rsprc
    #15recursion 15common 15gw 16geomoptim 16response 17suscep 18seqpars 21rdm
    #21paral_mds
    )
