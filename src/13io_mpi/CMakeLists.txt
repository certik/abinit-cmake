set(SRC
    chkexi.F90
    handle_ncerr.F90
    hdr_io.F90
    hdr_io_netcdf.F90
    hdr_skip.F90
    outxfhist.F90
    rwwf.F90
    wffclose.F90
    wffdelete.F90
    wffkg.F90
    wffoffset.F90
    wffopen.F90
    wffreaddatarec.F90
    wffreadnpwrec.F90
    wffreadskiprec.F90
    wffwritecg.F90
    wffwritedatarec.F90
    wffwritedatarecint.F90
    wffwritenpwrec.F90
)

if(WITH_MPI)
    set(SRC ${SRC} hdr_comm.F90)
endif(WITH_MPI)

include_directories(../defs)
add_library(13io_mpi ${SRC})
add_dependencies(13io_mpi defs)
