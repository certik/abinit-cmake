set(SRC
    xallgather_mpi.F90
    xallgatherv_mpi.F90
    xalltoallv_mpi.F90
    xalltoall_mpi.F90
    xcast_mpi.F90
    xdef_comm.F90
    xdefineoff.F90
    xexch_mpi.F90
    xmax_mpi.F90
    xmin_mpi.F90
    xsum_master.F90
    xsum_mpi.F90
    xderiveread.F90
    xderivereadval.F90
    xderivewrite.F90
    xderivewriteval.F90
)

include_directories(../defs)
add_library(lib01hidempi ${SRC})
add_dependencies(lib01hidempi defs)
