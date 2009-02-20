set(SRC
    clnmpi_band.F90
    clnmpi_fft.F90
    clnmpi_gs.F90
    distrb2.F90
    herald.F90
    initmpi_band.F90
    initmpi_grid.F90
    initmpi_fft.F90
    initmpi_gs.F90
    initmpi_respfn.F90
    initmpi_seq.F90
    leave_new.F90
    leave_test.F90
    pre_scatter.F90
    scalapack.F90
    split_work.F90
    split_work2.F90
    wrtout.F90
)

include_directories(../defs)
add_library(01manage_mpi ${SRC})
add_dependencies(01manage_mpi defs)
