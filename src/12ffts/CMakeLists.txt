set(SRC
    ccfft.F90
    fftpac.F90
    fftw.F90
    fourdp.F90
    fourwf.F90
    indfftrisc.F90
    kgindex.F90
    sg_ctrig.F90
    sg_fft.F90
    sg_fftpad.F90
    sg_fftpx.F90
    sg_fftrisc.F90
    sg_fftx.F90
    sg_ffty.F90
    sg_fftz.F90
    sg_fourwf.F90
    sphereboundary.F90
    sphere.F90
    sphere_fft.F90
    zerosym.F90
)

include_directories(../defs)
add_library(12ffts ${SRC})
add_dependencies(12ffts defs)
