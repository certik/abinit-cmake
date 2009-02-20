set(SRC
    dotprodm_v.F90
    dotprodm_vn.F90
    dotprod_g.F90
    dotprod_v.F90
    dotprod_vn.F90
    matrixelmt_g.F90
    mean_fftr.F90
    meanvalue_g.F90
    overlap_g.F90
    sqnormm_v.F90
    sqnorm_g.F90
    sqnorm_v.F90
)

include_directories(../defs)
add_library(12spacepar ${SRC})
add_dependencies(12spacepar defs)
