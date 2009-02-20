set(SRC
    bracketing.F90
    brent.F90
    cgpr.F90
    dotproduct.F90
    linmin.F90
)

include_directories(../defs)
add_library(lib01cg ${SRC})
add_dependencies(lib01cg defs)
