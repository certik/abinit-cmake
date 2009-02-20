set(SRC
    ass_leg_pol.F90
    coeffs_gausslegint.F90
    dzgedi.F90
    dzgefa.F90
    intrpl.F90
    interp.F90
    j0y0.F90
    j1y1.F90
    jacobi.F90
    k0.F90
    k1.F90
    ludcmp.F90
    polyn_coeff.F90
    smooth.F90
    sort_int.F90
    sort_dp.F90
    splfit.F90
    spline.F90
    splint.F90
    uniformrandom.F90
    zgedi.F90
    zgefa.F90
)

include_directories(../defs)
add_library(lib00numeric ${SRC})
add_dependencies(lib00numeric defs)
