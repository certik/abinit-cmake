set(SRC
    abi_etsf_init.F90
    abi_etsf_electrons_put.F90
    abi_etsf_geo_put.F90
    handle_err_netcdf.F90
    hdr_io_etsf.F90
    ini_wf_etsf.F90
    ini_wf_netcdf.F90
    write_header_moldynnetcdf.F90
    write_moldynvaluenetcdf.F90
)

include_directories(../defs)
add_library(13ionetcdf ${SRC})
add_dependencies(13ionetcdf defs)
