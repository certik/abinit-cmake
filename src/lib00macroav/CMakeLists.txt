set(SRC
    io.f
    iorho.f
    paste.f
    numeric.f
    surpla.f
    thetaft.f
    volcel.f
)

include_directories(../defs)
add_library(lib00macroav ${SRC})
add_dependencies(lib00macroav defs)
