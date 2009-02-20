include_directories(../defs)
add_executable(abinis abinit.F90)
add_dependencies(abinis 00basis)
