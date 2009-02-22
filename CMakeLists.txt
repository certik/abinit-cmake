project(abinit)

cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
enable_language(Fortran)

if(DEFINED ENV{FCFLAGS})
    set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
else(DEFINED ENV{FCFLAGS})
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-256")
endif(DEFINED ENV{FCFLAGS})

if(NOT DEFINED WITH_MPI)
    set(WITH_MPI NO)
endif(NOT DEFINED WITH_MPI)

if(WITH_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_Fortran_COMPILER ${MPI_COMPILER})
endif(WITH_MPI)

add_definitions(-DHAVE_CONFIG_H)

add_subdirectory(src)

message("\nConfiguration results (these options can be set on the commandline)")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
message("Build with MPI: ${WITH_MPI}")
message("\n")
