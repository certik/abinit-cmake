project(abinit)

cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "-ffree-line-length-256")

if(NOT DEFINED WITH_MPI)
    set(WITH_MPI NO)
endif(NOT DEFINED WITH_MPI)

if(WITH_MPI)
    set(CMAKE_Fortran_COMPILER mpif90)
    include_directories(/usr/lib/openmpi/lib)
endif(WITH_MPI)

add_definitions(-DHAVE_CONFIG_H)

add_subdirectory(src)

message("\nConfiguration results (these options can be set on the commandline)")
message("---------------------")
message("Build with MPI: ${WITH_MPI}")
message("\n")
