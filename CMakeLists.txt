project(abinit)

cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "-ffree-line-length-256")

add_subdirectory(src)
