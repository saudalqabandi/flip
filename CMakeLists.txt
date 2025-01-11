cmake_minimum_required(VERSION 3.5.0)
project(flip VERSION 1.0.0 LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

set(CMAKE_Fortran_FLAGS "-g -Wall -Wextra -Wconversion -Warray-temporaries -fcheck=all -fbacktrace -ffree-line-length-none")
# set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")

file(GLOB SOURCES "src/*.f90")


add_executable(flip ${SOURCES})


