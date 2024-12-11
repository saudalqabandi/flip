cmake_minimum_required(VERSION 3.5.0)
project(flip VERSION 1.0.0 LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# set(CMAKE_Fortran_FLAGS "-O3 -march=native -Wall -Wextra -Wconversion -Warray-temporaries -fcheck=all -fbacktrace -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS " -g -O3 -march=native")

set(SOURCES
    src/main.f90
    src/config.f90
    src/utils.f90
    src/particle.f90
    src/potential.f90
    src/accumalators.f90
    src/moves.f90
)

add_executable(flip ${SOURCES})


