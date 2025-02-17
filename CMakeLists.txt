cmake_minimum_required(VERSION 3.10)
project(flip VERSION 1.0.0 LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -Wextra -Wconversion -Warray-temporaries -fcheck=all -fbacktrace -ffree-line-length-none")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
endif()

message(STATUS "Build flags: ${CMAKE_Fortran_FLAGS}")
file(GLOB SOURCES "src/*.f90")

find_package(OpenMP REQUIRED)

add_executable(flip ${SOURCES})
target_link_libraries(flip PRIVATE OpenMP::OpenMP_Fortran)

set(OMP_NUM_THREADS 10)
set_property(TARGET flip PROPERTY ENVIRONMENT "OMP_NUM_THREADS=${OMP_NUM_THREADS}")

find_package(pFUnit REQUIRED)
enable_testing()
add_subdirectory(src)

add_pfunit_ctest (flip_random
  TEST_SOURCES tests/test_random.pf
  LINK_LIBRARIES flipLib
  )

  





