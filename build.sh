#!/bin/bash

if [ "$1" != "Release" ] && [ "$1" != "Debug" ]; then
  echo "Argument must be either 'release' or 'debug'"
  exit 1
fi

#build pfunit
# cd external/pFUnit

# if [ "$2" == "-r" ]; then
#   if [ -d "build" ]; then
#     rm -rf build
#   fi
# fi

# if [ ! -d "build" ]; then
#   mkdir build
# fi
# cd build
# cmake .. -DSKIP_MPI=YES -DSKIP_OPENMP=YES
# make install -j

# cd ../../..

# Configure the project
cmake -DCMAKE_BUILD_TYPE=$1 -S . -B build -DCMAKE_INSTALL_PREFIX=./external/pFUnit/build/installed

# Build the project
cmake --build build
