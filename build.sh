#!/bin/bash

if [ "$1" != "Release" ] && [ "$1" != "Debug" ]; then
  echo "Argument must be either 'release' or 'debug'"
  exit 1
fi

# Configure the project
cmake -DCMAKE_BUILD_TYPE=$1 -S . -B build

# Build the project
cmake --build build
