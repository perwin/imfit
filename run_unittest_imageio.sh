#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for image_io..."
$CXXTESTGEN --error-printer -o test_runner_imageio.cpp unit_tests/unittest_image_io.t.h 
$CPP -std=c++11 -o test_runner_imageio test_runner_imageio.cpp core/image_io.cpp -I. \
-I/usr/local/include -Icore -I$CXXTEST -lcfitsio -lfftw3
if [ $? -eq 0 ]
then
  echo "Running unit tests for image_io:"
  ./test_runner_imageio
  exit
else
  echo -e "${RED}Compilation of unit tests for add_functions.cpp failed.${NC}"
  exit 1
fi
