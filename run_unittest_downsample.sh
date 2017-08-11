#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for downsample..."
$CXXTESTGEN --error-printer -o test_runner_downsample.cpp unit_tests/unittest_downsample.t.h 
$CPP -std=c++11 -o test_runner_downsample test_runner_downsample.cpp core/downsample.cpp \
core/image_io.cpp -I. -Icore -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for downsample:"
  ./test_runner_downsample
  exit
else
  echo -e "${RED}Compilation of unit tests for add_functions.cpp failed.${NC}"
  exit 1
fi
