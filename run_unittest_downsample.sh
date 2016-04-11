#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for downsample..."
$CXXTESTGEN --error-printer -o test_runner_downsample.cpp unit_tests/unittest_downsample.t.h 
$CPP -o test_runner_downsample test_runner_downsample.cpp core/downsample.cpp core/image_io.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for downsample:"
  ./test_runner_downsample
  exit
else
  echo "Compilation of unit tests for downsample.cpp failed."
  exit 1
fi
