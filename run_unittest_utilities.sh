#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for utilities..."
$CXXTESTGEN --error-printer -o test_runner_utilities.cpp unit_tests/unittest_utilities.t.h 
$CPP -o test_runner_utilities test_runner_utilities.cpp core/utilities.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for utilities:"
  ./test_runner_utilities
  exit
else
  echo "Compilation of unit tests for utilities.cpp failed."
  exit 1
fi
