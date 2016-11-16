#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for options_makeimage..."
$CXXTESTGEN --error-printer -o test_runner_options.cpp unit_tests/unittest_options.t.h
$CPP -std=c++11 -o test_runner_options test_runner_options.cpp -I. -Icore -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for options_makeimage:"
  ./test_runner_options
  exit
else
  echo "Compilation of unit tests for options_makeimage.h, etc. failed."
  exit 1
fi
