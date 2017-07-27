#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for config-file parser..."
$CXXTESTGEN --error-printer -o test_runner_config.cpp unit_tests/unittest_config_parser.t.h
$CPP -std=c++11 -o test_runner_config test_runner_config.cpp core/config_file_parser.cpp core/utilities.cpp -I. -Icore -Isolvers -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for config_file_parser:"
  ./test_runner_config
  exit
else
  echo "Compilation of unit tests for config_file_parser.cpp failed."
  exit 1
fi
