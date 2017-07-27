#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for command-line parser..."
$CXXTESTGEN --error-printer -o test_runner_cmparser.cpp unit_tests/unittest_commandline_parser.t.h 
$CPP -std=c++11 -Wno-write-strings -o test_runner_cmparser test_runner_cmparser.cpp \
core/commandline_parser.cpp core/utilities.cpp -I. -Icore -Isolvers \
-I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for commandline_parser:"
  ./test_runner_cmparser
  exit
else
  echo "Compilation of unit tests for commandline_parser.cpp failed."
  exit 1
fi
