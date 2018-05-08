#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for command-line parser..."
$CXXTESTGEN --error-printer -o test_runner_cmparser.cpp unit_tests/unittest_commandline_parser.t.h 
$CPP -std=c++11 -Wno-write-strings -o test_runner_cmparser test_runner_cmparser.cpp \
core/commandline_parser.cpp core/utilities.cpp -I. -Icore -Isolvers -Ifunction_objects \
-I/usr/local/include -I$CXXTEST -L/usr/local/lib 
if [ $? -eq 0 ]
then
  echo "Running unit tests for commandline_parser:"
  ./test_runner_cmparser
  exit
else
  echo -e "${RED}Compilation of unit tests for commandline_parser.cpp failed.${NC}"
  exit 1
fi
