#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for config-file parser..."
$CXXTESTGEN --error-printer -o test_runner_config.cpp unit_tests/unittest_config_parser.t.h
$CPP -std=c++11 -o test_runner_config test_runner_config.cpp core/config_file_parser.cpp \
core/utilities.cpp -I. -Icore -Isolvers -Ifunction_objects -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for config_file_parser:"
  ./test_runner_config NewTestSuite
  exit
else
  echo -e "${RED}Compilation of unit tests for config_file_parser.cpp failed.${NC}"
  exit 1
fi
