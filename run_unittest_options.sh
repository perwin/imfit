#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for options_makeimage..."
$CXXTESTGEN --error-printer -o test_runner_options.cpp unit_tests/unittest_options.t.h
$CPP -std=c++11 -o test_runner_options test_runner_options.cpp -I. -Icore -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for OptionsBase and derived classes:"
  ./test_runner_options
  exit
else
  echo -e "${RED}Compilation of unit tests for OptionsBase and derived classes failed.${NC}"
  exit 1
fi
