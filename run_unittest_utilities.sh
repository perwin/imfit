#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for utilities..."
$CXXTESTGEN --error-printer -o test_runner_utilities.cpp unit_tests/unittest_utilities.t.h 
$CPP -std=c++11 -o test_runner_utilities test_runner_utilities.cpp core/utilities.cpp \
-I. -Icore -Isolvers -Ifunction_objects -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for utilities:"
  ./test_runner_utilities
  exit
else
  echo -e "${RED}Compilation of unit tests for utilities.cpp failed.${NC}"
  exit 1
fi
