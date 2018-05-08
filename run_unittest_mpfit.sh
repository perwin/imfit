#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for mpfit..."
$CXXTESTGEN --error-printer -o test_runner_mpfit.cpp unit_tests/unittest_mpfit.t.h 
$CPP -std=c++11 -o test_runner_mpfit test_runner_mpfit.cpp core/mp_enorm.cpp \
core/utilities.cpp solvers/mpfit.cpp -I. -Icore -Isolvers  -Ifunction_objects \
-I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for mpfit:"
  ./test_runner_mpfit
  exit
else
  echo -e "${RED}Compilation of unit tests for mpfit.cpp failed.${NC}"
  exit 1
fi
