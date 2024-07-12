#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for solver_results..."
$CXXTESTGEN --error-printer -o test_runner_solver_results.cpp unit_tests/unittest_solver_results.t.h
$CPP -std=c++11 -o test_runner_solver_results test_runner_solver_results.cpp \
solvers/solver_results.cpp -I. -Isolvers -Ifunction_objects -Icore -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for solver_results:"
  ./test_runner_solver_results
  exit
else
  echo -e "${RED}Compilation of unit tests for solver_results.cpp failed.${NC}"
  exit 1
fi
