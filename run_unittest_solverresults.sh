#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for solver_results..."
$CXXTESTGEN --error-printer -o test_runner_solver_results.cpp unit_tests/unittest_solver_results.t.h
$CPP -std=c++11 -o test_runner_solver_results test_runner_solver_results.cpp solvers/solver_results.cpp -I. -Isolvers -Icore -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for solver_results:"
  ./test_runner_solver_results
  exit
else
  echo "Compilation of unit tests for solver_results.cpp failed."
  exit 1
fi
