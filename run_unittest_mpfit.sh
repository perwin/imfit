#!/bin/bash

# determine if we're running on a Mac (specify specific version of
# gcc/g++ to avoid , or in a Travis CI VM
# (Travis CI defines TRAVIS=true)
if [[ $OSTYPE == darwin* ]]
then
  CPP=g++
  CC=gcc
#   CPP=g++-5
#   CC=gcc-5
else
  CPP=g++
  CC=gcc
fi

# Set the path to cxxtestgen depending on whether we're running under Travis or not
if env | grep -q ^TRAVIS=
then
  CXXTEST=/usr
else
  CXXTEST=/usr/local/cxxtest-4.4
fi
CXXTESTGEN=$CXXTEST/bin/cxxtestgen

echo
echo "Generating and compiling unit tests for mpfit..."
$CXXTESTGEN --error-printer -o test_runner_mpfit.cpp unit_tests/unittest_mpfit.t.h 
$CPP -o test_runner_mpfit test_runner_mpfit.cpp solvers/mpfit.cpp -I. -Icore -Ic_code -Isolvers -Isolvers -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for mpfit:"
  ./test_runner_mpfit
  exit
else
  echo "Compilation of unit tests for mpfit.cpp failed."
  exit 1
fi
