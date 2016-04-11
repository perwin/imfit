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
echo "Generating and compiling unit tests for downsample..."
$CXXTESTGEN --error-printer -o test_runner_downsample.cpp unit_tests/unittest_downsample.t.h 
$CPP -o test_runner_downsample test_runner_downsample.cpp core/downsample.cpp core/image_io.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for downsample:"
  ./test_runner_downsample
  exit
else
  echo "Compilation of unit tests for downsample.cpp failed."
  exit 1
fi
