#!/bin/bash

# determine if we're running locally, or in a Travis CI VM:
if env | grep -q ^TRAVIS=
then
  CXXTEST=/usr
  CPP=g++
  CC=gcc
else
  CXXTEST=/usr/local/cxxtest-4.4
  CPP=g++-5
  CC=gcc-5
fi
#CXXTEST=/usr/local/cxxtest-4.4
CXXTESTGEN=$CXXTEST/bin/cxxtestgen
#CPP=g++-5
#CC=gcc-5

echo
echo "Generating and compiling unit tests for image_io..."
$CXXTESTGEN --error-printer -o test_runner_imageio.cpp unit_tests/unittest_image_io.t.h 
$CPP -o test_runner_imageio test_runner_imageio.cpp core/image_io.cpp -I. -I/usr/local/include -Icore -I$CXXTEST -lcfitsio -lfftw3
if [ $? -eq 0 ]
then
  echo "Running unit tests for utilities:"
  ./test_runner_imageio
  exit
else
  echo "Compilation of unit tests for image_io.cpp failed."
  exit 1
fi
