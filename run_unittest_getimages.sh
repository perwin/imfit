#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for getimages..."
$CXXTESTGEN --error-printer -o test_runner_getimages.cpp unit_tests/unittest_getimages.t.h 
$CPP -std=c++11 -o test_runner_getimages test_runner_getimages.cpp \
core/getimages.cpp core/image_io.cpp -lcfitsio -lfftw3 \
-I. -I/usr/local/include -Icore -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for getimages:"
  ./test_runner_getimages TestGetMaskAndErrorImages
  exit
else
  echo "** Compilation of unit tests for getimages.cpp failed."
  exit 1
fi
