#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for getimages..."
$CXXTESTGEN --error-printer -o test_runner_getimages.cpp unit_tests/unittest_getimages.t.h 
$CPP -std=c++11 -o test_runner_getimages test_runner_getimages.cpp \
core/getimages.cpp core/image_io.cpp core/psf_oversampling_info.cpp core/utilities.cpp \
-L/usr/local/lib -lcfitsio -lfftw3 -I. -I/usr/local/include -Icore -Isolvers \
-Ifunction_objects -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for getimages:"
  ./test_runner_getimages
  exit
else
  echo -e "${RED}Compilation of unit tests for getimages.cpp failed.${NC}"
  exit 1
fi
