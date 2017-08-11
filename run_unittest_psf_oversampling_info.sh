#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for psf_oversampling_info..."
$CXXTESTGEN --error-printer -o test_runner_psf_oversampling_info.cpp \
unit_tests/unittest_psf_oversampling_info.t.h
$CPP -std=c++11 -o test_runner_psf_oversampling_info \
test_runner_psf_oversampling_info.cpp core/psf_oversampling_info.cpp core/utilities.cpp \
-lfftw3 -I. -Icore -Isolvers -I/usr/local/include -I$CXXTEST
if [ $? -eq 0 ]
then
  echo "Running unit tests for psf_oversampling_info:"
  ./test_runner_psf_oversampling_info
  exit
else
  echo -e "${RED}Compilation of unit tests for add_functions.cpp failed.${NC}"
  exit 1
fi
