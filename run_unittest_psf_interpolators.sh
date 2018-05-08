#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for PsfInterpolator classes..."
$CXXTESTGEN --error-printer -o test_runner_funcs.cpp unit_tests/unittest_psf_interpolators.t.h 
$CPP -std=c++11 -o test_runner_funcs test_runner_funcs.cpp \
function_objects/psf_interpolators.cpp core/image_io.cpp \
-I/usr/local/include -I$CXXTEST -I. -Icore -Isolvers -Ifunction_objects \
-L/usr/local/lib -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
if [ $? -eq 0 ]
then
  echo "Running unit tests for PsfInterpolator classes:"
  ./test_runner_funcs
  exit
else
  echo -e "${RED}Compilation of unit tests for PsfInterpolator classes failed.${NC}"
  exit 1
fi
