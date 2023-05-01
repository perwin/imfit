#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

echo
echo "Generating and compiling unit tests for function objects..."
$CXXTESTGEN --error-printer -o test_runner_funcs.cpp unit_tests/unittest_funcs.t.h 
$CPP -std=c++11 -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
function_objects/func_exp.cpp function_objects/func_flatsky.cpp \
function_objects/func_gaussian.cpp function_objects/func_moffat.cpp \
function_objects/func_sersic.cpp function_objects/func_king.cpp function_objects/func_king2.cpp \
function_objects/func_nuker.cpp \
function_objects/func_broken-exp.cpp function_objects/func_double-broken-exp.cpp \
function_objects/func_broken-exp2d.cpp function_objects/func_edge-on-disk.cpp \
function_objects/func_gauss_extraparams.cpp function_objects/func_ferrersbar3d.cpp \
function_objects/func_pointsource.cpp function_objects/psf_interpolators.cpp \
function_objects_1d/func1d_exp_test.cpp \
function_objects/helper_funcs.cpp function_objects/helper_funcs_3d.cpp \
function_objects/integrator.cpp core/utilities.cpp \
-I/usr/local/include -I$CXXTEST -I. -Icore -Isolvers -Ifunction_objects \
-L/usr/local/lib -lm -lgsl -lgslcblas
if [ $? -eq 0 ]
then
  echo "Running unit tests for function objects:"
  ./test_runner_funcs
  exit
else
  echo -e "${RED}Compilation of unit tests for FunctionObject classes failed.${NC}"
  exit 1
fi
