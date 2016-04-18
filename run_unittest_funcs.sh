#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for function objects..."
$CXXTESTGEN --error-printer -o test_runner_funcs.cpp unit_tests/unittest_funcs.t.h 
$CPP -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
function_objects/func_exp.cpp function_objects/func_flatsky.cpp \
function_objects/func_gaussian.cpp function_objects/func_moffat.cpp \
function_objects/func_sersic.cpp function_objects/func_king.cpp function_objects/func_king2.cpp \
function_objects/func_edge-on-disk.cpp \
-I/usr/local/include -I$CXXTEST -I. -Icore -Isolvers -lm -lgsl -lgslcblas
if [ $? -eq 0 ]
then
  echo "Running unit tests for function objects:"
  ./test_runner_funcs
  exit
else
  echo "Compilation of unit tests for function objects failed."
  exit 1
fi
