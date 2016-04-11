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
echo "Generating and compiling unit tests for function objects..."
$CXXTESTGEN --error-printer -o test_runner_funcs.cpp unit_tests/unittest_funcs.t.h 
$CPP -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
function_objects/func_exp.cpp function_objects/func_flatsky.cpp \
function_objects/func_gaussian.cpp function_objects/func_moffat.cpp \
function_objects/func_sersic.cpp function_objects/func_king.cpp function_objects/func_king2.cpp \
function_objects/func_edge-on-disk.cpp \
-I/usr/local/include -I$CXXTEST -I. -Icore -Ic_code -Isolvers -lm -lgsl -lgslcblas
if [ $? -eq 0 ]
then
  echo "Running unit tests for function objects:"
  ./test_runner_funcs
  exit
else
  echo "Compilation of unit tests for function objects failed."
  exit 1
fi
