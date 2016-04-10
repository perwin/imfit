#!/bin/bash
#
# Unit tests for imfit/makeimage:

# determine if we're running locally, or in a Travis CI VM
# (Travis CI defines TRAVIS=true)
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
CXXTESTGEN=$CXXTEST/bin/cxxtestgen


# counter to keep track of how many errors occur
declare -i RESULT=0


# Unit tests for utilities
echo
echo "Generating and compiling unit tests for utilities..."
$CXXTESTGEN --error-printer -o test_runner_utilities.cpp unit_tests/unittest_utilities.t.h 
$CPP -o test_runner_utilities test_runner_utilities.cpp core/utilities.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST
RESULT+=$?
echo "Running unit tests for utilities:"
./test_runner_utilities 2> temperror_unittests.log
RESULT+=$?

# Unit tests for command-line parser
echo
echo "Generating and compiling unit tests for command-line parser..."
$CXXTESTGEN --error-printer -o test_runner_cmparser.cpp unit_tests/unittest_commandline_parser.t.h 
$CPP -Wno-write-strings -o test_runner_cmparser test_runner_cmparser.cpp core/commandline_parser.cpp core/utilities.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST
RESULT+=$?
echo "Running unit tests for command-line parser:"
./test_runner_cmparser 2>> temperror.log
RESULT+=$?

# Unit tests for config-file parser
echo
echo "Generating and compiling unit tests for config-file parser..."
$CXXTESTGEN --error-printer -o test_runner_config.cpp unit_tests/unittest_config_parser.t.h
$CPP -o test_runner_config test_runner_config.cpp core/config_file_parser.cpp core/utilities.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST
RESULT+=$?
echo "Running unit tests for config-file parser:"
./test_runner_config 2>> temperror.log
RESULT+=$?

# NOTE: the following code will correctly set up and run unittest_oversampled_region.t.h;
# However, that "unit test" doesn't perform proper tests; instead, it generates test image
# output meant to be inspected with DS9. Therefore, it is currently commented out.
#
# Unit tests for oversampled_region
# echo
# echo "Generating and compiling unit tests for oversampled_region..."
# $CXXTESTGEN --error-printer -o test_runner_oversampled_region.cpp unit_tests/unittest_oversampled_region.t.h 
# $CPP -o test_runner_oversampled_region test_runner_oversampled_region.cpp core/downsample.cpp \
# core/oversampled_region.cpp core/convolver.cpp core/image_io.cpp function_objects/function_object.cpp \
# function_objects/func_gaussian.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
# echo "Running unit tests for oversampled_region:"
# ./test_runner_oversampled_region 2>> temperror.log

# Unit tests for downsample
echo
echo "Generating and compiling unit tests for downsample..."
$CXXTESTGEN --error-printer -o test_runner_downsample.cpp unit_tests/unittest_downsample.t.h 
$CPP -o test_runner_downsample test_runner_downsample.cpp core/downsample.cpp core/image_io.cpp -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
RESULT+=$?
echo "Running unit tests for downsample:"
./test_runner_downsample 2>> temperror.log
RESULT+=$?

echo
echo "Generating and compiling unit tests for mpfit..."
$CXXTESTGEN --error-printer -o test_runner_mpfit.cpp unit_tests/unittest_mpfit.t.h 
$CPP -o test_runner_mpfit test_runner_mpfit.cpp solvers/mpfit.cpp -I. -Icore -Ic_code -Isolvers -Isolvers -I/usr/local/include -I$CXXTEST
RESULT+=$?
echo "Running unit tests for utilities:"
./test_runner_mpfit 2>> temperror.log
RESULT+=$?


# Unit tests for function objects
echo
echo "Generating and compiling unit tests for function objects..."
$CXXTESTGEN --error-printer -o test_runner_funcs.cpp unit_tests/unittest_funcs.t.h 
$CPP -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
function_objects/func_exp.cpp function_objects/func_flatsky.cpp \
function_objects/func_gaussian.cpp function_objects/func_moffat.cpp \
function_objects/func_sersic.cpp function_objects/func_king.cpp function_objects/func_king2.cpp \
function_objects/func_edge-on-disk.cpp \
-I/usr/local/include -I$CXXTEST -I. -Icore -Ic_code -Isolvers -lm -lgsl
# $CPP -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
# function_objects/func_exp.o function_objects/func_sersic.o \
# function_objects/func_gaussian.o function_objects/func_edge-on-disk.o \
# -I. -Icore -Ic_code -Isolvers -I/usr/local/include -I$CXXTEST -lm -lgsl
RESULT+=$?
echo "Running unit tests for function objects:"
./test_runner_funcs 2>> ../temperror.log
RESULT+=$?

# Unit tests for add_functions
echo
echo "Generating and compiling unit tests for add_functions..."
$CXXTESTGEN --error-printer -o test_runner_add_functions.cpp unit_tests/unittest_add_functions.t.h
$CC -c c_code/mersenne_twister.c c_code/mp_enorm.c
RESULT+=$?
$CPP -fopenmp -o test_runner_add_functions test_runner_add_functions.cpp core/add_functions.cpp \
core/model_object.cpp core/utilities.cpp core/convolver.cpp \
core/config_file_parser.cpp c_code/mersenne_twister.o c_code/mp_enorm.o \
core/oversampled_region.cpp core/downsample.cpp core/image_io.cpp \
function_objects/function_object.cpp function_objects/func_gaussian.cpp \
function_objects/func_exp.cpp function_objects/func_gen-exp.cpp \
function_objects/func_sersic.cpp function_objects/func_gen-sersic.cpp \
function_objects/func_core-sersic.cpp function_objects/func_broken-exp.cpp \
function_objects/func_broken-exp2d.cpp function_objects/func_moffat.cpp \
function_objects/func_flatsky.cpp function_objects/func_gaussian-ring.cpp \
function_objects/func_gaussian-ring2side.cpp function_objects/func_edge-on-disk_n4762.cpp \
function_objects/func_edge-on-disk_n4762v2.cpp function_objects/func_edge-on-ring.cpp \
function_objects/func_edge-on-ring2side.cpp function_objects/func_edge-on-disk.cpp \
function_objects/integrator.cpp function_objects/func_expdisk3d.cpp function_objects/func_brokenexpdisk3d.cpp \
function_objects/func_gaussianring3d.cpp function_objects/func_king.cpp \
function_objects/func_king2.cpp \
-I. -Icore -Ic_code -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lfftw3 -lcfitsio -lgsl -lm
RESULT+=$?
echo "Running unit tests for add_functions:"
./test_runner_add_functions 2>> temperror.log
RESULT+=$?

# Unit tests for model_object (for speed, we'll assume things have already been compiled...)
echo
echo "Generating and compiling unit tests for model_object..."
$CXXTESTGEN --error-printer -o test_runner_modelobj.cpp unit_tests/unittest_model_object.t.h
$CPP -fopenmp -o test_runner_modelobj test_runner_modelobj.cpp core/model_object.cpp \
core/utilities.o core/convolver.o \
core/add_functions.o core/config_file_parser.o mersenne_twister.o mp_enorm.o \
core/oversampled_region.o core/downsample.o core/image_io.o \
function_objects/function_object.o function_objects/func_gaussian.o \
function_objects/func_exp.o function_objects/func_gen-exp.o \
function_objects/func_sersic.o function_objects/func_gen-sersic.o \
function_objects/func_core-sersic.o function_objects/func_broken-exp.o \
function_objects/func_broken-exp2d.o function_objects/func_moffat.o \
function_objects/func_flatsky.o function_objects/func_gaussian-ring.o \
function_objects/func_gaussian-ring2side.o function_objects/func_edge-on-disk_n4762.o \
function_objects/func_edge-on-disk_n4762v2.o function_objects/func_edge-on-ring.o \
function_objects/func_edge-on-ring2side.o function_objects/func_edge-on-disk.o \
function_objects/integrator.o function_objects/func_expdisk3d.o function_objects/func_brokenexpdisk3d.o \
function_objects/func_gaussianring3d.o function_objects/func_king.o \
function_objects/func_king2.o \
-I. -Icore -Ic_code -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lcfitsio -lfftw3 -lgsl -lm -std=c++11
RESULT+=$?
echo "Running unit tests for model_object:"
./test_runner_modelobj 2>> temperror.log
RESULT+=$?


echo
echo "Done with unit tests!"
echo

exit $RESULT

