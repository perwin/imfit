#!/bin/bash
#
# Unit tests for imfit/makeimage:

# define things explicitly (rather than relying on PATH), so we can call this from
# SConstruct
CXXTEST=/usr/local/cxxtest-4.4
CXXTESTGEN=$CXXTEST/bin/cxxtestgen
CPP=g++-4.9


# Unit tests for utilities
echo
echo "Generating and compiling unit tests for utilities..."
$CXXTESTGEN --error-printer -o test_runner_utilities.cpp unit_tests/unittest_utilities.t.h 
$CPP -o test_runner_utilities test_runner_utilities.cpp utilities.cpp -I. -Isolvers -I/usr/local/include -I$CXXTEST
echo "Running unit tests for utilities:"
./test_runner_utilities 2> temperror_unittests.log

# Unit tests for command-line parser
echo
echo "Generating and compiling unit tests for command-line parser..."
$CXXTESTGEN --error-printer -o test_runner_cmparser.cpp unit_tests/unittest_commandline_parser.t.h 
$CPP -Wno-write-strings -o test_runner_cmparser test_runner_cmparser.cpp commandline_parser.cpp utilities.cpp -I. -Isolvers -I/usr/local/include -I$CXXTEST
echo "Running unit tests for command-line parser:"
./test_runner_cmparser 2>> temperror.log

# Unit tests for config-file parser
echo
echo "Generating and compiling unit tests for config-file parser..."
$CXXTESTGEN --error-printer -o test_runner_config.cpp unit_tests/unittest_config_parser.t.h
$CPP -o test_runner_config test_runner_config.cpp config_file_parser.cpp utilities.cpp -I. -Isolvers -I/usr/local/include -I$CXXTEST
echo "Running unit tests for config-file parser:"
./test_runner_config 2>> temperror.log

# Unit tests for oversampled_region
echo
echo "Generating and compiling unit tests for oversampled_region..."
$CXXTESTGEN --error-printer -o test_runner_oversampled_region.cpp unit_tests/unittest_downsample.t.h 
$CPP -o test_runner_oversampled_region test_runner_oversampled_region.cpp downsample.cpp \
convolver.cpp image_io.cpp function_objects/function_object.cpp \
function_objects/func_gaussian.cpp -I. -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
echo "Running unit tests for oversampled_region:"
./test_runner_oversampled_region 2>> temperror.log

# Unit tests for downsample
echo
echo "Generating and compiling unit tests for downsample..."
$CXXTESTGEN --error-printer -o test_runner_downsample.cpp unit_tests/unittest_downsample.t.h 
$CPP -o test_runner_downsample test_runner_downsample.cpp downsample.cpp image_io.cpp -I. -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lm
echo "Running unit tests for downsample:"
./test_runner_downsample 2>> temperror.log

echo
echo "Generating and compiling unit tests for mpfit..."
$CXXTESTGEN --error-printer -o test_runner_mpfit.cpp unit_tests/unittest_mpfit.t.h 
$CPP -o test_runner_mpfit test_runner_mpfit.cpp solvers/mpfit.cpp -I. -Isolvers -Isolvers -I/usr/local/include -I$CXXTEST
echo "Running unit tests for utilities:"
./test_runner_mpfit 2>> temperror.log


# Unit tests for function objects
echo
echo "Generating and compiling unit tests for function objects..."
$CXXTESTGEN --error-printer -o test_runner_funcs.cpp unit_tests/unittest_funcs.t.h 
$CPP -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
function_objects/func_exp.cpp function_objects/func_gaussian.cpp function_objects/func_sersic.cpp \
-I/usr/local/include -I$CXXTEST -I. -Isolvers -lm -lgsl
# $CPP -o test_runner_funcs test_runner_funcs.cpp function_objects/function_object.cpp \
# function_objects/func_exp.o function_objects/func_sersic.o \
# function_objects/func_gaussian.o function_objects/func_edge-on-disk.o \
# -I. -Isolvers -I/usr/local/include -I$CXXTEST -lm -lgsl
echo "Running unit tests for function objects:"
./test_runner_funcs 2>> ../temperror.log

# Unit tests for add_functions
echo
echo "Generating and compiling unit tests for add_functions..."
$CXXTESTGEN --error-printer -o test_runner_add_functions.cpp unit_tests/unittest_add_functions.t.h
$CPP -fopenmp -o test_runner_add_functions test_runner_add_functions.cpp add_functions.cpp \
model_object.o utilities.cpp convolver.o \
config_file_parser.o mersenne_twister.o mp_enorm.o \
oversampled_region.o downsample.o image_io.o \
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
-I. -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lfftw3 -lcfitsio -lgsl -lm
echo "Running unit tests for add_functions:"
./test_runner_add_functions 2>> temperror.log

# Unit tests for model_object
echo
echo "Generating and compiling unit tests for model_object..."
$CXXTESTGEN --error-printer -o test_runner_modelobj.cpp unit_tests/unittest_model_object.t.h
$CPP -fopenmp -o test_runner_modelobj test_runner_modelobj.cpp model_object.cpp utilities.cpp convolver.o \
add_functions.o config_file_parser.o mersenne_twister.o mp_enorm.o \
oversampled_region.o downsample.o image_io.o \
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
-I. -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lcfitsio -lfftw3 -lgsl -lm -std=c++11
echo "Running unit tests for model_object:"
./test_runner_modelobj 2>> temperror.log


echo
echo "Done with unit tests!"
echo
