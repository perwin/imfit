#! /bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Unit tests for model_object (includes extra -lgslcblas in case we compile on Linux)
echo
echo "Generating and compiling unit tests for model_object..."
$CXXTESTGEN --error-printer -o test_runner_modelobj.cpp unit_tests/unittest_model_object.t.h
$CPP -DDEBUG -DUSE_TEST_FUNCS -o test_runner_modelobj test_runner_modelobj.cpp core/model_object.cpp \
core/utilities.cpp core/convolver.cpp core/add_functions.cpp core/config_file_parser.cpp \
core/mersenne_twister.cpp core/mp_enorm.cpp core/oversampled_region.cpp core/downsample.cpp \
core/image_io.cpp \
function_objects/function_object.cpp function_objects/func_gaussian.cpp \
function_objects/func_exp.cpp function_objects/func_gen-exp.cpp \
function_objects/func_sersic.cpp function_objects/func_gen-sersic.cpp \
function_objects/func_core-sersic.cpp function_objects/func_broken-exp.cpp \
function_objects/func_broken-exp2d.cpp function_objects/func_moffat.cpp \
function_objects/func_flatsky.cpp function_objects/func_gaussian-ring.cpp \
function_objects/func_gaussian-ring2side.cpp function_objects/func_edge-on-ring.cpp \
function_objects/func_edge-on-ring2side.cpp function_objects/func_edge-on-disk.cpp \
function_objects/integrator.cpp function_objects/func_expdisk3d.cpp function_objects/func_brokenexpdisk3d.cpp \
function_objects/func_gaussianring3d.cpp function_objects/func_king.cpp \
function_objects/func_king2.cpp function_objects/func_gauss_extraparams.cpp \
-I. -Icore -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST \
-lfftw3_threads -lcfitsio -lfftw3 -lgsl -lgslcblas -lm

echo
echo "Running unit tests for model_object:"
./test_runner_modelobj
