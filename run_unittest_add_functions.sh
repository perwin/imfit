#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

echo
echo "Generating and compiling unit tests for add_functions..."
$CXXTESTGEN --error-printer -o test_runner_add_functions.cpp unit_tests/unittest_add_functions.t.h
$CPP -o test_runner_add_functions test_runner_add_functions.cpp core/add_functions.cpp \
core/model_object.cpp core/utilities.cpp core/convolver.cpp core/config_file_parser.cpp  \
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
core/mersenne_twister.cpp core/mp_enorm.cpp \
-I. -Icore -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lfftw3 -lcfitsio -lgsl -lgslcblas -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for add_functions:"
  ./test_runner_add_functions
  exit
else
  echo "Compilation of unit tests for add_functions.cpp failed."
  exit 1
fi
