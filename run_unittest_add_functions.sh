#!/bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color


echo
echo "Generating and compiling unit tests for add_functions..."
$CXXTESTGEN --error-printer -o test_runner_add_functions.cpp unit_tests/unittest_add_functions.t.h
$CPP -std=c++11 -o test_runner_add_functions test_runner_add_functions.cpp core/add_functions.cpp \
core/model_object.cpp core/utilities.cpp core/convolver.cpp core/config_file_parser.cpp  \
core/oversampled_region.cpp core/downsample.cpp core/image_io.cpp core/psf_oversampling_info.cpp \
function_objects/function_object.cpp function_objects/func_gaussian.cpp \
function_objects/func_exp.cpp function_objects/func_gen-exp.cpp \
function_objects/func_sersic.cpp function_objects/func_gen-sersic.cpp \
function_objects/func_core-sersic.cpp function_objects/func_broken-exp.cpp \
function_objects/func_broken-exp2d.cpp function_objects/func_moffat.cpp \
function_objects/func_flatsky.cpp function_objects/func_tilted-sky-plane.cpp \
function_objects/func_flatbar.cpp \
function_objects/func_gaussian-ring.cpp function_objects/func_gaussian-ring2side.cpp \
function_objects/func_gaussian-ring-az.cpp function_objects/func_edge-on-disk_n4762.cpp \
function_objects/func_edge-on-disk_n4762v2.cpp function_objects/func_edge-on-ring.cpp \
function_objects/func_edge-on-ring2side.cpp function_objects/func_edge-on-disk.cpp \
function_objects/integrator.cpp function_objects/func_expdisk3d.cpp function_objects/func_brokenexpdisk3d.cpp \
function_objects/func_gaussianring3d.cpp function_objects/func_ferrersbar3d.cpp \
function_objects/func_ferrersbar2d.cpp function_objects/func_king.cpp \
function_objects/func_king2.cpp function_objects/func_pointsource.cpp \
function_objects/func_pointsource-rot.cpp \
function_objects/helper_funcs.cpp function_objects/helper_funcs_3d.cpp \
function_objects/psf_interpolators.cpp \
core/mersenne_twister.cpp core/mp_enorm.cpp \
-I. -Icore -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -L/usr/local/lib \
-lfftw3_threads -lfftw3 -lcfitsio -lgsl -lgslcblas -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for add_functions:"
  ./test_runner_add_functions
  exit
else
  echo -e "${RED}Compilation of unit tests for add_functions.cpp failed.${NC}"
  exit 1
fi
