#! /bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

# Unit tests for model_object (includes extra -lgslcblas in case we compile on Linux)
echo
echo "Generating and compiling unit tests for model_object..."
$CXXTESTGEN --error-printer -o test_runner_modelobj.cpp unit_tests/unittest_model_object.t.h
$CPP -std=c++11 -fsanitize=address -Wl,-no_compact_unwind -DDEBUG -DUSE_TEST_FUNCS \
-o test_runner_modelobj \
test_runner_modelobj.cpp core/model_object.cpp core/utilities.cpp core/convolver.cpp \
core/add_functions.cpp core/config_file_parser.cpp core/mersenne_twister.cpp \
core/mp_enorm.cpp core/oversampled_region.cpp core/downsample.cpp \
core/image_io.cpp core/psf_oversampling_info.cpp \
function_objects/function_object.cpp function_objects/func_gaussian.cpp \
function_objects/func_exp.cpp function_objects/func_gen-exp.cpp \
function_objects/func_sersic.cpp function_objects/func_gen-sersic.cpp \
function_objects/func_core-sersic.cpp function_objects/func_broken-exp.cpp \
function_objects/func_broken-exp2d.cpp function_objects/func_moffat.cpp \
function_objects/func_flatsky.cpp function_objects/func_gaussian-ring.cpp \
function_objects/func_gaussian-ring2side.cpp function_objects/func_edge-on-ring.cpp \
function_objects/func_edge-on-ring2side.cpp function_objects/func_edge-on-disk.cpp \
function_objects/integrator.cpp function_objects/func_expdisk3d.cpp \
function_objects/func_brokenexpdisk3d.cpp function_objects/func_gaussianring3d.cpp \
function_objects/func_ferrersbar3d.cpp function_objects/func_king.cpp \
function_objects/func_king2.cpp function_objects/func_gauss_extraparams.cpp \
function_objects/func_pointsource.cpp \
function_objects/helper_funcs.cpp function_objects/helper_funcs_3d.cpp \
function_objects/psf_interpolators.cpp \
-I. -Icore -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST \
-L/usr/local/lib -lfftw3_threads -lcfitsio -lfftw3 -lgsl -lgslcblas -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for model_object:"
  ./test_runner_modelobj
  exit
else
  echo -e "${RED}Compilation of unit tests for model_object.cpp failed.${NC}"
  exit 1
fi
