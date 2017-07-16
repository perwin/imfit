#! /bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Unit tests for setup_model_object
echo
echo "Generating and compiling unit tests for setup_model_object..."
$CXXTESTGEN --error-printer -o test_runner_setup_modelobj.cpp unit_tests/unittest_setup_model_object.t.h
$CPP -DDEBUG -fsanitize=address -o test_runner_setup_modelobj test_runner_setup_modelobj.cpp core/model_object.cpp \
core/setup_model_object.cpp core/utilities.cpp core/convolver.cpp core/config_file_parser.cpp \
core/mersenne_twister.cpp core/mp_enorm.cpp core/oversampled_region.cpp core/downsample.cpp \
core/image_io.cpp core/psf_oversampling_info.cpp \
-I. -Icore -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST \
-lfftw3_threads -lcfitsio -lfftw3 -lgsl -lgslcblas -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for setup_model_object:"
  ./test_runner_setup_modelobj NewTestSuite
  exit
else
  echo "** Compilation of unit tests for setup_model_object failed."
  exit 1
fi
