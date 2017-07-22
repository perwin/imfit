#!/bin/bash
#
# Unit tests for imfit/makeimage, using individual shell scripts for each
# subset of unit tests.

# The following has been moved into the individual run_unittest_*.sh shell scripts

# determine if we're running on a Mac (specify specific version of
# gcc/g++ to avoid , or in a Travis CI VM
# (Travis CI defines TRAVIS=true)
# if [[ $OSTYPE == darwin* ]]
# then
#   CPP=g++
#   CC=gcc
#   CPP=g++-5
#   CC=gcc-5
# else
#   CPP=g++
#   CC=gcc
# fi
# 
# Set the path to cxxtestgen depending on whether we're running under Travis or not
# if env | grep -q ^TRAVIS=
# then
#   CXXTEST=/usr
# else
#   CXXTEST=/usr/local/cxxtest-4.4
# fi
# CXXTESTGEN=$CXXTEST/bin/cxxtestgen


# integer-variable counter to keep track of how many errors occur (if no errors, 
# then value will remain = 0)
declare -i RESULT=0


# Unit tests for add_functions
./run_unittest_add_functions.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for command-line parser
./run_unittest_cmlineparser.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for config-file parser
./run_unittest_configfileparser.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for downsample
./run_unittest_downsample.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for standard image-function classes
./run_unittest_funcs.sh 2>> ../temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for getimages
./run_unittest_getimages.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for image_io
./run_unittest_imageio.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for mpfit
./run_unittest_mpfit.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for solver_results
./run_unittest_solverresults.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for utilities
./run_unittest_utilities.sh 2>> temperror.log
RESULT+=$?
echo $RESULT

# Unit tests for model_object
./run_unittest_model_object.sh
RESULT+=$?
echo $RESULT

# Unit tests for options classes
./run_unittest_options.sh
RESULT+=$?
echo $RESULT

# Unit tests for setup_model_object
./run_unittest_setup_model_object.sh
RESULT+=$?
echo $RESULT

# Unit tests for psf_oversampling_info
./run_unittest_psf_oversampling_info.sh
RESULT+=$?
echo $RESULT

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
# function_objects/func_gaussian.cpp -I. -Icore -Isolvers -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
# echo "Running unit tests for oversampled_region:"
# ./test_runner_oversampled_region 2>> temperror.log



# Unit tests for model_object (for speed, we'll assume things have already been compiled...)
# echo
# echo "Generating and compiling unit tests for model_object..."
# $CXXTESTGEN --error-printer -o test_runner_modelobj.cpp unit_tests/unittest_model_object.t.h
# $CPP -fopenmp -o test_runner_modelobj test_runner_modelobj.cpp core/model_object.cpp \
# core/utilities.o core/convolver.o \
# core/add_functions.o core/config_file_parser.o c_code/mersenne_twister.o c_code/mp_enorm.o \
# core/oversampled_region.o core/downsample.o core/image_io.o \
# function_objects/function_object.o function_objects/func_gaussian.o \
# function_objects/func_exp.o function_objects/func_gen-exp.o \
# function_objects/func_sersic.o function_objects/func_gen-sersic.o \
# function_objects/func_core-sersic.o function_objects/func_broken-exp.o \
# function_objects/func_broken-exp2d.o function_objects/func_moffat.o \
# function_objects/func_flatsky.o function_objects/func_gaussian-ring.o \
# function_objects/func_gaussian-ring2side.o function_objects/func_edge-on-disk_n4762.o \
# function_objects/func_edge-on-disk_n4762v2.o function_objects/func_edge-on-ring.o \
# function_objects/func_edge-on-ring2side.o function_objects/func_edge-on-disk.o \
# function_objects/integrator.o function_objects/func_expdisk3d.o function_objects/func_brokenexpdisk3d.o \
# function_objects/func_gaussianring3d.o function_objects/func_king.o \
# function_objects/func_king2.o \
# -I. -Icore -Ic_code -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lcfitsio -lfftw3 -lgsl -lm -std=c++11
# RESULT+=$?
# echo "Running unit tests for model_object:"
# ./test_runner_modelobj 2>> temperror.log
# RESULT+=$?


echo
echo "Done with unit tests!"
echo

exit $RESULT

