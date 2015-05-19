#! /bin/bash

# define things explicitly (rather than relying on PATH), so we can call this from
# SConstruct
CXXTEST=/usr/local/cxxtest-4.4
CXXTESTGEN=$CXXTEST/bin/cxxtestgen
CPP=g++-4.9

# Unit tests for model_object
echo
echo "Generating and compiling unit tests for model_object..."
$CXXTESTGEN --error-printer -o test_runner_modelobj.cpp unittest_model_object.t.h
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
function_objects/func_gaussianring3d.o \
-I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lcfitsio -lfftw3 -lgsl -lm

echo
echo "Running unit tests for model_object:"
./test_runner_modelobj
