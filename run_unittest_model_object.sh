#! /bin/bash
cxxtestgen --error-printer -o test_runner.cpp unittest_model_object.h
g++-4.8 -o test_runner test_runner.cpp model_object.cpp utilities.cpp convolver.o \
add_functions.o config_file_parser.o mersenne_twister.o mp_enorm.o \
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
-I/usr/local/include -Ifunction_objects -I$CXXTEST -lfftw3_threads -lfftw3 -lgsl -lm && ./test_runner
