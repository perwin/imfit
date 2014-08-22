#!/bin/bash
#
# Unit tests for imfit/makeimage:

# Unit tests for utilities
echo
echo "Generating and compiling unit tests for utilities..."
cxxtestgen --error-printer -o test_runner_utilities.cpp unittest_utilities.h 
g++ -o test_runner_utilities test_runner_utilities.cpp utilities.cpp -I/usr/local/include -I$CXXTEST
echo
echo "Running unit tests for utilities:"
./test_runner_utilities

# Unit tests for command-line parser
echo
echo "Generating and compiling unit tests for command-line parser..."
cxxtestgen --error-printer -o test_runner_cmparser.cpp unittest_commandline_parser.h 
g++ -Wno-write-strings -o test_runner_cmparser test_runner_cmparser.cpp commandline_parser.cpp utilities.cpp -I/usr/local/include -I$CXXTEST
echo
echo "Running unit tests for command-line parser:"
./test_runner_cmparser

# Unit tests for config-file parser
echo
echo "Generating and compiling unit tests for config-file parser..."
cxxtestgen --error-printer -o test_runner_config.cpp unittest_config_parser.h
g++ -o test_runner_config test_runner_config.cpp config_file_parser.cpp utilities.cpp -I/usr/local/include -I$CXXTEST
echo
echo "Running unit tests for config-file parser:"
./test_runner_config

# Unit tests for function objects
cd function_objects
echo
echo "Generating and compiling unit tests for function objects..."
cxxtestgen --error-printer -o test_runner_funcs.cpp unittest_funcs.h 
g++ -o test_runner_funcs test_runner_funcs.cpp function_object.cpp func_exp.cpp \
       func_sersic.cpp func_gaussian.cpp func_edge-on-disk.cpp -I/usr/local/include -I$CXXTEST \
       -lm -lgsl
echo
echo "Running unit tests for function objects:"
./test_runner_funcs
cd ..

echo
echo "Done with unit tests!"
echo
