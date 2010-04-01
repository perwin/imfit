#!/bin/bash
#
# Unit tests for imfit/makeimage:

# Unit tests for utilities
echo
echo "Generating and compiling unit tests for utilities..."
cxxtestgen.py --error-printer -o test_runner_utilities.cpp unittest_utilities.h 
g++ -o test_runner_utilities test_runner_utilities.cpp utilities.cpp -I/usr/local/include
echo "Running unit tests for utilities:"
./test_runner_utilities

# Unit tests for function objects
echo
echo "Generating and compiling unit tests for function objects..."
cxxtestgen.py --error-printer -o test_runner_funcs.cpp unittest_funcs.h 
g++ -o test_runner_funcs test_runner.cpp function_object.cpp func_exp.cpp \
       func_sersic.cpp -I/usr/local/include
echo "Running unit tests for function objects:"
./test_runner_funcs

echo
echo "Done with unit tests!"
echo
