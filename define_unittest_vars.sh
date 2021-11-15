#!/bin/bash

# This script defines the correct environment variables for running unit tests.
# Specifically, it defines the path to the CxxTest installation and the CxxTest
# binary (cxxtestgen); it also defines which compilers are used for compiling
# C and C++ code.

# optional specification if we're running on a Mac (currently, we allow use
# of Clang/Clang++ because it provides nice error messages).
# "-w" suppresses all warnings
# "-Wl,-no_compact_unwind" is for suppressing the silly "no compact unwind"
# warnings from the clang linker.
if [[ $OSTYPE == darwin* ]]
then
  CPP="g++ -w -Wl,-no_compact_unwind"
  CC=gcc
#   CPP=g++-5
#   CC=gcc-5
else
  CPP="g++ -w"
  CC=gcc
fi

# attempt to suppress annoying, pointless "compact unwind" warnings
LDFLAGS="-Wl,-no_compact_unwind"

# Set the path to CxxTest
if [[ $OSTYPE == darwin* ]]
then
  CXXTEST=/usr/local
else
  # assume we're running on some version of Ubuntu
  CXXTEST=/usr
# simplest thing is to assume PATH includes location of cxxtestgen
fi
CXXTESTGEN=cxxtestgen


export CPP
export CC
export CXXTEST
export CXXTESTGEN
export LDFLAGS
