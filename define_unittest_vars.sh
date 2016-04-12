#!/bin/bash

# This script defines the correct environment variables for running unit tests.
# Specifically, it defines the path to the CxxTest installation and the CxxTest
# binary (cxxtestgen); it also defines which compilers are used for compiling
# C and C++ code.

# optional specification if we're running on a Mac (e.g., specify specific version 
# of gcc/g++ to avoid getting LLVM/clang by default)
if [[ $OSTYPE == darwin* ]]
then
  CPP=g++
  CC=gcc
#   CPP=g++-5
#   CC=gcc-5
else
  CPP=g++
  CC=gcc
fi

# Set the path to CxxTest (and thus cxxtestgen)
# (Travis CI defines TRAVIS=true)
if env | grep -q ^TRAVIS=
then
  # OK, running on Travis CI, so it's the apt-get install location
  CXXTEST=/usr
else
  # Not Travis CI; use path to local CxxTest installation
  # (change this to the appropriate path if yours is different!)
  CXXTEST=/usr/local/cxxtest-4.4
fi
CXXTESTGEN=$CXXTEST/bin/cxxtestgen

export CPP
export CC
export CXXTEST
export CCTESTGEN
