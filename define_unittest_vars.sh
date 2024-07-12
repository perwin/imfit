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


# Locations of external libraries (CFITSIO, FFTW3, GSL, NLopt)
# Figure out which type of macOS architecture we're running under

# defaults for Linux and for x86-64 macOS
EXTERNAL_INCLUDE_PATH="/usr/local/include"
EXTERNAL_LIB_PATH="/usr/local/lib"
ARCH=$(uname -m)
# if we're running on an arm64 (Apple silicon) Mac, use Homebrew's arm64 location
if [[ $OSTYPE == darwin* ]] && [[ "$ARCH" -eq "arm64" ]]
then
  EXTERNAL_INCLUDE_PATH="/opt/homebrew/include"
  EXTERNAL_LIB_PATH="/opt/homebrew/lib"
fi


export CPP
export CC
export CXXTEST
export CXXTESTGEN
export LDFLAGS
export EXTERNAL_INCLUDE_PATH
export EXTERNAL_LIB_PATH
