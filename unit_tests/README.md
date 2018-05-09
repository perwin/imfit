## Unit Tests

This directory contains code for unit tests for (some of) Imfit's various functions and classes.

The shell script `run_unit_tests.sh` in the top-level directory will compile and execute all
the tests; it does so by calling the individual `run_unittest_*` shell scripts in the
top-level directory. (You can also execute the individual scripts by themselves if you
just want to run one set of tests.)

To run the tests, you will need to install the
[CxxTest](http://cxxtest.com) unit-test framework. You should also
modify the shell script `define_unittest_vars.sh` (in the top-level
directory), which specifies the paths for CxxTest and also which
compilers are used. (Unlike the case for compiling imfit, imfit-mcmc, or makeimage,
the unit tests can be compiled with LLVM/clang as well as with GCC, since
they don't use OpenMP.)

Note that the "unit tests" for oversampled_region.cpp (unittest\_oversampled\_region.t.h) are not
proper *tests* -- they merely generate output FITS files intended to be inspected in, e.g.,
SAOimage DS9. Accordingly, they are *not* executed by the `run_unit_tests.sh`
script.
