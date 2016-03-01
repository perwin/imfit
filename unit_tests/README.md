## unit_tests Directory

This directory contains code for unit tests for (some of) Imfit's various functions and classes.
(The shell script run_unit_tests.sh in the top-level directory will compile and execute all
the unit tests.)

To run these tests, you will need to install the [CxxTest](http://cxxtest.com) unit-test framework.

Note that the "unit tests" for oversampled_region (unittest_oversampled_region.t.h) are not
proper *tests* -- they merely generate output FITS files intended to be inspected in, e.g.,
SAOimage DS9.
