## unit_tests Directory

This directory contains code for unit tests for (some of) Imfit's various functions and classes.
(The shell script run\_unit\_tests.sh in the top-level directory will compile and execute all
the unit tests.)

To run these tests, you will need to install the [CxxTest](http://cxxtest.com) unit-test framework.

Note that the "unit tests" for oversampled_region.cpp (unittest\_oversampled\_region.t.h) are not
proper *tests* -- they merely generate output FITS files intended to be inspected in, e.g.,
SAOimage DS9. Accordingly, they are *not* automatically executed by the run\_unit\_tests.sh
script.
