image_io.h
==========

Overview
--------

These are utility functions for inspecting and reading in FITS files;
they are wrappers around calls to CFITSIO functions.

Note that ReadImageAsVector allocates the memory for the image it reads
(via the FFTW3 utility routine fftw_alloc_real); the caller is
responsible for freeing the memory (via fftw_free) when it is no longer
needed.


API
---

**Files:** ``core/image_io.h``, ``core/image_io.cpp``


.. doxygenfile:: image_io.h
