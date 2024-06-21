## tests Directory

This directory contains files for use by the regression tests for imfit,
imfit-mcmc, and makeimage (`do_makeimage_tests`, `do_imfit_tests`, and
`do_mcmc_tests` in the top-level directory). This includes input
configuration files, FITS images, and reference output files.

(See `run_unit_tests.sh` and the various `run_unittest_*.sh` scripts in
the parent directory for how to run unit tests; the actual unit-test header
files are found in the `../unit_test` directory.)

**Note:** The reference text files and FITS images are mostly tuned for
development work on a Mac. Minor differences in numerical output will
occur on a Linux system, and on x86-64 vs arm64 versions of macOS, due
to the vagaries of how floating-point numbers are handled. A partial
attempt to deal with this involves separate images in the `osx/` and
`linux/` subdirectories, and the different os+architecture
subdirectories inside `imfit_reference/`.
