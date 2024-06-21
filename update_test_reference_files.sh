#!/bin/bash

# Simple script to copy current output files from running do_imfit_tests to the
# appropriate OS+CPU-specific subdirectory of tests/imfit_reference, for those tests
# where numerical output is (slightly) difference from one system to another.

# This should ideally only need to be run when using a new OS or CPU architecture.

# replace as appropriate ("macos_os86-64", "macos_arm64", "linux_x86-64")
ARCH_TEST_DIR="macos_arm64"

cp temptest/test_dump2 tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout2
cp temptest/test_dump3 tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout3
cp temptest/test_dump3d tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout3d
cp temptest/test_dump3d2 tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout3d2
cp temptest/test_dump3e tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout3e
cp temptest/test_dump3e2 tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout3e2
cp temptest/test_dump4 tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout4
cp temptest/test_dump4b tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout4b
cp temptest/test_dump4e3 tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout4e3
cp temptest/test_dump5d_tail tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout5d_tail
cp temptest/temp_bootstrap_output_tail tests/imfit_reference/${ARCH_TEST_DIR}/bootstrap_output_seed10_tail.dat
cp temptest/test_dump5e_tail tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout5e_tail
cp temptest/temp_bootstrap_output2_tail tests/imfit_reference/${ARCH_TEST_DIR}/bootstrap_output_seed10_2_tail.dat
cp temptest/test_dump8a tests/imfit_reference/${ARCH_TEST_DIR}/imfit_textout8a
