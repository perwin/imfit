#!/bin/bash
#
# Regression tests for multimfit

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

# create an integer variable for status (& counting number of failed tests)
declare -i STATUS=0

FAILED_TEXT=" ** Test failed. **"

do_fits_tests=$(./python/py_startup_test.py)
if [ $do_fits_tests = "0" ]
then
  echo -n "Unable to load numpy and/or pyfits Python modules"
  echo
  echo -n "Tests will skip comparison of FITS files"
fi


if [ ! -f ./multimfit ]
  then
    echo "*** Unable to find local copy of multimfit! ***"
    echo
    exit 2
fi

# Determine OS and processor type, so we can compare against the right references files 
# (this accounts for the fact that different processor-OS combinations produce slightly 
# different numerical output)
OS_NAME=$(uname -s)
ARCH=$(uname -m)

if [[ $OS_NAME == "Darwin" ]]
then
    OS_NAME_BETTER="macos"
else
    OS_NAME_BETTER="linux"
fi
ARCH_TEST_DIR="${OS_NAME_BETTER}_${ARCH}"
echo "Using OS+architecture = ${ARCH_TEST_DIR}"


# create output directory for test_dump*, etc. files (will do nothing if
# directory already exists)
mkdir -p temptest


echo ""
echo "Running tests for multimfit..."
# do we correctly print list of available functions?
./multimfit --list-functions > temptest/test_dump0

# 2 tiny (2x2) images, one function block (FlatSky)
#    --fitstat-only test (correct value: CHI-SQUARE = 6.75)
./multimfit -c tests/multimfit_reference/config_imfit_smalldataimages.dat --image-info=tests/multimfit_reference/imageinfo_multimfit0.txt --fitstat-only > temptest/test_dump1a
#    minimization using L-M
./multimfit -c tests/multimfit_reference/config_imfit_smalldataimages.dat --image-info=tests/multimfit_reference/imageinfo_multimfit0.txt > temptest/test_dump1b
#    minimization using N-M
./multimfit -c tests/multimfit_reference/config_imfit_smalldataimages.dat --image-info=tests/multimfit_reference/imageinfo_multimfit0.txt --nm > temptest/test_dump1c
#    L-M minimization with config-file parameter limits
./multimfit -c tests/multimfit_reference/config_imfit_smalldataimages_with-limits.dat --image-info=tests/multimfit_reference/imageinfo_multimfit0.txt > temptest/test_dump1d
#    L-M minimization with image-info-file parameter limits
./multimfit -c tests/multimfit_reference/config_imfit_smalldataimages.dat --image-info=tests/multimfit_reference/imageinfo_multimfit0_with-limits.txt > temptest/test_dump1e

# SLOW: 2 small (200x200) images, one function block (Gaussian), with multiple 
# PSF oversampling regions + PMLR (use NM, since L-M doesn't really converge)
if [ "$1" == "--all" ]
then
  echo -n "   (now running fits with oversampled PSF ...)"
  ./multimfit -c tests/imfit_reference/config_imfit_gauss-oversample-test2.dat --image-info tests/multimfit_reference/imageinfo_multimfit+osample.txt --mlr --nm  > temptest/test_dump2a
  echo ""
else
  echo "   (-- skipping oversampled PSF tests --)"
fi

# 5x5 images, one function block (Gaussian), identical noise image used for both images
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3.txt --nm > temptest/test_dump3a
# same, but now with image subsections for 1st (reference) image
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3b.txt --nm > temptest/test_dump3b
# same, but now with image subsection for 2nd image
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3c.txt --nm > temptest/test_dump3c

# Same as previous, but w/o noise image
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3_noerror.txt --nm > temptest/test_dump3d
# Now using model errors
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3_noerror.txt --nm --model-errors > temptest/test_dump3e
# Now using Poisson MLR statistic (works OK with L-M minimization)
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3_noerror.txt --mlr > temptest/test_dump3f

# test bootstrap resampling
./multimfit -c tests/multimfit_reference/config_imfit_smallgauss.dat -i tests/multimfit_reference/imageinfo_multimfit3d.txt --mlr --seed 10 --bootstrap 5 --save-bootstrap temptest/temp_multimfit_bootstrap_output.dat > temptest/test_dump3g





echo ""

echo -n "*** Diff comparison with archives: printing function list... "
if (diff --brief temptest/test_dump0 tests/multimfit_reference/multimfit_textout0)
then
  echo " OK"
else
  echo "Diff output:"
  diff temptest/test_dump0 tests/multimfit_reference/multimfit_textout0
  STATUS+=1
fi


echo -n "*** Diff comparison with archives: simplest 2-image run (--fitstat-only)... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout1a temptest/test_dump1a --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: simplest 2-image run (L-M)... "
./python/diff_printouts.py tests/multimfit_reference/${ARCH_TEST_DIR}/multimfit_textout1b temptest/test_dump1b --skip-last=3
STATUS+=$?


echo -n "*** Diff comparison with archives: simplest 2-image run (N-M simplex)... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout1c temptest/test_dump1c --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: simplest 2-image run (config-file limits)... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout1d temptest/test_dump1d --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: simplest 2-image run (image-info-file limits)... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout1e temptest/test_dump1e --skip-last=3
STATUS+=$?


if [ "$1" == "--all" ]
then
  echo -n "*** Diff comparison with archives: fit using oversampled PSF... "
  ./python/diff_printouts.py tests/multimfit_reference/multimfit_textout2a temptest/test_dump2a --skip-last=3
  STATUS+=$?
fi


echo -n "*** Diff comparison with archives: 5x5 2-image run w/ sigma images... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout3a temptest/test_dump3a --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: same as previous + image 1 subsection... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout3b temptest/test_dump3b --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: same as previous + image 2 subsection... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout3c temptest/test_dump3c --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: same as previous, no sigma images... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout3d temptest/test_dump3d --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: same as previous, model errors... "
./python/diff_printouts.py tests/multimfit_reference/multimfit_textout3e temptest/test_dump3e --skip-last=3
STATUS+=$?

echo -n "*** Diff comparison with archives: same as previous, Poisson MLR errors... "
./python/diff_printouts.py tests/multimfit_reference/${ARCH_TEST_DIR}/multimfit_textout3f temptest/test_dump3f --skip-last=3
STATUS+=$?


echo -n "*** Diff comparison with archives: bootstrap resampling... "
./python/diff_printouts.py tests/multimfit_reference/${ARCH_TEST_DIR}/multimfit_textout3g temptest/test_dump3g --skip-last=3
STATUS+=$?




echo ""
if [ $STATUS -eq 1 ]
then
  echo -e "${RED}One test failed!${NC}"
elif [ $STATUS -gt 1 ]
then
  echo -e "${RED}${STATUS} tests failed!${NC}"
else
  echo -e "${GREEN}All tests passed.${NC}"
fi

echo ""
echo "Done."
echo ""

exit $STATUS
