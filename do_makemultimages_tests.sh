#!/bin/bash
#
# Regression tests for makemultimages

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


if [ ! -f ./makemultimages ]
  then
    echo "*** Unable to find local copy of makemultimages! ***"
    echo
    exit 2
fi

case "$OSTYPE" in
  darwin*)  osname="osx" ;; 
  linux*)   osname="linux" ;;
  *)        echo "unknown: $OSTYPE" ;;
esac


# create output directory for test_dump*, etc. files (will do nothing if
# directory already exists)
mkdir -p temptest


echo ""
echo "Running tests for makemultimages..."
# do we correctly print list of available functions?
./makemultimages --list-functions > temptest/test_dump0

# simple tests using small images
rm temptest/multimage_diff-sizes*.fits
rm temptest/multimage-with-imageinfo_*.fits

# 3 images, one function block, using --image-info
./makemultimages tests_for_makemultimage/config_for_globalparams.txt -o temptest/multimage-with-imageinfo --image-info=tests_for_makemultimage/imageinfo_3image.txt > temptest/test_dump2

./makemultimages tests_for_makemultimage/config_for_globalparams.txt -o temptest/multimage-with-imageinfo_b --image-info=tests_for_makemultimage/imageinfo_3image_b.txt > temptest/test_dump2b

./makemultimages tests_for_makemultimage/config_for_globalparams.txt -o temptest/multimage-with-imageinfo_c --image-info=tests_for_makemultimage/imageinfo_3image_c.txt > temptest/test_dump2c

./makemultimages tests_for_makemultimage/config_for_globalparams.txt -o temptest/multimage-with-imageinfo_d --image-info=tests_for_makemultimage/imageinfo_3image_d.txt > temptest/test_dump2d

./makemultimages tests_for_makemultimage/config_for_globalparams.txt -o temptest/multimage-with-imageinfo_e --image-info=tests_for_makemultimage/imageinfo_3image_e.txt > temptest/test_dump2d


# model with two function blocks, 3 images (offsets only), using --image-info
./makemultimages tests_for_makemultimage/config_for_globalparams_2blocks.txt -o temptest/multimage-with-imageinfo_2blocks --image-info=tests_for_makemultimage/imageinfo_3image.txt > temptest/test_dump3


# model with 2 images, different sizes, using --image-info
./makemultimages tests_for_makemultimage/config_for_globalparams_mult-sizes.txt -o temptest/multimage_diff-sizes --image-info=tests_for_makemultimage/imageinfo_2sizes.txt > temptest/test_dump5

# model with 3 images, different sizes, 1st and 3rd with PSF, using --image-info
./makemultimages tests_for_makemultimage/config_for_globalparams_mult-sizes.txt -o temptest/multimage_diff-sizes-psf --image-info=tests_for_makemultimage/imageinfo_3sizes-with-psf.txt > temptest/test_dump6

# tests of PSF oversampling [NO COMPARISON WITH REFERENCE YET]
./makemultimages tests/makeimage_reference/config_makeimage_gauss-oversample.dat --image-info tests_for_makemultimage/imageinfo_2sizes-with-oversampled-psf.txt -o temptest/multimage_osamp  > temptest/test_dump7a

# tests of multi-region PSF oversampling [NO COMPARISON WITH REFERENCE YET]
./makemultimages tests/makeimage_reference/config_makeimage_gauss-oversample.dat --image-info tests_for_makemultimage/imageinfo_2sizes-with-2oversampled-psf.txt -o temptest/multimage_2osamp  > temptest/test_dump7b




echo ""

echo -n "*** Diff comparison with archives: printing function list... "
if (diff --brief temptest/test_dump0 tests/makemultimages_textout0)
then
  echo " OK"
else
  echo "Diff output:"
  diff temptest/test_dump0 tests/makemultimages_textout0
  STATUS+=1
fi


echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_0... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_0.fits tests/multimagetest_0.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_1... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_1.fits tests/multimagetest_1.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_2.fits tests/multimagetest_2.fits
  STATUS+=$?
fi


echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_b_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_b_2.fits tests/multimagetest_2b.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_c_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_c_2.fits tests/multimagetest_2c.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_d_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_d_2.fits tests/multimagetest_2d.fits
  STATUS+=$?
fi

# note that the following two use comparisons with reference images from previous tests
echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_e_1... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_e_1.fits tests/multimagetest_1.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_e_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_e_2.fits tests/multimagetest_2d.fits
  STATUS+=$?
fi


echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_2blocks_0... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_2blocks_0.fits tests/multimagetest_2blocks_0.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_2blocks_1... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_2blocks_1.fits tests/multimagetest_2blocks_offset_1.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage-with-imageinfo_2blocks_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage-with-imageinfo_2blocks_2.fits tests/multimagetest_2blocks_offset_2.fits
  STATUS+=$?
fi


# archive images generated by
# makeimage[v1.4] -c tests_for_makemultimage/config_for_globalparams_2sizes.txt
echo -n "*** Diff comparison of image with archives: multimage_diff-sizes_0... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage_diff-sizes_0.fits tests/multimagetest_diff-sizes_0.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage_diff-sizes_1... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage_diff-sizes_1.fits tests/multimagetest_diff-sizes_1.fits
  STATUS+=$?
fi


# archive images generated by
# makeimage[v1.4] -c tests_for_makemultimage/config_for_globalparams_2sizes.txt --psf=tests/psf_moffat_35_n4699z.fits
echo -n "*** Diff comparison of image with archives: multimage_diff-sizes-psf_0... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage_diff-sizes-psf_0.fits tests/multimagetest_diff-sizes-psf_0.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage_diff-sizes-psf_1... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage_diff-sizes-psf_1.fits tests/multimagetest_diff-sizes-psf_1.fits
  STATUS+=$?
fi

echo -n "*** Diff comparison of image with archives: multimage_diff-sizes-psf_2... "
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py temptest/multimage_diff-sizes-psf_2.fits tests/multimagetest_diff-sizes-psf_2.fits
  STATUS+=$?
fi


# WARNING: We currently don't have a reference image for this yet!
# Note that we chose an OS-specific path for the reference file (image output is
# marginally different for Linux vs OS X)
# archive image generated by
# makimage[v1.3.1] tests/config_makeimage_gauss-oversample.dat --psf=tests/psf_standard.fits
echo -n "*** Diff comparison of PSF-oversampled image with archives: multimage_diff-sizes-psf_2... "
echo -n "[WARNING: NO REFERENCE IMAGE YET!]"
if [ $do_fits_tests = "1" ]
then
  echo ""
  ./python/compare_fits_files.py --min-value=1.0e-7 temptest/multimage_osamp_1.fits tests/${osname}/oversampled_orig.fits
  STATUS+=$?
fi




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

