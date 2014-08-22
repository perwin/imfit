// Unit tests for downsample.cpp
//
// Includes notes (in comments) on how references images were made in Python
//
// cxxtestgen --error-printer -o test_runner.cpp unittest_downsample.h
// g++ -o test_runner test_runner.cpp downsample.cpp image_io.cpp -I/usr/local/include -I$CXXTEST -lcfitsio -lm
// ./test_runner

#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>

using namespace std;

#include "image_io.h"
#include "downsample.h"

#define DELTA  1.0e-10

// * 20x20 image of zeros
// >>> da = np.zeros((20,20))
// >>> pyfits.writeto(fimfit+"simpleimage_20x20_zeros.fits",da)
string  simpleNullImage_filename = string("simpleimage_20x20_zeros.fits");

// * 6x6 image of ones
// >>> db = np.ones((6,6))
// >>> pyfits.writeto(fimfit+"simpleimage_6x6_ones.fits",db)
string  osampOnesImage_filename = string("simpleimage_6x6_ones.fits");
// * 6x6 image of ones, but with LL corner pixel = 10
// >>> db[0:1,0:1] = 1.0
// >>> pyfits.writeto(fimfit+"simpleimage_6x6_ones+ten.fits",db)
string  osampOnesImage3_filename = string("simpleimage_6x6_ones+ten.fits");
// * 10x10 image with inner 6x6 composed of ones, outer 2 rows & columns on each side = 0
// e.g., 6x6 all-ones image with padding for 2x2 PSF (padding set = 0 to help catch
// possible downsampling errors -- should be ignored by downsampling)
// >>> dc = np.zeros((10,10))
// >>> dc[2:8,2:8] = 1.0
// >>> pyfits.writeto(fimfit+"simpleimage_10x10_inner6x6ones.fits",dc)
string  osampOnesPlusZeroBorderImage_filename = string("simpleimage_10x10_inner6x6ones.fits");
// * Same, but now with LL corner of central image (interior to padding) = 1.0
// >>> dc[2:3,2:3] = 10.0
// >>> pyfits.writeto(fimfit+"simpleimage_10x10_inner6x6ones+llten.fits",dc)
string  osampOnesPlusZeroBorderImage2_filename = string("simpleimage_10x10_inner6x6ones+llten.fits");
// * Same, but now with UR corner of central image (interior to padding) = 1.0
// >>> dc[2:8,2:8] = 1.0
// >>> dc[7:8,7:8] = 10.0
// >>> pyfits.writeto(fimfit+"simpleimage_10x10_inner6x6ones+urten.fits",dc)
string  osampOnesPlusZeroBorderImage3_filename = string("simpleimage_10x10_inner6x6ones+urten.fits");

// Desired output images (no PSF padding)
// * desired output of downsample-and-replace with overSample = 1 -- has [5:11,11:16] = 1.0
string  modifiedNullImage1_filename = string("simpleimage_20x20_zeros+ones6x6.fits");
// * desired output of downsample-and-replace with overSample = 3 -- has [5:6,11:13] = 1.0
string  modifiedNullImage2_filename = string("simpleimage_20x20_zeros+ones2x2.fits");
// * desired output of downsample-and-replace with overSample = 3, using osampOnesImage3_filename -- has [5:6,11:13] = 1.0,
// except [5,11] = 2.0
string  modifiedNullImage3_filename = string("simpleimage_20x20_zeros+3x3mix.fits");

// * Desired output images (PSF padding) -- 
// desired output of downsample-and-replace with overSample = 3, assumes main-image with 2x2 PSF padding -- has [7:8,13:14] = 1.0 (full-image coords)
string  modifiedNullImage_with_psf1_filename = string("simpleimage_20x20-with-psf_zeros+3x3ones.fits");



class TestDownsample : public CxxTest::TestSuite 
{
  // data members
  int  nColsMain, nRowsMain;
  int  nColsOsamp, nRowsOsamp;
  int  debug0;

  
public:
  void setUp()
  {
    debug0 = 0;
  }
  

  // and now the actual tests

  void test1x1Downsample_noPSF( void )
  {
    double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
    double *osampImage = ReadImageAsVector(osampOnesImage_filename, &nColsOsamp, &nRowsOsamp);
    double *refFinalMainImage1 = ReadImageAsVector(modifiedNullImage1_filename, &nColsMain, &nRowsMain);

    int correct_nColsMain, correct_nRowsMain, correct_nColsOsamp, correct_nRowsOsamp;
    correct_nColsMain = correct_nRowsMain = 20;
    correct_nColsOsamp = correct_nRowsOsamp = 6;

    // downsample and replace, using oversampleSCale = 1 and start coord = [5,11]; no PSF padding
    DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,0,0, mainImage, nColsMain,nRowsMain,0,0,
     					5,11, 1, debug0);

// We'd like to use the following, but it actually doesn't work: same-sized images
// with different data are claimed to be the same.
//    TS_ASSERT_SAME_DATA(mainImage, refImageMain1, nColsMain*nRowsMain);
// So instead we'll do it the hard way:
    for (int k = 0; k < nColsMain*nRowsMain; k++) {
      TS_ASSERT_DELTA(mainImage[k], refFinalMainImage1[k], DELTA);
    }
  }

  void test3x3Downsample_noPSF( void )
  {
    double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
    double *osampImage = ReadImageAsVector(osampOnesImage_filename, &nColsOsamp, &nRowsOsamp);
    double *refFinalMainImage2 = ReadImageAsVector(modifiedNullImage2_filename, &nColsMain, &nRowsMain);

    int correct_nColsMain, correct_nRowsMain, correct_nColsOsamp, correct_nRowsOsamp;
    correct_nColsMain = correct_nRowsMain = 20;
    correct_nColsOsamp = correct_nRowsOsamp = 6;

    // downsample and replace, using oversampleSCale = 3 and start coord = [5,11]; no PSF padding
    DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,0,0, mainImage, nColsMain,nRowsMain,0,0,
     					5,11, 3, debug0);

    for (int k = 0; k < nColsMain*nRowsMain; k++) {
      TS_ASSERT_DELTA(mainImage[k], refFinalMainImage2[k], DELTA);
    }
  }

  void test3x3Downsample2_noPSF( void )
  {
    double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
    double *osampImage = ReadImageAsVector(osampOnesImage3_filename, &nColsOsamp, &nRowsOsamp);
    double *refFinalMainImage = ReadImageAsVector(modifiedNullImage3_filename, &nColsMain, &nRowsMain);

    // downsample and replace, using oversampleSCale = 3 and start coord = [5,11]; no PSF padding
    DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,0,0, mainImage, nColsMain,nRowsMain,0,0,
     					5,11, 3, debug0);

    for (int k = 0; k < nColsMain*nRowsMain; k++) {
      TS_ASSERT_DELTA(mainImage[k], refFinalMainImage[k], DELTA);
    }
  }


  // case of PSF padding in main image, but not in oversampled image
  // take base 20x20 null image and assume it's a 16x16 image with 2x2 PSF
  void test3x3Downsample_mainPSFOnly( void )
  {
    double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
    double *osampImage = ReadImageAsVector(osampOnesImage_filename, &nColsOsamp, &nRowsOsamp);
    double *refFinalMainImage = ReadImageAsVector(modifiedNullImage_with_psf1_filename, &nColsMain, &nRowsMain);

    // downsample and replace, using oversampleSCale = 3 and start coord = [5,11]; assume
    // main image (only) has PSF padding for 2x2 PSF
    DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,0,0, mainImage, nColsMain,nRowsMain,2,2,
     					5,11, 3, debug0);
    vector<string>  dummy;
    int status = SaveVectorAsImage(mainImage, string("bob.fits"), nColsMain, nRowsMain, dummy);

    for (int k = 0; k < nColsMain*nRowsMain; k++) {
      TS_ASSERT_DELTA(mainImage[k], refFinalMainImage[k], DELTA);
    }
  }


  // case of PSF padding in both main & oversampled images
  // take base 20x20 null image and assume it's a 16x16 image with 2x2 PSF
  void test3x3Downsample_bothPSF1( void )
  {
    double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
    double *osampImage = ReadImageAsVector(osampOnesPlusZeroBorderImage_filename, &nColsOsamp, &nRowsOsamp);
    double *refFinalMainImage = ReadImageAsVector(modifiedNullImage_with_psf1_filename, &nColsMain, &nRowsMain);

    // downsample and replace, using oversampleSCale = 3 and start coord = [5,11]; assume
    // main and oversampled images both have PSF padding for 2x2 PSF
    DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,2,2, mainImage, nColsMain,nRowsMain,2,2,
     					5,11, 3, debug0);
    vector<string>  dummy;
    int status = SaveVectorAsImage(mainImage, string("bob.fits"), nColsMain, nRowsMain, dummy);

    for (int k = 0; k < nColsMain*nRowsMain; k++) {
      TS_ASSERT_DELTA(mainImage[k], refFinalMainImage[k], DELTA);
    }
  }

  void test3x3Downsample_bothPSF2( void )
  {
    double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
    double *osampImage = ReadImageAsVector(osampOnesPlusZeroBorderImage2_filename, &nColsOsamp, &nRowsOsamp);
    double *refFinalMainImage = ReadImageAsVector(modifiedNullImage3_filename, &nColsMain, &nRowsMain);

    // downsample and replace, using oversampleSCale = 3 and start coord = [3,9] (so it should
    // end up the same as the non-PSF case with start coord = [5,11], above); assume
    // main and oversampled images both have PSF padding for 2x2 PSF
    DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,2,2, mainImage, nColsMain,nRowsMain,2,2,
     					3,9, 3, debug0);
    vector<string>  dummy;
    int status = SaveVectorAsImage(mainImage, string("bob.fits"), nColsMain, nRowsMain, dummy);

    for (int k = 0; k < nColsMain*nRowsMain; k++) {
      TS_ASSERT_DELTA(mainImage[k], refFinalMainImage[k], DELTA);
    }
  }
  
};
