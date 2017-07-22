#include <cxxtest/TestSuite.h>

#include <stdio.h>
#include <string>
#include <string>
#include <vector>
using namespace std;
#include "fftw3.h"
#include "getimages.h"
#include "options_base.h"
#include "definitions.h"


const string NONEXISTENT_IMAGENAME("bobobbo.fits");
const string NONEXISTENT_IMAGENAME2("bobobbo2.fits");
const string NONIMAGE_FILE("tests/test_table.fits");
const string ONES_IMAGE_3x3("tests/testimage_3x3_ones.fits");
const string IMAGE_3x3("tests/testimage_3x3_onezero.fits");
const int N_PIXELS_ONES_IMAGE = 9;


// std::tuple<double *, int> GetAndCheckImage( const string imageName, const string imageType,
// 											int nColumns_ref, int nRows_ref );
// std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
// 										OptionsBase *options,  bool &maskPixelsAllocated, 
// 										bool &errorPixelsAllocated );

const double  DELTA = 1.0e-9;
const double  DELTA_e7 = 1.0e-7;

class TestGetAndCheckImage : public CxxTest::TestSuite
{
public:

  void testBadName( void )
  {
    int  status;
    double  *pixels;
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    
    std::tie(pixels, status) = GetAndCheckImage(NONEXISTENT_IMAGENAME, "test", nCols_ref, nRows_ref);
    TS_ASSERT_EQUALS(status, -1);
    TS_ASSERT_EQUALS(pixels, nullptr);
  }

  void testBadFile( void )
  {
    int  status;
    double  *pixels;
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    
    std::tie(pixels, status) = GetAndCheckImage(NONIMAGE_FILE, "test", nCols_ref, nRows_ref);
    TS_ASSERT_EQUALS(status, -1);
    TS_ASSERT_EQUALS(pixels, nullptr);
  }

  void testSimpleImage_badSize( void )
  {
    int  status;
    double  *pixels;
    int  nCols_ref = 100;
    int  nRows_ref = 100;
    
    std::tie(pixels, status) = GetAndCheckImage(ONES_IMAGE_3x3, "test", nCols_ref, nRows_ref);
    TS_ASSERT_EQUALS(status, -2);
  }

  void testSimpleImage_good( void )
  {
    int  status;
    double  *pixels;
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    
    std::tie(pixels, status) = GetAndCheckImage(ONES_IMAGE_3x3, "test", nCols_ref, nRows_ref);
    TS_ASSERT_EQUALS(status, 0);
    for (int i = 0; i < N_PIXELS_ONES_IMAGE; i++)
      TS_ASSERT_EQUALS(pixels[i], 1.0);
    
    fftw_free(pixels);
  }

  // Don't test for matching dimensions if input dimensions = (0,0)
  void testSimpleImage_ignoreDimensions( void )
  {
    int  status;
    double  *pixels;
    int  nCols_ref = 0;
    int  nRows_ref = 0;
    
    std::tie(pixels, status) = GetAndCheckImage(ONES_IMAGE_3x3, "psf", nCols_ref, nRows_ref);
    TS_ASSERT_EQUALS(status, 0);
    for (int i = 0; i < N_PIXELS_ONES_IMAGE; i++)
      TS_ASSERT_EQUALS(pixels[i], 1.0);
    
    fftw_free(pixels);
  }
};


class TestGetMaskAndErrorImages : public CxxTest::TestSuite
{
public:

  void testNoImagesRequested( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, 0);
    
    delete options;
  }

  void testBadMaskName( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->maskImagePresent = true;
    options->maskFileName = NONEXISTENT_IMAGENAME;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, -1);
    
    delete options;
  }

  void testBadMaskFile( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->maskImagePresent = true;
    options->maskFileName = NONIMAGE_FILE;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, -1);
    
    delete options;
  }

  void testBadMaskSize( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 10;
    int  nRows_ref = 10;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->maskImagePresent = true;
    options->maskFileName = ONES_IMAGE_3x3;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, -1);
    TS_ASSERT_EQUALS(maskAllocated, false);
    TS_ASSERT_EQUALS(errAllocated, false);
    
    delete options;
  }

  void testGoodMaskFile( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->maskImagePresent = true;
    options->maskFileName = ONES_IMAGE_3x3;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, 1);
    TS_ASSERT_EQUALS(maskAllocated, true);
    TS_ASSERT_EQUALS(errAllocated, false);
    for (int i = 0; i < N_PIXELS_ONES_IMAGE; i++)
      TS_ASSERT_EQUALS(maskPixels[i], 1.0);
    TS_ASSERT_EQUALS(errPixels, nullptr);

    delete options;
    fftw_free(maskPixels);
  }

  void testBadErrSize( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 10;
    int  nRows_ref = 10;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->noiseImagePresent = true;
    options->noiseFileName = ONES_IMAGE_3x3;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, -1);
    TS_ASSERT_EQUALS(maskAllocated, false);
    TS_ASSERT_EQUALS(errAllocated, false);
    
    delete options;
  }

  void testGoodErrFile( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->noiseImagePresent = true;
    options->noiseFileName = ONES_IMAGE_3x3;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, 2);
    TS_ASSERT_EQUALS(maskAllocated, false);
    TS_ASSERT_EQUALS(errAllocated, true);
    for (int i = 0; i < N_PIXELS_ONES_IMAGE; i++)
      TS_ASSERT_EQUALS(errPixels[i], 1.0);
    TS_ASSERT_EQUALS(maskPixels, nullptr);

    delete options;
    fftw_free(errPixels);
  }

  void testGetBoth( void )
  {
    int  status;
    double  *maskPixels;
    double  *errPixels;
    OptionsBase *options = new OptionsBase();
    int  nCols_ref = 3;
    int  nRows_ref = 3;
    bool  maskAllocated = false;
    bool  errAllocated = false;

    options->maskImagePresent = true;
    options->noiseImagePresent = true;
    options->maskFileName = ONES_IMAGE_3x3;
    options->noiseFileName = ONES_IMAGE_3x3;
    std::tie(maskPixels, errPixels, status) = GetMaskAndErrorImages(nCols_ref, nRows_ref,
    				 							options, maskAllocated, errAllocated);
    TS_ASSERT_EQUALS(status, 3);
    TS_ASSERT_EQUALS(maskAllocated, true);
    TS_ASSERT_EQUALS(errAllocated, true);
    for (int i = 0; i < N_PIXELS_ONES_IMAGE; i++)
      TS_ASSERT_EQUALS(maskPixels[i], 1.0);
    for (int i = 0; i < N_PIXELS_ONES_IMAGE; i++)
      TS_ASSERT_EQUALS(errPixels[i], 1.0);

    delete options;
    fftw_free(maskPixels);
    fftw_free(errPixels);
  }
};
