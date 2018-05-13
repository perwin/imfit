// Unit tests for PsfInterpolator classes
//

// See run_unittest_psf_interpolators.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>
#include <vector>
using namespace std;

#include "core/image_io.h"
#include "function_objects/psf_interpolators.h"


#define DELTA  1.0e-8

// These files are assumed to exist
string  psfImage_filename = string("unit_tests/psf_gauss_sigma0.5_5.fits");


//  PsfInterpolator( double *inputImage, int nCols_image, int nRows_image ) { ; };


class TestPsfInterpolator_bicubic : public CxxTest::TestSuite 
{
  // data members
  int  nColsPsf, nRowsPsf;
  double  *psfPixels;
  PsfInterpolator *psfInterp;
  
public:
  void setUp()
  {
    psfPixels = ReadImageAsVector(psfImage_filename, &nColsPsf, &nRowsPsf);
    psfInterp = new PsfInterpolator_bicubic(psfPixels, nColsPsf, nRowsPsf);
  }

  void tearDown()
  {
    delete psfInterp;
    free(psfPixels);
  }


  // and now the actual tests

  void testGetInterpolatorType( void )
  {
    int returnVal = psfInterp->GetInterpolatorType();
    TS_ASSERT_EQUALS( returnVal, kInterpolator_bicubic );
  }

  void testGetValues_noshift( void )
  {
    double returnVal0, returnVal1, returnVal2;
    
    // central pixel
    returnVal0 = psfInterp->GetValue(0.0,0.0);
    TS_ASSERT_DELTA( returnVal0, 0.73212016, DELTA );
    // 1 pixel to right of center
    returnVal1 = psfInterp->GetValue(1.0,0.0);
    TS_ASSERT_DELTA( returnVal1, 0.16868566, DELTA );
    // 1 pixel above center
    returnVal1 = psfInterp->GetValue(0.0,1.0);
    TS_ASSERT_DELTA( returnVal1, 0.16868566, DELTA );
    // 2 pixels below center
    returnVal2 = psfInterp->GetValue(0.0,-2.0);
    TS_ASSERT_DELTA( returnVal2, 0.0014417765, DELTA );
  }

  void testGetValues_shifted( void )
  {
    double returnVal0, returnVal1, returnVal2;
    
    // 0.5 pixels to right of central pixel
    returnVal0 = psfInterp->GetValue(0.5,0.0);
    TS_ASSERT_DELTA( returnVal0, 0.51973038, DELTA );
    // 1.5 pixels to right of center, 0.5 above
    returnVal1 = psfInterp->GetValue(1.5,0.5);
    TS_ASSERT_DELTA( returnVal1, 0.00880937, DELTA );
    // 1.5 pixels above center
    returnVal1 = psfInterp->GetValue(0.0,1.5);
    TS_ASSERT_DELTA( returnVal1, 0.01243073, DELTA );
  }
};


class TestPsfInterpolator_lanczos2 : public CxxTest::TestSuite 
{
  // data members
  int  nColsPsf, nRowsPsf;
  double  *psfPixels;
  PsfInterpolator *psfInterp;
  
public:
  void setUp()
  {
    psfPixels = ReadImageAsVector(psfImage_filename, &nColsPsf, &nRowsPsf);
    psfInterp = new PsfInterpolator_lanczos2(psfPixels, nColsPsf, nRowsPsf);
  }

  void tearDown()
  {
    delete psfInterp;
    free(psfPixels);
  }


  // and now the actual tests

  // Test to see if we can successfully add PSF vector to OversampledRegion object
  void testGetInterpolatorType( void )
  {
    int returnVal = psfInterp->GetInterpolatorType();
    TS_ASSERT_EQUALS( returnVal, kInterpolator_lanczos2 );
  }

  void testGetValues_noshift( void )
  {
    double returnVal0, returnVal1, returnVal2;
    
    // central pixel
    returnVal0 = psfInterp->GetValue(0.0,0.0);
    TS_ASSERT_DELTA( returnVal0, 0.73212016, DELTA );
    // 1 pixel to right of center
    returnVal1 = psfInterp->GetValue(1.0,0.0);
    TS_ASSERT_DELTA( returnVal1, 0.16868566, DELTA );
    // 1 pixel above center
    returnVal1 = psfInterp->GetValue(0.0,1.0);
    TS_ASSERT_DELTA( returnVal1, 0.16868566, DELTA );
    // 2 pixels below center
    returnVal2 = psfInterp->GetValue(0.0,-2.0);
    TS_ASSERT_DELTA( returnVal2, 0.0014417765, DELTA );
  }


// Python code for generating reference values
// In [0]: from astropy.io import fits
// In [1]: import lanczos_test
// In [2]: psfImage_filename = <path-to-imfit> + "unit_tests/psf_gauss_sigma0.5_5.fits"
// In [3]: psfim = fits.getdata(psfImage_filename)
// In [4]: ny,nx = psfim.shape
// In [5]: x_arr = np.arange(1,nx + 1)
// In [6]: y_arr = np.arange(1,ny + 1)
// # center of image is at [2,2] in numpy access, but at x=3,y=3 in our code
// # center of PSF image (no shift)
// In [7]: lanczos_test.LanczosInterp2D(x_arr,y_arr,psfim, 2, 3.0, 3.0)
// Out[8]: 0.732120156288147
// # interpolated value at center + 0.5 pixels in x
// In [9]: lanczos_test.LanczosInterp2D(x_arr,y_arr,psfim, 2, 3.5, 3.0)
// Out[10]: 0.5054706567431267
// # interpolated value at center + 1.5 pixels in x, + 0.5 pixels in y
// In [11]: lanczos_test.LanczosInterp2D(x_arr,y_arr,psfim, 2, 4.5, 3.5)
// Out[12]: 0.035119419221815946
// # interpolated value at center + 1.5 pixels in y
// In [13]: lanczos_test.LanczosInterp2D(x_arr,y_arr,psfim, 2, 3.0, 4.5)
// Out[14]: 0.050885502118477165
// # interpolated value at center - 1.5 pixels in x and y
// In [15]: lanczos_test.LanczosInterp2D(x_arr,y_arr,psfim, 2, 1.5, 1.5)
// Out[16]: 0.00352206909108757
// # interpolated value at center - 1.9 pixels in x and y
// In [92]: lanczos_test.LanczosInterp2D(x_arr,y_arr,psfim, 2, 1.1, 1.1)
// Out[92]: 0.0002065470426474115

  void testGetValues_shifted( void )
  {
    double returnVal0, returnVal1, returnVal2;
    
    // 0.5 pixels to right of central pixel
    returnVal0 = psfInterp->GetValue(0.5,0.0);
    TS_ASSERT_DELTA( returnVal0, 0.5054706567431267, DELTA );
    // 1.5 pixels to right of center, 0.5 above
    returnVal1 = psfInterp->GetValue(1.5,0.5);
    TS_ASSERT_DELTA( returnVal1, 0.035119419221815946, DELTA );
    // 1.5 pixels above center
    returnVal1 = psfInterp->GetValue(0.0,1.5);
    TS_ASSERT_DELTA( returnVal1, 0.050885502118477165, DELTA );
    // 1.5 pixels to left of center, 1.5 pixels below center
    returnVal1 = psfInterp->GetValue(-1.5,-1.5);
    TS_ASSERT_DELTA( returnVal1, 0.00352206909108757, DELTA );
    // 1.9 pixels to left of center, 1.9 pixels below center
    returnVal1 = psfInterp->GetValue(-1.9,-1.9);
    TS_ASSERT_DELTA( returnVal1, 0.0002065470426474115, DELTA );
  }

};



// Tests for auxiliary functions

class TestLanczosFunction : public CxxTest::TestSuite
{
  // data members

public:

  void testLanczosFunction_n2( void )
  {
    const int n = 2;
    const int nVals = 6;
    double returnVal1, returnVal2;
    double inputVals[nVals] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0};
    double correctVals[nVals] = {1.0, 0.5731591682507563, 0.0, -0.06368435202786181,
    						0.0, 0.0};
    
    // test that we get correct answers, including when x < 0 (should be same
    // as |x|)
    for (int i = 0; i < nVals; i++) {
      returnVal1 = Lanczos(inputVals[i], n);
      TS_ASSERT_DELTA( returnVal1, correctVals[i], DELTA );
      returnVal2 = Lanczos(-inputVals[i], n);
      TS_ASSERT_DELTA( returnVal2, correctVals[i], DELTA );
    }
  }

  void testLanczosFunction_n3( void )
  {
    const int n = 3;
    const int nVals = 8;
    double returnVal1, returnVal2;
    double inputVals[nVals] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0};
    double correctVals[nVals] = {1.0, 0.6079271018540265, 0.0, -0.13509491152311703,
    						0.0, 0.02431708407416106, 0.0, 0.0};
    
    // test that we get correct answers, including when x < 0 (should be same
    // as |x|)
    for (int i = 0; i < nVals; i++) {
      returnVal1 = Lanczos(inputVals[i], n);
      TS_ASSERT_DELTA( returnVal1, correctVals[i], DELTA );
      returnVal2 = Lanczos(-inputVals[i], n);
      TS_ASSERT_DELTA( returnVal2, correctVals[i], DELTA );
    }
  }
};


class testFindIndex : public CxxTest::TestSuite
{
  // data members

public:

  void testFindIndex_general( void )
  {
    double xArray[5] = {-2.0, -1.0, 0.0, 1.0, 2.0};
    int returnVal;
    double inputXvals[9] = {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
    int correctIndices[9] = {0, 0, 1, 1, 2, 2, 3, 3, 4};
    for (int i = 0; i < 9; i++) {
      returnVal = FindIndex(xArray, inputXvals[i]);
      TS_ASSERT_EQUALS(returnVal, correctIndices[i]);
    }
  }
};
