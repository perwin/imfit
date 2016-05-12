// Unit tests for code in image_io.cpp

// See run_unittest_imageio.sh for how to compile and run these tests.

// older compilation notes:
// $CXXTESTGEN --error-printer -o test_runner_imageio.cpp unit_tests/unittest_image_io.t.h 
// $CPP -o test_runner_imageio test_runner_imageio.cpp core/image_io.cpp -I. -I/usr/local/include -Icore -I$CXXTEST -lcfitsio -lfftw3

// makeimage
// 	1. GetImageSize(options.referenceImageName, &nColumns, &nRows)
// 	2. ReadImageAsVector(options.psfFileName, &nColumns_psf, &nRows_psf)
// 	3. same as 2, for oversampled PSF
// 
// imfit
// 	Numerous calls to ReadImageAsVector
// // [ ] CheckImage
// [ ] GetImageSize with various HDUs



#include <cxxtest/TestSuite.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

//#include "utilities_pub.h"
#include "image_io.h"

// small, along with values of lower-left and upper-right pixels
const string  TEST_IMAGE_32x32("tests/testimage_expdisk32.fits");
const string  TARGET_FOR_NULL_IMAGE("nullimage_temp.fits");
const double  VALUE_LL = 20.903625;
const double  VALUE_UR = 14.947027;
// test images with no image data, or none in primary HDU
const string  TEST_FILE_EMPTYHDU("tests/test_emptyhdu.fits");
const string  TEST_FILE_TABLE("tests/test_table.fits");
// possibly good image (empty primary HDU, good image in 2nd HDU)
const string  TEST_IMAGE_MULTI_EMPTY_PRIMARYHDU("tests/test_multiextension_hdu0empty.fits");
// files with images in both primary HDU and 2nd HDU
const string  TEST_IMAGE_MULTI_2IMAGEHDUS("tests/test_multiextension_2images.fits");
const string  TEST_IMAGE_MULTI_2IMAGEHDUS2("tests/multiextension_ic3478rss_64x64_twice.fits");
// nonexistent file (for error-checking)
const string  BAD_IMAGE_NAME("no_image_with_this_name.fits");
// filename for an image file we can't possibly write (permission denied or
// at least directory doesn't exist)
const string  IMPOSSIBLE_IMAGE_FILENAME("/nopath/nohow/ZZZ/arglebargle.fits");

// simple image data and filename we will save it under
double  tinyImage[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
const string  TINY_IMAGE_FILENAME("tinyimage_temp.fits");



class BasicTestSuite : public CxxTest::TestSuite 
{
public:

  // Test for CheckForImage() with good file
  void testCheckForImage_good( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_IMAGE_32x32);    
    TS_ASSERT_EQUALS(response, 1);
  }

  // Test for CheckForImage() with nonexistent file
  void testCheckForImage_nofile( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(BAD_IMAGE_NAME);    
    TS_ASSERT_EQUALS(response, -1);
  }

  // Test for CheckForImage() with file containing only empty HDU
  void testCheckForImage_emptyfiles( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_FILE_EMPTYHDU);    
    TS_ASSERT_EQUALS(response, -1);
  }

  // Test for CheckForImage() with files containing only a table
  void testCheckForImage_tablefile( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_FILE_TABLE);    
    TS_ASSERT_EQUALS(response, -1);
  }

  // Test for CheckForImage() with multi-ext. containing empty primary HDU
  void testCheckForImage_multiext_emptyprimary( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_IMAGE_MULTI_EMPTY_PRIMARYHDU);    
    TS_ASSERT_EQUALS(response, -1);
  }

  // Test for CheckForImage() with multi-ext. containing empty primary HDU, specifying
  // 2nd HDU (valid image)
  void testCheckForImage_multiext_emptyprimary_hdu1( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_IMAGE_MULTI_EMPTY_PRIMARYHDU + "[1]");    
    TS_ASSERT_EQUALS(response, 1);
  }

  // Test for CheckForImage() with multi-ext. file (primary HDU = valid image)
  void testCheckForImage_multiext_hdu1( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_IMAGE_MULTI_2IMAGEHDUS);    
    TS_ASSERT_EQUALS(response, 1);
  }

  // Test for CheckForImage() with multi-ext. file (explicitly specifying primary HDU = valid image)
  void testCheckForImage_multiext_hdu1v2( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_IMAGE_MULTI_2IMAGEHDUS + "[0]");    
    TS_ASSERT_EQUALS(response, 1);
  }

  // Test for CheckForImage() with multi-ext. file (explicitly specifying second HDU = valid image)
  void testCheckForImage_multiext_hdu2( void )
  {
    int  response;
    int  nCols, nRows;

    response = CheckForImage(TEST_IMAGE_MULTI_2IMAGEHDUS + "[1]");    
    TS_ASSERT_EQUALS(response, 1);
  }



  // Test for GetImageSize() with good image file
  void testGetImageSize_good( void )
  {
    int  status;
    int  nCols, nRows;

    status = GetImageSize(TEST_IMAGE_32x32, &nCols, &nRows);    
    TS_ASSERT_EQUALS(nCols, 32);
    TS_ASSERT_EQUALS(nRows, 32);
  }

  // Test for GetImageSize() with nonexistent file
  void testGetImageSize_nofile( void )
  {
    int  status;
    int  nCols, nRows;

    status = GetImageSize(BAD_IMAGE_NAME, &nCols, &nRows);
    TS_ASSERT_EQUALS(-1, status);
  }

  // Test for GetImageSize() with empty FITS files (or empty primary HDU)
  void testGetImageSize_emptyfile( void )
  {
    int  status;
    int  nCols, nRows;

    status = GetImageSize(TEST_FILE_EMPTYHDU, &nCols, &nRows);
    TS_ASSERT_EQUALS(-1, status);

    status = GetImageSize(TEST_IMAGE_MULTI_EMPTY_PRIMARYHDU, &nCols, &nRows);
    TS_ASSERT_EQUALS(-1, status);
  }

  // Test for GetImageSize() with FITS table, not image
  void testGetImageSize_tablefile( void )
  {
    int  status;
    int  nCols, nRows;

    status = GetImageSize(TEST_FILE_TABLE, &nCols, &nRows);
    TS_ASSERT_EQUALS(-1, status);
  }


  // Test to see that ReadImageAsVector() succeeds
  void testReadImageAsVector_good( void )
  {
    int  nCols, nRows;
    double  *pixelData;

	pixelData = ReadImageAsVector(TEST_IMAGE_32x32, &nCols, &nRows);
	// check image size
    TS_ASSERT_EQUALS(nCols, 32);
    TS_ASSERT_EQUALS(nRows, 32);
    // check lower-left and upper-right pixel values
    TS_ASSERT_DELTA(pixelData[0], VALUE_LL, 1.0e-6);
    TS_ASSERT_DELTA(pixelData[32*32 - 1], VALUE_UR, 1.0e-6);
    
    free(pixelData);
  }

  // Test to see that ReadImageAsVector() correctly reports failure
  void testReadImageAsVector_failure( void )
  {
    int  nCols, nRows;
    double  *pixelData;

	pixelData = ReadImageAsVector(BAD_IMAGE_NAME, &nCols, &nRows);
	// check that error value (NULL) was returned
    TS_ASSERT(pixelData == NULL);
  }


  // Test round-trip of image data with SaveVectorAsImage() and ReadImageAsVector()
  void testReadAndWrite_noComments( void )
  {
    int nColsOrig, nRowsOrig, nCols, nRows, status;
    vector<string>  comments;
    double  *pixelData;
    
    nColsOrig = 3;
    nRowsOrig = 3;
    
    status = SaveVectorAsImage(tinyImage, TINY_IMAGE_FILENAME, nColsOrig, nRowsOrig, comments);
    TS_ASSERT_EQUALS(status, 0);

	pixelData = ReadImageAsVector(TINY_IMAGE_FILENAME, &nCols, &nRows);
	// check image size
    TS_ASSERT_EQUALS(nCols, 3);
    TS_ASSERT_EQUALS(nRows, 3);
    // check lower-left and upper-right pixel values
    TS_ASSERT_DELTA(pixelData[0], tinyImage[0], 1.0e-6);
    TS_ASSERT_DELTA(pixelData[8], tinyImage[8], 1.0e-6);
    
    free(pixelData);
  }

  // Same as previous, but now with comments for FITS header as well
  void testReadAndWrite_withComments( void )
  {
    int nColsOrig, nRowsOrig, nCols, nRows, status;
    vector<string>  comments;
    string  aComment;
    double  *pixelData;
    
    nColsOrig = 3;
    nRowsOrig = 3;
    aComment = "Here is a comment";
    comments.push_back(aComment);
    aComment = "this is a second comment";
    comments.push_back(aComment);
    
    status = SaveVectorAsImage(tinyImage, TINY_IMAGE_FILENAME, nColsOrig, nRowsOrig, comments);
    TS_ASSERT_EQUALS(status, 0);

	pixelData = ReadImageAsVector(TINY_IMAGE_FILENAME, &nCols, &nRows);
	// check image size
    TS_ASSERT_EQUALS(nCols, 3);
    TS_ASSERT_EQUALS(nRows, 3);
    // check lower-left and upper-right pixel values
    TS_ASSERT_DELTA(pixelData[0], tinyImage[0], 1.0e-6);
    TS_ASSERT_DELTA(pixelData[8], tinyImage[8], 1.0e-6);
    
    free(pixelData);
  }

  // Test to see that SaveVectorAsImage() correctly reports failure
  void testSaveImageAsVector_failure( void )
  {
    int  nCols, nRows, status;
    double  *nullImage = NULL;
    vector<string>  comments;
    
    nCols = 3;
    nRows = 3;

    status = SaveVectorAsImage(tinyImage, IMPOSSIBLE_IMAGE_FILENAME, nCols, nRows, comments);
	// check that error value was returned
    TS_ASSERT_EQUALS(status, -1);

    status = SaveVectorAsImage(nullImage, TARGET_FOR_NULL_IMAGE, nCols, nRows, comments);
	// check that error value was returned
    TS_ASSERT_EQUALS(status, -2);
  }

};


class MultiExtTestSuite : public CxxTest::TestSuite 
{
public:

  void testCheckHDUForImage_valid1( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_IMAGE_32x32.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 1, &status);
    TS_ASSERT_EQUALS(result, true);
    fits_close_file(imfile_ptr, &status);
  }

  void testCheckHDUForImage_valid2( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_IMAGE_MULTI_2IMAGEHDUS.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 1, &status);
    TS_ASSERT_EQUALS(result, true);
    fits_close_file(imfile_ptr, &status);
  }

  void testCheckHDUForImage_valid3( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_IMAGE_MULTI_EMPTY_PRIMARYHDU.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 2, &status);
    TS_ASSERT_EQUALS(result, true);
    fits_close_file(imfile_ptr, &status);
  }

  void testCheckHDUForImage_valid4( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_IMAGE_MULTI_2IMAGEHDUS.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 2, &status);
    TS_ASSERT_EQUALS(result, true);
    fits_close_file(imfile_ptr, &status);
  }


  void testCheckHDUForImage_invalid1( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_FILE_EMPTYHDU.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 1, &status);
    TS_ASSERT_EQUALS(result, false);    
    fits_close_file(imfile_ptr, &status);
  }

  void testCheckHDUForImage_invalid2( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_IMAGE_MULTI_EMPTY_PRIMARYHDU.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 1, &status);
    TS_ASSERT_EQUALS(result, false);    
    fits_close_file(imfile_ptr, &status);
  }

  void testCheckHDUForImage_invalid3( void )
  {
    fitsfile  *imfile_ptr;
    int  problems = 0;
    int  status = 0;
    int  imageHDU_num = -1;
    bool  result;

    problems = fits_open_file(&imfile_ptr, TEST_FILE_TABLE.c_str(), READONLY, &status);
    result = CheckHDUForImage(imfile_ptr, 1, &status);
    TS_ASSERT_EQUALS(result, false);    
    fits_close_file(imfile_ptr, &status);
  }
  
  
  
  // Test reading of 1st HDU in multi-ext image via "imagename.fits[1]"
  // Test for GetImageSize() with good file
  void testGetImageSize_hdu2image( void )
  {
    int  status;
    int  nCols, nRows;
    string  imageName;

	imageName = TEST_IMAGE_MULTI_2IMAGEHDUS + "[0]";
    status = GetImageSize(imageName, &nCols, &nRows);    
    TS_ASSERT_EQUALS(nCols, 5);
    TS_ASSERT_EQUALS(nRows, 5);

	imageName = TEST_IMAGE_MULTI_2IMAGEHDUS + "[1]";
    status = GetImageSize(imageName, &nCols, &nRows);    
    TS_ASSERT_EQUALS(nCols, 10);
    TS_ASSERT_EQUALS(nRows, 10);

	imageName = TEST_IMAGE_MULTI_2IMAGEHDUS2 + "[0]";
    status = GetImageSize(imageName, &nCols, &nRows);    
    TS_ASSERT_EQUALS(nCols,64);
    TS_ASSERT_EQUALS(nRows, 64);

	imageName = TEST_IMAGE_MULTI_2IMAGEHDUS2 + "[1]";
    status = GetImageSize(imageName, &nCols, &nRows);    
    TS_ASSERT_EQUALS(nCols, 64);
    TS_ASSERT_EQUALS(nRows, 64);
  }


};

