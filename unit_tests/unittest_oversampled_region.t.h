// Unit tests for oversampled_region.cpp
//
// This is not very useful except for testing the basic interface, because most of the
// real work is done internally, or requires a lot of surrounding setup first.
//
// cxxtestgen --error-printer -o test_runner_oversampled_region.cpp.cpp unit_tests/unittest_oversampled_region.h
// g++ -o test_runner test_runner_oversampled_region.cpp.cpp oversampled_region.cpp downsample.cpp convolver.cpp \
// image_io.cpp function_objects/function_object.cpp function_objects/func_gaussian.cpp \
// -I/usr/local/include -I$CXXTEST -lcfitsio -lfftw3 -lm
// ./test_runner

#include <cxxtest/TestSuite.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

#include "core/image_io.h"
#include "function_objects/function_object.h"
#include "function_objects/func_gaussian.h"
#include "core/oversampled_region.h"

#define DELTA  1.0e-10


// Need multiple PSFs, one for each oversampling scale
// These files are assumed to exist
string  psfImage_scale1_filename = string("unit_tests/psf_gauss_sigma3_35.fits");
string  psfImage_scale5_filename = string("unit_tests/psf_gauss_sigma15_175.fits");
string  psfImage_scale10_filename = string("unit_tests/psf_gauss_sigma30_351.fits");
// "target" image into which oversampled-and-downsampled regions will be copied
string  simpleNullImage51_filename = string("unit_tests/simpleimage_51x51_zeros.fits");



class TestOversampledRegion : public CxxTest::TestSuite 
{
  // data members
  int  nColsMain, nRowsMain;
  int  nColsOsamp, nRowsOsamp;
  int  oversamplingScale, nPSFColumns, nPSFRows;
  int  nOversampledModelColumns, nOversampledModelRows, nOversampledModelVals;
  int  nFunctions;
  bool  psfAllocated;
  double  *psfPixels;
  FunctionObject *gaussFunc;
  vector<FunctionObject *> functionObjects;
  
  OversampledRegion *oversampledRegion0;  // for simple testing
  
  OversampledRegion *oversampledRegion1a;  // for more complex testing w/o PSF
  OversampledRegion *oversampledRegion1b;  // for more complex testing w/o PSF
  OversampledRegion *oversampledRegion1c;  // for more complex testing w/o PSF

  OversampledRegion *oversampledRegion2a;  // for more complex testing w PSF
  OversampledRegion *oversampledRegion2b;  // for more complex testing w PSF
  OversampledRegion *oversampledRegion2c;  // for more complex testing w PSF

  
public:
  void setUp()
  {
    // get PSF pixel vector
    psfAllocated = false;
    psfPixels = ReadImageAsVector(psfImage_scale1_filename, &nPSFColumns, &nPSFRows);
    if (psfPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read PSF image file \"%s\"!\n\n", 
      			psfImage_scale1_filename.c_str());
      exit(-1);
    }
    psfAllocated = true;

    oversampledRegion0 = new OversampledRegion();
    oversampledRegion1a = new OversampledRegion();
    oversampledRegion1b = new OversampledRegion();
    oversampledRegion1c = new OversampledRegion();
    oversampledRegion2a = new OversampledRegion();
    oversampledRegion2b = new OversampledRegion();
    oversampledRegion2c = new OversampledRegion();
    
    gaussFunc = new Gaussian();
    functionObjects.push_back(gaussFunc);
    nFunctions = 1;
  }

  void tearDown()
  {
    delete oversampledRegion0;
    delete oversampledRegion1a;
    delete oversampledRegion1b;
    delete oversampledRegion1c;
    delete oversampledRegion2a;
    delete oversampledRegion2b;
    delete oversampledRegion2c;
    if (psfAllocated)
      free(psfPixels);
  }


  // and now the actual tests

  // Test to see if we can successfully add PSF vector to OversampledRegion object
  void testAddPSFVector( void )
  {
    oversampledRegion0->AddPSFVector(psfPixels, nPSFColumns, nPSFRows);
  }

  void testSetMaxThreads( void )
  {
    oversampledRegion0->SetMaxThreads(10);
  }

  void testSetupModelImage( void )
  {
    oversampledRegion0->SetupModelImage(1, 1, 10, 10, 500, 500,
    						31, 31, 3);
  }


  // And now the hard part: testing the main method
  // Currently these are a bit of a kludge, since we write the modified main image
  // to disk so we can inspect it with DS9 (instead of doing a proper "test" inside
  // the method here)
  
  // First, assume no PSFs at all
  void testMainStuff_NoPSF_osamp1( void )
  {
    int  oversampleScale = 1;
    int  x1 = 20, y1 = 20;
    int  delXmain = 11, delYmain = 11;
    double  params[] = {0.0, 0.0, 1.0, 2.0};
    
    double *mainImage = ReadImageAsVector(simpleNullImage51_filename, &nColsMain, &nRowsMain);
    oversampledRegion1a->SetupModelImage(x1, y1, delXmain, delYmain, nColsMain, nRowsMain, 0, 0, oversampleScale);

    // set up functions
    double  x0 = 25.0, y0 = 25.0;
    functionObjects[0]->Setup(params, 0, x0, y0);
    
    oversampledRegion1a->ComputeRegionAndDownsample(mainImage, functionObjects, 1);
    
    // Save mainImage to file
    int  status;
    vector<string>  imageCommentsList;
    string  imageName = string("simpleimage_51x51_pastedin_noPSF_osamp1.fits");
    printf("\nSaving output model image (\"%s\") ...\n", imageName.c_str());
    status = SaveVectorAsImage(mainImage, imageName, nColsMain, nRowsMain, imageCommentsList);

    free(mainImage);
  }

  // Same, but with higher oversampling
  void testMainStuff_NoPSF_osamp5( void )
  {
    int  oversampleScale = 5;
    int  x1 = 20, y1 = 20;
    int  delXmain = 11, delYmain = 11;
    double  params[] = {0.0, 0.0, 1.0, 2.0};
    
    double *mainImage = ReadImageAsVector(simpleNullImage51_filename, &nColsMain, &nRowsMain);
    oversampledRegion1b->SetupModelImage(x1, y1, delXmain, delYmain, nColsMain, nRowsMain, 0, 0, oversampleScale);

    // set up functions
    double  x0 = 25.0, y0 = 25.0;
    functionObjects[0]->Setup(params, 0, x0, y0);
    
    oversampledRegion1b->ComputeRegionAndDownsample(mainImage, functionObjects, 1);
    
    // Save mainImage to file
    int  status;
    vector<string>  imageCommentsList;
    string  imageName = string("simpleimage_51x51_pastedin_noPSF_osamp5.fits");
    printf("\nSaving output model image (\"%s\") ...\n", imageName.c_str());
    status = SaveVectorAsImage(mainImage, imageName, nColsMain, nRowsMain, imageCommentsList);

    free(mainImage);
  }

  // Same, but with higher oversampling
  void testMainStuff_NoPSF_osamp10( void )
  {
    int  oversampleScale = 10;
    int  x1 = 20, y1 = 20;
    int  delXmain = 11, delYmain = 11;
    double  params[] = {0.0, 0.0, 1.0, 2.0};
    
    printf("About to read image...\n");
    double *mainImage = ReadImageAsVector(simpleNullImage51_filename, &nColsMain, &nRowsMain);
    printf("Done.\n");
    oversampledRegion1c->SetupModelImage(x1, y1, delXmain, delYmain, nColsMain, nRowsMain, 0, 0, oversampleScale);

    // set up functions
    double  x0 = 25.0, y0 = 25.0;
    functionObjects[0]->Setup(params, 0, x0, y0);
    
    oversampledRegion1c->ComputeRegionAndDownsample(mainImage, functionObjects, 1);
    
    // Save mainImage to file
    int  status;
    vector<string>  imageCommentsList;
    string  imageName = string("simpleimage_51x51_pastedin_noPSF_osamp10.fits");
    printf("\nSaving output model image (\"%s\") ...\n", imageName.c_str());
    status = SaveVectorAsImage(mainImage, imageName, nColsMain, nRowsMain, imageCommentsList);

    free(mainImage);
  }


  // Now try it with PSF
  void testMainStuff_WithPSF_osamp1( void )
  {
    int  oversampleScale = 1;
    int  x1 = 20, y1 = 20;
    int  delXmain = 11, delYmain = 11;
    double  params[] = {0.0, 0.0, 1.0, 2.0};
    
    if (psfAllocated) {
      free(psfPixels);
      psfAllocated = false;
    }
    psfPixels = ReadImageAsVector(psfImage_scale1_filename, &nPSFColumns, &nPSFRows);
    if (psfPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read PSF image file \"%s\"!\n\n", 
      			psfImage_scale1_filename.c_str());
      exit(-1);
    }
    psfAllocated = true;

    double *mainImage = ReadImageAsVector(simpleNullImage51_filename, &nColsMain, &nRowsMain);
    oversampledRegion2a->AddPSFVector(psfPixels, nPSFColumns, nPSFRows);
    oversampledRegion2a->SetupModelImage(x1, y1, delXmain, delYmain, nColsMain, nRowsMain, 0, 0, oversampleScale);

    // set up functions
    double  x0 = 25.0, y0 = 25.0;
    functionObjects[0]->Setup(params, 0, x0, y0);
    
    string  outName = "oversampled_region_testoutput_withPSF_osamp1";
    oversampledRegion2a->SetDebugImageName(outName);
    oversampledRegion2a->SetDebugLevel(1);
    oversampledRegion2a->ComputeRegionAndDownsample(mainImage, functionObjects, 1);
    
    // Save mainImage to file
    int  status;
    vector<string>  imageCommentsList;
    string  imageName = string("simpleimage_51x51_pastedin_withPSF_osamp1.fits");
    printf("\nSaving output model image (\"%s\") ...\n", imageName.c_str());
    status = SaveVectorAsImage(mainImage, imageName, nColsMain, nRowsMain, imageCommentsList);

    free(mainImage);
  }

  void testMainStuff_WithPSF_osamp5( void )
  {
    int  oversampleScale = 5;
    int  x1 = 20, y1 = 20;
    int  delXmain = 11, delYmain = 11;
    double  params[] = {0.0, 0.0, 1.0, 2.0};
    
    if (psfAllocated) {
      free(psfPixels);
      psfAllocated = false;
    }
    psfPixels = ReadImageAsVector(psfImage_scale5_filename, &nPSFColumns, &nPSFRows);
    if (psfPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read PSF image file \"%s\"!\n\n", 
      			psfImage_scale5_filename.c_str());
      exit(-1);
    }
    psfAllocated = true;
    
    double *mainImage = ReadImageAsVector(simpleNullImage51_filename, &nColsMain, &nRowsMain);
    oversampledRegion2b->AddPSFVector(psfPixels, nPSFColumns, nPSFRows);
    oversampledRegion2b->SetupModelImage(x1, y1, delXmain, delYmain, nColsMain, nRowsMain, 0, 0, oversampleScale);

    // set up functions
    double  x0 = 25.0, y0 = 25.0;
    functionObjects[0]->Setup(params, 0, x0, y0);
    
    string  outName = "oversampled_region_testoutput_withPSF_osamp5";
    oversampledRegion2b->SetDebugImageName(outName);
    oversampledRegion2b->SetDebugLevel(1);
    oversampledRegion2b->ComputeRegionAndDownsample(mainImage, functionObjects, 1);
    
    // Save mainImage to file
    int  status;
    vector<string>  imageCommentsList;
    string  imageName = string("simpleimage_51x51_pastedin_withPSF_osamp5.fits");
    printf("\nSaving output model image (\"%s\") ...\n", imageName.c_str());
    status = SaveVectorAsImage(mainImage, imageName, nColsMain, nRowsMain, imageCommentsList);

    free(mainImage);
  }

  void testMainStuff_WithPSF_osamp10( void )
  {
    int  oversampleScale = 10;
    int  x1 = 20, y1 = 20;
    int  delXmain = 11, delYmain = 11;
    double  params[] = {0.0, 0.0, 1.0, 2.0};
    
    if (psfAllocated) {
      free(psfPixels);
      psfAllocated = false;
    }
    psfPixels = ReadImageAsVector(psfImage_scale10_filename, &nPSFColumns, &nPSFRows);
    if (psfPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read PSF image file \"%s\"!\n\n", 
      			psfImage_scale10_filename.c_str());
      exit(-1);
    }
    psfAllocated = true;
    
    double *mainImage = ReadImageAsVector(simpleNullImage51_filename, &nColsMain, &nRowsMain);
    oversampledRegion2c->AddPSFVector(psfPixels, nPSFColumns, nPSFRows);
    oversampledRegion2c->SetupModelImage(x1, y1, delXmain, delYmain, nColsMain, nRowsMain, 0, 0, oversampleScale);

    // set up functions
    double  x0 = 25.0, y0 = 25.0;
    functionObjects[0]->Setup(params, 0, x0, y0);
    
    string  outName = "oversampled_region_testoutput_withPSF_osamp10";
    oversampledRegion2c->SetDebugImageName(outName);
    oversampledRegion2c->SetDebugLevel(1);
    oversampledRegion2c->ComputeRegionAndDownsample(mainImage, functionObjects, 1);
    
    // Save mainImage to file
    int  status;
    vector<string>  imageCommentsList;
    string  imageName = string("simpleimage_51x51_pastedin_withPSF_osamp10.fits");
    printf("\nSaving output model image (\"%s\") ...\n", imageName.c_str());
    status = SaveVectorAsImage(mainImage, imageName, nColsMain, nRowsMain, imageCommentsList);

    free(mainImage);
  }

};

