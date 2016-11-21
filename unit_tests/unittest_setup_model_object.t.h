// Unit tests for code in setup_model_object.cpp

// Evidence for the idea that we're now doing integration tests is the large
// number of other modules that have to be compiled along with model_object.cpp
// to get these tests to work...

// See run_unittest_setup_model_object.sh for how to compile & run this


// ModelObject* SetupModelObject( OptionsBase *options, vector<int> nColumnsRowsVector, 
// 					double *dataPixels, double *psfPixels, double *maskPixels, 
// 					double *errorPixels, double *psfOversampledPixels, 
// 					vector<int> xyOversamplePos )

// We assume that the input nColumnsRowsVector has the following entries:
// nColumnsRowsVector[0] = nColumns
// nColumnsRowsVector[1] = nRows
// nColumnsRowsVector[2] = nColumns_psf  [optional]
// nColumnsRowsVector[3] = nRows_psf  [optional]
// nColumnsRowsVector[4] = nColumns_psf_oversampled  [optional]
// nColumnsRowsVector[5] = nRows_psf_oversampled  [optional]


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
using namespace std;
#include "definitions.h"
//#include "function_objects/function_object.h"
#include "setup_model_object.h"
#include "model_object.h"
//#include "add_functions.h"
#include "config_file_parser.h"
#include "options_base.h"
#include "options_makeimage.h"
#include "options_imfit.h"
#include "options_mcmc.h"


#define SIMPLE_CONFIG_FILE "tests/config_imfit_flatsky.dat"
#define CONFIG_FILE "tests/config_imfit_poisson.dat"


// Reference things
const string  headerLine_correct = "# X0_1		Y0_1		PA_1	ell_1	I_0_1	h_1	I_sky_2	";


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  int  status;
  configOptions  userConfigOptions1, userConfigOptions3;
  double  *smallDataImage;
  double  *smallErrorImage;
  double  *smallVarianceImage;
  double  *smallWeightImage;
  double  *smallMaskImage;
  double  *smallPSFImage;
  double  *oversampledPSFImage;
  int nSmallDataCols, nSmallDataRows;
  int nSmallPSFCols, nSmallPSFRows;
  int nOsampPSFCols, nOsampPSFRows;


  // Note that setUp() gets called prior to *each* individual test function!
  void setUp()
  {
    int  status;
    string  filename1 = CONFIG_FILE;
    string  filename3 = SIMPLE_CONFIG_FILE;
    
    nSmallDataCols = nSmallDataRows = 2;
    nSmallPSFCols = nSmallPSFRows = 2;
    nOsampPSFCols = nOsampPSFRows = 4;
    
    smallDataImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallDataImage[0] = 0.25;
    smallDataImage[1] = 0.25;
    smallDataImage[2] = 0.25;
    smallDataImage[3] = 1.0;

    // small image corresponding to square root of smallDataImage
    smallErrorImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallErrorImage[0] = 0.5;
    smallErrorImage[1] = 0.5;
    smallErrorImage[2] = 0.5;
    smallErrorImage[3] = 1.0;

    smallMaskImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallMaskImage[0] = 0.0;
    smallMaskImage[1] = 1.0;
    smallMaskImage[2] = 0.0;
    smallMaskImage[3] = 0.0;

    smallPSFImage = (double *)calloc(nSmallPSFCols*nSmallPSFRows, sizeof(double));
    smallPSFImage[0] = 0.0;
    smallPSFImage[1] = 1.0;
    smallPSFImage[2] = 0.0;
    smallPSFImage[3] = 0.0;
    
    // just use a zero-filled image to keep things simple
    oversampledPSFImage = (double *)calloc(nOsampPSFCols*nOsampPSFRows, sizeof(double));
  }

  void tearDown()
  {
    free(smallDataImage);
    free(smallErrorImage);
    free(smallWeightImage);
    free(smallMaskImage);
    free(smallPSFImage);
    free(oversampledPSFImage);
  }
  
  
  void testSetupMakeimage_simple( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new MakeimageOptions();
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupMakeimage_withPSF( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;

    long nDataVals_true = 4;
    long nDataVals;
  
    optionsPtr = new MakeimageOptions();
    optionsPtr->psfImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL, smallPSFImage);
  
    nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupMakeimage_withOversampledPSF( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
    vector<int> xyOsamplePos;

    long nDataVals_true = 4;
    long nDataVals;
  
    optionsPtr = new MakeimageOptions();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampledImagePresent = true;
    optionsPtr->psfOversamplingScale = 2;

    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
    nColumnsRowsVect.push_back(nOsampPSFCols);
    nColumnsRowsVect.push_back(nOsampPSFRows);
    xyOsamplePos.push_back(1);
    xyOsamplePos.push_back(2);
    xyOsamplePos.push_back(1);
    xyOsamplePos.push_back(2);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL, smallPSFImage,
    					NULL, NULL, oversampledPSFImage, xyOsamplePos);
  
    nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);

    delete optionsPtr;
    delete theModel;
  }
 
 
  void testSetupImfit_simple( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new ImfitOptions();
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupImfit_withPSF( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new ImfitOptions();
    optionsPtr->psfImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			smallPSFImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupImfit_withOversampledPSF( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
    vector<int> xyOsamplePos;

    optionsPtr = new ImfitOptions();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampledImagePresent = true;
    optionsPtr->psfOversamplingScale = 2;

    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
    nColumnsRowsVect.push_back(nOsampPSFCols);
    nColumnsRowsVect.push_back(nOsampPSFRows);
    xyOsamplePos.push_back(1);
    xyOsamplePos.push_back(2);
    xyOsamplePos.push_back(1);
    xyOsamplePos.push_back(2);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage, smallPSFImage,
    					NULL, NULL, oversampledPSFImage, xyOsamplePos);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupImfit_withMask( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new ImfitOptions();
    optionsPtr->maskImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, smallMaskImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, true);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupImfit_withError( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new ImfitOptions();
    optionsPtr->noiseImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, NULL, smallErrorImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete optionsPtr;
    delete theModel;
  }


  void testSetupMCMC_simple( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new MCMCOptions();
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);
	
    delete optionsPtr;
    delete theModel;
  }

  void testSetupMCMC_withPSF( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new MCMCOptions();
    optionsPtr->psfImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			smallPSFImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupMCMC_withOversampledPSF( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
    vector<int> xyOsamplePos;

    optionsPtr = new MCMCOptions();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampledImagePresent = true;
    optionsPtr->psfOversamplingScale = 2;

    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
    nColumnsRowsVect.push_back(nOsampPSFCols);
    nColumnsRowsVect.push_back(nOsampPSFRows);
    xyOsamplePos.push_back(1);
    xyOsamplePos.push_back(2);
    xyOsamplePos.push_back(1);
    xyOsamplePos.push_back(2);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage, smallPSFImage,
    					NULL, NULL, oversampledPSFImage, xyOsamplePos);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupMCMC_withMask( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new MCMCOptions();
    optionsPtr->maskImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, smallMaskImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, true);

    delete optionsPtr;
    delete theModel;
  }

  void testSetupMCMC_withError( void )
  {
    ModelObject *theModel = nullptr;
    OptionsBase *optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = new MCMCOptions();
    optionsPtr->noiseImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, NULL, smallErrorImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete optionsPtr;
    delete theModel;
  }

};
