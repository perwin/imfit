// Unit tests for code in model_object.cpp
// Note that this could be considered a case of unit tests merging into the
// realm of "integration testing", since so much of ModeObject's code requires
// successful interaction with code from other modules.

// Evidence for the idea that we're now doing integration tests is the large
// number of other modules that have to be compiled along with model_object.cpp
// to get these tests to work...

// See run_unittest_model_object.sh for how to compile & run this


// Things to test:


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
#include <stdlib.h>
using namespace std;
#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object.h"
#include "add_functions.h"
#include "config_file_parser.h"
#include "param_struct.h"


#define SIMPLE_CONFIG_FILE "tests/config_imfit_flatsky.dat"
#define CONFIG_FILE "tests/config_imfit_poisson.dat"


// Reference things
const string  headerLine_correct = "# X0_1		Y0_1		PA_1	ell_1	I_0_1	h_1	I_sky_2	";


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  ModelObject *modelObj1;
  ModelObject *modelObj2a;
  ModelObject *modelObj2b;
  ModelObject *modelObj2c;
  ModelObject *modelObj2d;
  ModelObject *modelObj2e;
  ModelObject *modelObj2f;
  ModelObject *modelObj3;
  ModelObject *modelObj4a;
  ModelObject *modelObj4b;
  vector<string>  functionList1, functionList3;
  vector<double>  parameterList1, parameterList3;
  vector<mp_par>  paramLimits1, paramLimits3;
  vector<int>  FunctionBlockIndices1, FunctionBlockIndices3;
  bool  paramLimitsExist1, paramLimitsExist3;
  mp_par  *parameterInfo;
  int  status;
  configOptions  userConfigOptions1, userConfigOptions3;
  double  *smallDataImage;
  double  *smallErrorImage;
  double  *smallVarianceImage;
  double  *smallWeightImage;
  double  *smallMaskImage;
  int nSmallDataCols, nSmallDataRows;


  void setUp()
  {
    int  status;
    string  filename1 = CONFIG_FILE;
    string  filename3 = SIMPLE_CONFIG_FILE;
    
    nSmallDataCols = nSmallDataRows = 2;
    
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

	// variance image
    smallVarianceImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallVarianceImage[0] = 0.25;
    smallVarianceImage[1] = 0.25;
    smallVarianceImage[2] = 0.25;
    smallVarianceImage[3] = 1.0;

	// weight = inverse variance image
    smallWeightImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallWeightImage[0] = 4.0;
    smallWeightImage[1] = 4.0;
    smallWeightImage[2] = 4.0;
    smallWeightImage[3] = 1.0;

    smallMaskImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallMaskImage[0] = 0.0;
    smallMaskImage[1] = 1.0;
    smallMaskImage[2] = 0.0;
    smallMaskImage[3] = 0.0;


    modelObj1 = new ModelObject();    
    // Initialize modelObj1; set up internal FunctionObjects vector (Exp + FlatSky) inside
    status = ReadConfigFile(filename1, true, functionList1, parameterList1, 
  								paramLimits1, FunctionBlockIndices1, paramLimitsExist1, userConfigOptions1);
    status = AddFunctions(modelObj1, functionList1, FunctionBlockIndices1, true, -1);

    // Initialize modelObj2a, etc.
    modelObj2a = new ModelObject();
    modelObj2b = new ModelObject();
    modelObj2c = new ModelObject();
    modelObj2d = new ModelObject();
    modelObj2e = new ModelObject();
    modelObj2f = new ModelObject();
    
    
    // Initialize modelObj3 and add model function & params (FlatSky)
    status = ReadConfigFile(filename3, true, functionList3, parameterList3, 
  								paramLimits3, FunctionBlockIndices3, paramLimitsExist3, userConfigOptions3);
    modelObj3 = new ModelObject();
    status = AddFunctions(modelObj3, functionList3, FunctionBlockIndices3, true, -1);

    // Initialize modelObj4a and add model function & params (FlatSky)
    modelObj4a = new ModelObject();
    status = AddFunctions(modelObj4a, functionList3, FunctionBlockIndices3, true, -1);
    // Initialize modelObj4b and add model function & params (FlatSky)
    modelObj4b = new ModelObject();
    status = AddFunctions(modelObj4b, functionList3, FunctionBlockIndices3, true, -1);
    
  }

  void tearDown()
  {
    free(smallDataImage);
    free(smallErrorImage);
    free(smallVarianceImage);
    free(smallWeightImage);
    free(smallMaskImage);
    delete modelObj1;
    delete modelObj2a;
    delete modelObj2b;
    delete modelObj2c;
    delete modelObj2d;
    delete modelObj2e;
    delete modelObj2f;
    delete modelObj3;
    delete modelObj4a;
    delete modelObj4b;
  }
  
  
  void testGetParamHeader( void )
  {
    std::string  headerLine;
    headerLine = modelObj1->GetParamHeader();
    TS_ASSERT_EQUALS(headerLine, headerLine_correct);

  }
  
  void testSetAndGetStatisticType( void )
  {
    int  whichStat, status;
    bool cashStatUsed = false;
    
    // default should be chi^2
    whichStat = modelObj1->WhichFitStatistic();
    TS_ASSERT_EQUALS(whichStat, FITSTAT_CHISQUARE);
    // special case of verbose output
    whichStat = modelObj1->WhichFitStatistic(true);
    TS_ASSERT_EQUALS(whichStat, FITSTAT_CHISQUARE_DATA);
    // switch to using model-based chi^2
    status = modelObj1->UseModelErrors();
    whichStat = modelObj1->WhichFitStatistic(true);
    TS_ASSERT_EQUALS(whichStat, FITSTAT_CHISQUARE_MODEL);
    
    modelObj1->UseCashStatistic();
    whichStat = modelObj1->WhichFitStatistic();
    TS_ASSERT_EQUALS(whichStat, FITSTAT_CASH);
    cashStatUsed = modelObj1->UsingCashStatistic();
    TS_ASSERT_EQUALS(cashStatUsed, true);

    // the following will generate a couple of warnings from within UseCashStatistic,
    // which is OK
    modelObj1->UsePoissonMLR();
    whichStat = modelObj1->WhichFitStatistic();
    TS_ASSERT_EQUALS(whichStat, FITSTAT_POISSON_MLR);
    cashStatUsed = modelObj1->UsingCashStatistic();
    TS_ASSERT_EQUALS(cashStatUsed, true);
  }
 
 
  void testStoreAndRetrieveDataImage( void )
  {
    double *outputVect;
    double trueVals[4] = {0.25, 0.25, 0.25, 1.0};   // data values
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2a->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    outputVect = modelObj2a->GetDataVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputVect[i], trueVals[i]);
  }

  void testGenerateAndRetrieveErrorImage( void )
  {
    double *outputVect;
    double trueVals[4] = {4.0, 4.0, 4.0, 1.0};   // 1/sigma^2
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2b->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    modelObj2b->GenerateErrorVector();
    outputVect = modelObj2b->GetWeightImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputVect[i], trueVals[i]);
  }
 
  void testStoreAndRetrieveErrorImage( void )
  {
    double *outputVect;
    double trueVals[4] = {4.0, 4.0, 4.0, 1.0};   // 1/sigma^2
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2c->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    modelObj2c->AddErrorVector(nDataVals, nSmallDataCols, nSmallDataRows, smallErrorImage, WEIGHTS_ARE_SIGMAS);
    outputVect = modelObj2c->GetWeightImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputVect[i], trueVals[i]);
  }
 
  void testStoreAndRetrieveVarianceImage( void )
  {
    double *outputVect;
    double trueVals[4] = {4.0, 4.0, 4.0, 1.0};   // 1/sigma^2
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2d->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    modelObj2d->AddErrorVector(nDataVals, nSmallDataCols, nSmallDataRows, smallVarianceImage, WEIGHTS_ARE_VARIANCES);
    outputVect = modelObj2d->GetWeightImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputVect[i], trueVals[i]);
  }
 
  void testStoreAndRetrieveWeightImage( void )
  {
    double *outputVect;
    double trueVals[4] = {4.0, 4.0, 4.0, 1.0};   // 1/sigma^2
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2e->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    modelObj2e->AddErrorVector(nDataVals, nSmallDataCols, nSmallDataRows, smallWeightImage, WEIGHTS_ARE_WEIGHTS);
    outputVect = modelObj2e->GetWeightImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputVect[i], trueVals[i]);
  }

  void testStoreAndApplyMaskImage( void )
  {
    double *outputVect;
    double trueVals[4] = {4.0, 0.0, 4.0, 1.0};   // 1/sigma^2, with second pixel masked out
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2f->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    modelObj2f->AddErrorVector(nDataVals, nSmallDataCols, nSmallDataRows, smallWeightImage, WEIGHTS_ARE_WEIGHTS);
	modelObj2f->AddMaskVector(nDataVals, nSmallDataCols, nSmallDataRows, smallMaskImage, MASK_ZERO_IS_GOOD);
	modelObj2f->ApplyMask();
    outputVect = modelObj2f->GetWeightImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputVect[i], trueVals[i]);
  }

 
 
   void testModelImageGeneration( void )
  {
    // Simple model image: 4x4 pixels, FlatSky function with I_sky = 100.0
    double *outputModelVect;
    double params[3] = {26.0, 26.0, 100.0};   // X0, Y0, I_sky
    double trueVals[4] = {100.0, 100.0, 100.0, 100.0};
    int  nDataVals = nSmallDataCols*nSmallDataRows;
  
    modelObj3->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj3->CreateModelImage(params);
    outputModelVect = modelObj3->GetModelImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputModelVect[i], trueVals[i]);
  }
 
  // make sure ModelObject complains if we add a PSF *after* we've done the setup
  void testCatchOutOfOrderPSF( void )
  {
    int  nColumns = 10;
    int  nRows = 10;
    int  nColumns_psf = 3;
    int  nRows_psf = 3;
    int  nPixels_psf = 9;
    double  smallPSFImage[9] = {0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
    int  status4a, status4b;
    
    // This is the correct order
    // add PSF pixels first
    status4a = modelObj4a->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, smallPSFImage);
    // final setup for modelObj4a
    modelObj4a->SetupModelImage(nColumns, nRows);
    TS_ASSERT_EQUALS(status4a, 0);

	// This is the WRONG order
    // final setup for modelObj4b *first*
    modelObj4b->SetupModelImage(nColumns, nRows);
    // then add PSF pixels = wrong order!
    status4b = modelObj4b->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, smallPSFImage);
    TS_ASSERT_EQUALS(status4b, -1);
  }
};
