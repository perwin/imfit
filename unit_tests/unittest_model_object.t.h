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
#include <math.h>
using namespace std;
#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object.h"
#include "add_functions.h"
#include "config_file_parser.h"
#include "param_struct.h"


#define SIMPLE_CONFIG_FILE "tests/imfit_reference/config_imfit_flatsky.dat"
#define CONFIG_FILE "tests/imfit_reference/config_imfit_poisson.dat"
#define CONFIG_FILE_POINTSOURCE "tests/makeimage_reference/config_imfit_pointsource.dat"
#define CONFIG_FILE3b "tests/config_imfit_gauss-extra-params.dat"
const string CONFIG_FILE_2BLOCKS("tests/config_imfit_2sets.dat");


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
  ModelObject *modelObj3a;
  ModelObject *modelObj3b;
  ModelObject *modelObj3c;
  ModelObject *modelObj3d;
  ModelObject *modelObj4;
  vector<string>  functionList1, functionList3, functionList3b;
  vector<string>  functionLabelList1, functionLabelList3, functionLabelList3b;
  vector<double>  parameterList1, parameterList3, parameterList3b;
  vector<mp_par>  paramLimits1, paramLimits3, paramLimits3b;
  vector<int>  FunctionBlockIndices1, FunctionBlockIndices3, FunctionBlockIndices3b;
  bool  paramLimitsExist1, paramLimitsExist3, paramLimitsExist3b;
  mp_par  *parameterInfo;
  int  status;
  configOptions  userConfigOptions1, userConfigOptions3, userConfigOptions3b;
  double  *smallDataImage;
  double  *smallErrorImage;
  double  *smallVarianceImage;
  double  *smallWeightImage;
  double  *smallMaskImage;
  int nSmallDataCols, nSmallDataRows;


  // Note that setUp() gets called prior to *each* individual test function!
  void setUp()
  {
    int  status;
    string  filename1 = CONFIG_FILE;
    string  filename3 = SIMPLE_CONFIG_FILE;
    string  filename3b = CONFIG_FILE3b;
    
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
    status = ReadConfigFile(filename1, true, functionList1, functionLabelList1, parameterList1, 
  								paramLimits1, FunctionBlockIndices1, paramLimitsExist1, 
  								userConfigOptions1);
    status = AddFunctions(modelObj1, functionList1, functionLabelList1, FunctionBlockIndices1, 
    						true, -1);

    // Initialize modelObj2a, etc.
    modelObj2a = new ModelObject();
    modelObj2b = new ModelObject();
    modelObj2c = new ModelObject();
    modelObj2d = new ModelObject();
    modelObj2e = new ModelObject();
    modelObj2f = new ModelObject();
    
    
    // Initialize modelObj3a and add model function & params (FlatSky)
    status = ReadConfigFile(filename3, true, functionList3, functionLabelList3, parameterList3, 
  								paramLimits3, FunctionBlockIndices3, paramLimitsExist3, 
  								userConfigOptions3);
    modelObj3a = new ModelObject();
    status = AddFunctions(modelObj3a, functionList3, functionLabelList3, FunctionBlockIndices3, 
    						true, -1);
    // Initialize modelObj3b and add model function & params (GaussianExtraParams)
    // turn off subsampling to make things quicker and calculations simpler
    status = ReadConfigFile(filename3b, true, functionList3b, functionLabelList3b, parameterList3b, 
  								paramLimits3b, FunctionBlockIndices3b, paramLimitsExist3b, 
  								userConfigOptions3b);
    modelObj3b = new ModelObject();
    status = AddFunctions(modelObj3b, functionList3b, functionLabelList3b, 
    						FunctionBlockIndices3b, false, -1);
    modelObj3c = new ModelObject();
    status = AddFunctions(modelObj3c, functionList3, functionLabelList3, 
    						FunctionBlockIndices3, true, -1);

    // Initialize modelObj4 and add model function & params (Exp + FlatSky)
    modelObj4 = new ModelObject();
    status = AddFunctions(modelObj4, functionList1, functionLabelList1, FunctionBlockIndices1, 
    						true, -1);
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
    delete modelObj3a;
    delete modelObj3b;
    delete modelObj3c;
    delete modelObj4;
  }
  
  
  void testGetParamHeader( void )
  {
    std::string  headerLine;
    headerLine = modelObj1->GetParamHeader();
    TS_ASSERT_EQUALS(headerLine, headerLine_correct);
  }
  

  void testGetNFunctions( void )
  {
    int  nFuncs_true = 2;   // Exponential,FlatSky
    int  nFuncs = modelObj1->GetNFunctions();
    TS_ASSERT_EQUALS(nFuncs, nFuncs_true);
    
    nFuncs_true = 1;   // GaussianExtraParams
    nFuncs = modelObj3b->GetNFunctions();
    TS_ASSERT_EQUALS(nFuncs, nFuncs_true);
  }
  
  void testGetNParams( void )
  {
    int  nParams_true = 7;   // X0,Y0,Exponential,FlatSky
    int  nParams = modelObj1->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParams_true);
    
    nParams_true = 6;   // GaussianExtraParams
    nParams = modelObj3b->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParams_true);
  }
  
  void testStoreAndRetrieveDataImage( void )
  {
    double *outputVect;
    double trueVals[4] = {0.25, 0.25, 0.25, 1.0};   // data values
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    
    modelObj2a->AddImageDataVector(smallDataImage, nSmallDataCols, nSmallDataRows);
    outputVect = modelObj2a->GetDataVector();

    int nData = modelObj2a->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nData);
    
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

    bool maskPresent = modelObj2f->HasMask();
    TS_ASSERT_EQUALS(maskPresent, false);
	modelObj2f->AddMaskVector(nDataVals, nSmallDataCols, nSmallDataRows, smallMaskImage, MASK_ZERO_IS_GOOD);
    maskPresent = modelObj2f->HasMask();
    TS_ASSERT_EQUALS(maskPresent, true);

	modelObj2f->ApplyMask();
    outputVect = modelObj2f->GetWeightImageVector();

    int  n = modelObj2f->GetNDataValues();
    TS_ASSERT_EQUALS(n, nDataVals);
    
    int  nValidPix_true = nDataVals - 1;
    n = modelObj2f->GetNValidPixels();
    TS_ASSERT_EQUALS(n, nValidPix_true);
    
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
  
    modelObj3a->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj3a->CreateModelImage(params);
    outputModelVect = modelObj3a->GetModelImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputModelVect[i], trueVals[i]);
  }


   void testResidualImageGeneration( void )
  {
    // Simple model image: 4x4 pixels, FlatSky function with I_sky = 100.0
    double *outputModelVect;
    double params[3] = {26.0, 26.0, 100.0};   // X0, Y0, I_sky
    double trueVals[4] = {100.0, 100.0, 100.0, 100.0};
    double dataVals[4] = {101.0, 99.0, 100.1, 100.0};
    double trueResidualVals[4] = {101.0 - 100.0, 99.0 - 100.0, 100.1 - 100.0, 100.0 - 100.0};
    int  nDataVals = nSmallDataCols*nSmallDataRows;
  
    int status = modelObj3c->AddImageDataVector(dataVals, nSmallDataCols, nSmallDataRows);
    TS_ASSERT_EQUALS(status, 0);
    
    modelObj3c->CreateModelImage(params);

    outputModelVect = modelObj3c->GetModelImageVector();
    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputModelVect[i], trueVals[i]);

    outputModelVect = modelObj3c->GetResidualImageVector();
    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_EQUALS(outputModelVect[i], trueResidualVals[i]);
  }
 
 
  void testSetExtraParams( void )
  {
    // model image: 4x4 pixels, GaussianExtraParams function with center at (x,y) = (1,1)
    ModelObject *modelObj;
    double *outputModelVect;
    double params[6] = {1.0, 1.0, 0.0, 0.0, 100.0, 1.0};   // X0, Y0, PA, ell, I_0, sigma
    double r1val = 60.653065971263345;
    double r2val = 36.787944117144235;
    double trueVals[4] = {100.0, r1val, r1val, r2val};
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  status;
 
    // sanity check: do we generate proper model when optional parameter
    // floor = 0?
    modelObj3b->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj3b->CreateModelImage(params);
    outputModelVect = modelObj3b->GetModelImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_DELTA(outputModelVect[i], trueVals[i], 1e-7);
    

    // test with floor = 100.0 [added to each output pixel]
    vector< map<string, string> > optionalParamsVect;
    map<string, string> optionalParamsMap;
    string  keyword = "floor";
    string  value = "100";
    double  floorVal = 100.0;
    optionalParamsMap[keyword] = value;
    optionalParamsVect.push_back(optionalParamsMap);

    vector<string> funcList;
    vector<string> funcLabelList;
    vector<double> paramList;
    vector<mp_par> paramLimits;
    vector<int> funcBlockIndices;
    bool paramLimitsExist;
    configOptions userConfigOptions;
    status = ReadConfigFile(CONFIG_FILE3b, true, funcList, funcLabelList, paramList, paramLimits, 
    							funcBlockIndices, paramLimitsExist, userConfigOptions);

    modelObj = new ModelObject();
    status = AddFunctions(modelObj, funcList, funcLabelList, funcBlockIndices, false, -1,
    						optionalParamsVect);
    modelObj->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj->CreateModelImage(params);
    outputModelVect = modelObj->GetModelImageVector();

    for (int i = 0; i < nDataVals; i++)
      TS_ASSERT_DELTA(outputModelVect[i], trueVals[i] + floorVal, 1e-7);
    delete modelObj;
 }
};



class TestPrintParams : public CxxTest::TestSuite
{
public:
  ModelObject * modelObj1;
  ModelObject * modelObj2a;
  ModelObject * modelObj2b;
  ModelObject * modelObj3;
  vector<string>  functionList1, functionList2a, functionList3;
  vector<string>  functionLabelList1, functionLabelList2a, functionLabelList3;
  vector<double>  parameterList1, parameterList2a, parameterList3;
  vector<mp_par>  paramLimits1, paramLimits2a, paramLimits3;
  vector<int>  FunctionBlockIndices1, FunctionBlockIndices2a, FunctionBlockIndices3;
  bool  paramLimitsExist1, paramLimitsExist2a, paramLimitsExist3;
  configOptions  userConfigOptions1, userConfigOptions2a, userConfigOptions3;
  int nSmallDataCols, nSmallDataRows;
  
  void setUp( )
  {
    int  status;
    string  filename1 = CONFIG_FILE;
    string  filename2a = SIMPLE_CONFIG_FILE;
    string  filename3 = CONFIG_FILE_2BLOCKS;

    nSmallDataCols = nSmallDataRows = 2;
    
    status = ReadConfigFile(filename1, true, functionList1, functionLabelList1, parameterList1, 
    						paramLimits1, FunctionBlockIndices1, paramLimitsExist1, 
    						userConfigOptions1);
    // Initialize modelObj1 and add model function & params (Exp + FlatSky)
    modelObj1 = new ModelObject();
    status = AddFunctions(modelObj1, functionList1, functionLabelList1, FunctionBlockIndices1, 
    						true, -1);

    status = ReadConfigFile(filename2a, true, functionList2a, functionLabelList2a, parameterList2a, 
    						paramLimits2a, FunctionBlockIndices2a, paramLimitsExist2a, 
    						userConfigOptions2a);
    // Initialize modelObj2a and add model function & params (FlatSky)
    modelObj2a = new ModelObject();
    status = AddFunctions(modelObj2a, functionList2a, functionLabelList2a, FunctionBlockIndices2a, 
    						true, -1);

    modelObj2b = new ModelObject();
    status = AddFunctions(modelObj2b, functionList1, functionLabelList1, FunctionBlockIndices1, 
    						true, -1);

    // Exp + FlatSky, each in its own function set
    status = ReadConfigFile(filename3, true, functionList3, functionLabelList3, parameterList3, 
    						paramLimits3, FunctionBlockIndices3, paramLimitsExist3, 
    						userConfigOptions3);
    modelObj3 = new ModelObject();
    status = AddFunctions(modelObj3, functionList3, functionLabelList3, FunctionBlockIndices3, 
    						true, -1);
  }
  
  void tearDown( )
  {
    delete modelObj1;
    delete modelObj2a;
    delete modelObj2b;
    delete modelObj3;
  }
  
  void testPrintParamsToString( void )
  {
    int nParamsTot = 7;
    double params[7] = {21.0, 22.0, 0.0, 0.5, 50.0, 10.0, 100.0};   // X0, Y0, PA, ell, I_0, h, I_sky
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  retVal;
    vector<string> outputVect;
  
    vector<string> correctStrings1, correctStrings2;
    correctStrings1.push_back("#\n");
    correctStrings1.push_back("#X0\t\t21.0000\n");
    correctStrings1.push_back("#Y0\t\t22.0000\n");
    correctStrings1.push_back("#FUNCTION Exponential\n");
    correctStrings1.push_back("#PA\t\t      0\n");
    correctStrings1.push_back("#ell\t\t    0.5\n");
    correctStrings1.push_back("#I_0\t\t     50\n");
    correctStrings1.push_back("#h\t\t     10\n");
    correctStrings1.push_back("#FUNCTION FlatSky\n");
    correctStrings1.push_back("#I_sky\t\t    100\n");
    correctStrings2.push_back("#\n");
    correctStrings2.push_back("#X0\t\t21.0000\t\t0,200\n");
    correctStrings2.push_back("#Y0\t\t22.0000\t\tfixed\n");
    correctStrings2.push_back("#FUNCTION Exponential\n");
    correctStrings2.push_back("#PA\t\t      0\t\t2,202\n");
    correctStrings2.push_back("#ell\t\t    0.5\t\t3,203\n");
    correctStrings2.push_back("#I_0\t\t     50\t\t4,204\n");
    correctStrings2.push_back("#h\t\t     10\t\tfixed\n");
    correctStrings2.push_back("#FUNCTION FlatSky\n");
    correctStrings2.push_back("#I_sky\t\t    100\t\t6,206\n");

    string prefix = "#";
    vector<mp_par> parameterInfo_empty;
    vector<mp_par> parameterInfo;

    // kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 1;
      currentParameterInfo.limited[1] = 1;
      currentParameterInfo.limits[0] = 0.0 + i;
      currentParameterInfo.limits[1] = 200.0 + i;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }
    // specify that Y0 and Exponential h are fixed
    parameterInfo[1].fixed = 1;
    parameterInfo[5].fixed = 1;

    modelObj1->SetupModelImage(nSmallDataCols, nSmallDataRows);
  
    // Check that we correctly catch empty parameterInfo when we request printing
    // of parameter limits
    printf("size of parameterInfo = %d\n", (int)parameterInfo.size());
    retVal = modelObj1->PrintModelParamsToStrings(outputVect, params, NULL, prefix.c_str(), true);
    TS_ASSERT_EQUALS(retVal, -1);
    

    modelObj1->AddParameterInfo(parameterInfo);

    // output without parameter limits
    retVal = modelObj1->PrintModelParamsToStrings(outputVect, params, NULL, prefix.c_str(), false);
    TS_ASSERT_EQUALS(retVal, 0);
    for (int i = 0; i < 9; i++) {
      TS_ASSERT_EQUALS(outputVect[i], correctStrings1[i]);
    }

    // output *with* parameter limits
    outputVect.clear();
    retVal = modelObj1->PrintModelParamsToStrings(outputVect, params, NULL, prefix.c_str(), true);
    TS_ASSERT_EQUALS(retVal, 0);
    for (int i = 0; i < 9; i++) {
      TS_ASSERT_EQUALS(outputVect[i], correctStrings2[i]);
    }
  }

  void testPrintModelParamsHorizontalString_noOffset1( )
  {
    int nParamsTot = 3;
    double params[3] = {21.0, 22.0, 100.0};   // X0, Y0, I_sky
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  retVal;
    vector<mp_par> parameterInfo;

    string  correctString = "21.00000000\t22.00000000\t100.0000000";
    string  outputString;

    // OK, now kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 0;
      currentParameterInfo.limited[1] = 0;
      currentParameterInfo.limits[0] = 0.0;
      currentParameterInfo.limits[1] = 0.0;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }

    modelObj2a->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj2a->AddParameterInfo(parameterInfo);

    int nParams = modelObj2a->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParamsTot);

    outputString = modelObj2a->PrintModelParamsHorizontalString(params);
    TS_ASSERT_EQUALS(outputString, correctString);
  }

  void testPrintModelParamsHorizontalString_noOffset2( )
  {
    int nParamsTot = 7;
    double params[7] = {21.0, 22.0, 0.0, 0.5, 50.0, 10.0, 100.0};   // X0, Y0, PA, ell, I_0, h, I_sky
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  retVal;
    vector<mp_par> parameterInfo;

    string  correctString = "21.00000000\t22.00000000\t0.000000000\t0.5000000000\t50.00000000\t10.00000000\t100.0000000";
    string  outputString;

    // OK, now kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 0;
      currentParameterInfo.limited[1] = 0;
      currentParameterInfo.limits[0] = 0.0;
      currentParameterInfo.limits[1] = 0.0;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }

    modelObj2b->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj2b->AddParameterInfo(parameterInfo);

    int nParams = modelObj2b->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParamsTot);

    outputString = modelObj2b->PrintModelParamsHorizontalString(params);
    TS_ASSERT_EQUALS(outputString, correctString);
  }

  void testPrintParamsToString_with_errors( void )
  {
    int nParamsTot = 7;
    double params[7] = {21.0, 22.0, 0.0, 0.5, 50.0, 10.0, 100.0};   // X0, Y0, PA, ell, I_0, h, I_sky
    double errors[7] = {0.1,  0.2 , 0.001, 0.002, 2.5, 1.09, 3.3};   // X0, Y0, PA, ell, I_0, h, I_sky
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  retVal;
    vector<string> outputVect;
  
    vector<string> correctStrings1, correctStrings2;
    correctStrings1.push_back("#\n");
    correctStrings1.push_back("#X0\t\t21.0000\n");
    correctStrings1.push_back("#Y0\t\t22.0000\n");
    correctStrings1.push_back("#FUNCTION Exponential\n");
    correctStrings1.push_back("#PA\t\t      0\n");
    correctStrings1.push_back("#ell\t\t    0.5\n");
    correctStrings1.push_back("#I_0\t\t     50\n");
    correctStrings1.push_back("#h\t\t     10\n");
    correctStrings1.push_back("#FUNCTION FlatSky\n");
    correctStrings1.push_back("#I_sky\t\t    100\n");
    correctStrings2.push_back("#\n");
    correctStrings2.push_back("#X0\t\t21.0000 # +/- 0.1000\n");
    correctStrings2.push_back("#Y0\t\t22.0000 # +/- 0.2000\n");
    correctStrings2.push_back("#FUNCTION Exponential\n");
    correctStrings2.push_back("#PA\t\t      0 # +/- 0.001\n");
    correctStrings2.push_back("#ell\t\t    0.5 # +/- 0.002\n");
    correctStrings2.push_back("#I_0\t\t     50 # +/- 2.5\n");
    correctStrings2.push_back("#h\t\t     10 # +/- 1.09\n");
    correctStrings2.push_back("#FUNCTION FlatSky\n");
    correctStrings2.push_back("#I_sky\t\t    100 # +/- 3.3000\n");

    string prefix = "#";
    vector<mp_par> parameterInfo_empty;
    vector<mp_par> parameterInfo;

    // kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 1;
      currentParameterInfo.limited[1] = 1;
      currentParameterInfo.limits[0] = 0.0 + i;
      currentParameterInfo.limits[1] = 200.0 + i;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }
    // specify that Y0 and Exponential h are fixed
    parameterInfo[1].fixed = 1;
    parameterInfo[5].fixed = 1;

    modelObj1->SetupModelImage(nSmallDataCols, nSmallDataRows);  
    modelObj1->AddParameterInfo(parameterInfo);

    // output without parameter limits
    retVal = modelObj1->PrintModelParamsToStrings(outputVect, params, NULL, prefix.c_str(), false);
    TS_ASSERT_EQUALS(retVal, 0);
    for (int i = 0; i < 9; i++) {
      TS_ASSERT_EQUALS(outputVect[i], correctStrings1[i]);
    }

    // output *with* erros
    outputVect.clear();
    retVal = modelObj1->PrintModelParamsToStrings(outputVect, params, errors, prefix.c_str(), false);
    TS_ASSERT_EQUALS(retVal, 0);
    for (int i = 0; i < 9; i++) {
      TS_ASSERT_EQUALS(outputVect[i], correctStrings2[i]);
    }
  }

  void testPrintModelParamsHorizontalString_withOffset( )
  {
    int nParamsTot = 7;
    double params[7] = {21.0, 22.0, 0.0, 0.5, 50.0, 10.0, 100.0};   // X0, Y0, PA, ell, I_0, h, I_sky
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  retVal;
    vector<mp_par> parameterInfo;

    string  correctString = "121.0000000\t222.0000000\t0.000000000\t0.5000000000\t50.00000000\t10.00000000\t100.0000000";
    string  outputString;

    // OK, now kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 0;
      currentParameterInfo.limited[1] = 0;
      currentParameterInfo.limits[0] = 0.0;
      currentParameterInfo.limits[1] = 0.0;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }

    modelObj2b->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj2b->AddParameterInfo(parameterInfo);
    modelObj2b->AddImageOffsets(100.0, 200.0);
    
    int nParams = modelObj2b->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParamsTot);

    outputString = modelObj2b->PrintModelParamsHorizontalString(params);
    TS_ASSERT_EQUALS(outputString, correctString);
  }

  void testPrintModelParamsHorizontalString_withOffset_2sets( )
  {
    int nParamsTot = 9;
    // X0, Y0, PA, ell, I_0, h; X0, Y0, I_sky
    double params[9] = {21.0, 22.0, 0.0, 0.5, 50.0, 10.0, 35.0, 45.0, 1.0};
    int  nDataVals = nSmallDataCols*nSmallDataRows;
    int  retVal;
    vector<mp_par> parameterInfo;

    string  outputString;
    string  correctString = "121.0000000\t222.0000000\t0.000000000\t0.5000000000\t";
    correctString += string("50.00000000\t10.00000000\t");
    correctString += string("135.0000000\t245.0000000\t1.000000000");

    // OK, now kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 0;
      currentParameterInfo.limited[1] = 0;
      currentParameterInfo.limits[0] = 0.0;
      currentParameterInfo.limits[1] = 0.0;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }

    modelObj3->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj3->AddParameterInfo(parameterInfo);
    modelObj3->AddImageOffsets(100.0, 200.0);
    
    int nParams = modelObj3->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParamsTot);

    outputString = modelObj3->PrintModelParamsHorizontalString(params);
    TS_ASSERT_EQUALS(outputString, correctString);    
  }

  void testGetImageOffsets_zeroOffset( )
  {
    int nParamsTot = 7;
    vector<mp_par> parameterInfo;
    double paramOffsets[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double paramOffsets_correct[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // OK, now kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 0;
      currentParameterInfo.limited[1] = 0;
      currentParameterInfo.limits[0] = 0.0;
      currentParameterInfo.limits[1] = 0.0;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }

    modelObj2b->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj2b->AddParameterInfo(parameterInfo);
    
    int nParams = modelObj2b->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParamsTot);

    modelObj2b->GetImageOffsets(paramOffsets);
    for (int i = 0; i < nParamsTot; i++)
      TS_ASSERT_EQUALS(paramOffsets[i], paramOffsets_correct[i]);
  }

  void testGetImageOffsets_withOffset( )
  {
    int nParamsTot = 7;
    vector<mp_par> parameterInfo;
    double paramOffsets[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double paramOffsets_correct[7] = {100.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // OK, now kit out parameterInfo
    mp_par currentParameterInfo;
    for (int i = 0; i < nParamsTot; i++) {
      currentParameterInfo.fixed = 0;
      currentParameterInfo.limited[0] = 0;
      currentParameterInfo.limited[1] = 0;
      currentParameterInfo.limits[0] = 0.0;
      currentParameterInfo.limits[1] = 0.0;
      currentParameterInfo.offset = 0.0;
      parameterInfo.push_back(currentParameterInfo);
    }

    modelObj2b->SetupModelImage(nSmallDataCols, nSmallDataRows);
    modelObj2b->AddParameterInfo(parameterInfo);
    modelObj2b->AddImageOffsets(100.0, 200.0);
    
    int nParams = modelObj2b->GetNParams();
    TS_ASSERT_EQUALS(nParams, nParamsTot);

    modelObj2b->GetImageOffsets(paramOffsets);
    for (int i = 0; i < nParamsTot; i++)
      TS_ASSERT_EQUALS(paramOffsets[i], paramOffsets_correct[i]);
  }

}; 



class TestPsfUsage : public CxxTest::TestSuite
{
public:
  vector<string>  functionList1, functionList2, functionList3, functionList3b;
  vector<string>  functionLabelList1, functionLabelList2, functionLabelList3, functionLabelList3b;
  vector<double>  parameterList1, parameterList2, parameterList3, parameterList3b;
  vector<mp_par>  paramLimits1, paramLimits2, paramLimits3, paramLimits3b;
  vector<int>  FunctionBlockIndices1, FunctionBlockIndices2, FunctionBlockIndices3, FunctionBlockIndices3b;
  bool  paramLimitsExist1, paramLimitsExist2, paramLimitsExist3, paramLimitsExist3b;
  mp_par  *parameterInfo;
  int  status;
  configOptions  userConfigOptions1, userConfigOptions2, userConfigOptions3, userConfigOptions3b;
  double  *smallDataImage;
  double  *smallErrorImage;
  double  *smallVarianceImage;
  double  *smallWeightImage;
  double  *smallMaskImage;
  int nSmallDataCols, nSmallDataRows;


  // Note that setUp() gets called prior to *each* individual test function!
  void setUp()
  {
    int  status;
    string  filename1 = CONFIG_FILE;
    string  filename3 = SIMPLE_CONFIG_FILE;
    string  filename3b = CONFIG_FILE3b;
    
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
  }

  void tearDown()
  {
    free(smallDataImage);
    free(smallErrorImage);
    free(smallVarianceImage);
    free(smallWeightImage);
    free(smallMaskImage);
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

    ModelObject *modelObj4a;
    // Initialize modelObj4a and add model function & params (FlatSky)
    modelObj4a = new ModelObject();
    status = AddFunctions(modelObj4a, functionList3, functionLabelList3, FunctionBlockIndices3, 
    						true, -1);
    ModelObject *modelObj4b;
    // Initialize modelObj4b and add model function & params (FlatSky)
    modelObj4b = new ModelObject();
    status = AddFunctions(modelObj4b, functionList3, functionLabelList3, FunctionBlockIndices3, 
    						true, -1);

    
    // This is the correct order
    // add PSF pixels first
    status4a = modelObj4a->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, smallPSFImage);
    // final setup for modelObj4a
    modelObj4a->SetupModelImage(nColumns, nRows);
    TS_ASSERT_EQUALS(status4a, 0);
    
    bool  psfPresent = modelObj4a->HasPSF();
    TS_ASSERT_EQUALS(psfPresent, true);

	// This is the WRONG order
    // final setup for modelObj4b *first*
    modelObj4b->SetupModelImage(nColumns, nRows);
    // then add PSF pixels = wrong order!
    status4b = modelObj4b->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, smallPSFImage);
    TS_ASSERT_EQUALS(status4b, -1);

    delete modelObj4a;
    delete modelObj4b;
  }

  // make sure ModelObject complains if we add a PSF with NaN pixel values
  void testCatchBadPSF( void )
  {
    int  nColumns = 10;
    int  nRows = 10;
    int  nColumns_psf = 3;
    int  nRows_psf = 3;
    int  nPixels_psf = 9;
    double  badPSFImage[9] = {sqrt(-1.0), 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
    int  status;

    ModelObject *modelObj5a;
    // Initialize modelObj5a and add model function & params (Exp + FlatSky)
    modelObj5a = new ModelObject();
    status = AddFunctions(modelObj5a, functionList1, functionLabelList1, FunctionBlockIndices1, 
    						true, -1);

    // This is the correct order
    // add PSF pixels first
    status = modelObj5a->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, badPSFImage);
    TS_ASSERT_EQUALS(status, -1);
    
    delete modelObj5a;
  }

  // make sure ModelObject complains if we add a PSF too small for bicubic interpolation
  // in a PointSource object (PSF must be 4x4 or larger)
  void testCatchBadPSF_forPointSource( void )
  {
    int  nColumns = 10;
    int  nRows = 10;
    int  nColumns_psf = 3;
    int  nRows_psf = 3;
    int  nPixels_psf = 9;
    double  tooSmallPsfImage[9] = {0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
    int  dummy, status;
    string  filename2 = CONFIG_FILE_POINTSOURCE;

    status = ReadConfigFile(filename2, true, functionList2, functionLabelList2, parameterList2, 
    						paramLimits2, FunctionBlockIndices2, paramLimitsExist2, userConfigOptions2);

    ModelObject *modelObj5a2;
    // Initialize modelObj5a2 and add model function & params (Exp + PointSource)
    modelObj5a2 = new ModelObject();
    status = AddFunctions(modelObj5a2, functionList2, functionLabelList2, FunctionBlockIndices2, 
    						true, -1);

    // First, see if we catch case of missing PSF
    status = modelObj5a2->SetupPsfInterpolation();
    TS_ASSERT_EQUALS(status, -1);

    // Second, see if we catch case of PSF that's present but too small
    status = modelObj5a2->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, tooSmallPsfImage);
    TS_ASSERT_EQUALS(status, 0);
    status = modelObj5a2->SetupPsfInterpolation();
    TS_ASSERT_EQUALS(status, -2);

    delete modelObj5a2;
  }

  // make sure ModelObject complains if we add oversampled PSF with NaN pixel values
//   void testCatchBadOversampledPSF( void )
//   {
//     int  nColumns = 10;
//     int  nRows = 10;
//     int  nColumns_psf = 3;
//     int  nRows_psf = 3;
//     int  nPixels_psf = 9;
//     double  goodPSFImage[9] = {0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
//     double  badPSFImage[9] = {sqrt(-1.0), 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
//     int  status;
// 
//     ModelObject *modelObj5b;
//     // Initialize modelObj5b and add model function & params (Exp + FlatSky)
//     modelObj5b = new ModelObject();
//     status = AddFunctions(modelObj5b, functionList1, FunctionBlockIndices1, true, -1);
//     
//     // This is the correct order
//     // add PSF pixels first
//     status = modelObj5b->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, goodPSFImage);
//     TS_ASSERT_EQUALS(status, 0);
//     // final setup for modelObj5b
//     modelObj5b->SetupModelImage(nColumns, nRows);
// 
//     bool  overPsfPresent = modelObj5b->HasOversampledPSF();
//     TS_ASSERT_EQUALS(overPsfPresent, false);
//     
//     // now add the bad PSF as an oversampled PSF
//     status = modelObj5b->AddOversampledPSFVector(nPixels_psf, nColumns_psf, nRows_psf, 
//     					badPSFImage, 1, 1,2, 1,2);
//     TS_ASSERT_EQUALS(status, -1);
// 
//     delete modelObj5b;
//   }
  
  
  void testOversampledPSF_newMethod( void )
  {
    int  nColumns = 10;
    int  nRows = 10;
    int  nColumns_psf = 3;
    int  nRows_psf = 3;
    int  nPixels_psf = 9;
    double  goodPSFImage[9] = {0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
    double  overPSFImage_ref[9] = {0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 0.0};
    double *  overPSFImage_malloc;
    PsfOversamplingInfo *osampleInfo_ptr;
    string  oversampleRegion = "1:2,1:2";
    int  status;

    ModelObject *modelObj5d;
    // Initialize modelObj5d and add model function & params (Exp + FlatSky)
    modelObj5d = new ModelObject();
    status = AddFunctions(modelObj5d, functionList1, functionLabelList1, FunctionBlockIndices1, 
    						true, -1);

    // allocate memory from heap for overPSFImage_malloc, since PsfOversamplingInfo will
    // try to free it later.
    overPSFImage_malloc = (double *)malloc(9*sizeof(double));
    for (int i = 0; i < 9; i++)
      overPSFImage_malloc[i] = overPSFImage_ref[i];
  
    // This is the correct order
    // add PSF pixels first
    status = modelObj5d->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, goodPSFImage);
    TS_ASSERT_EQUALS(status, 0);
    // final setup for modelObj5d
    modelObj5d->SetupModelImage(nColumns, nRows);

    osampleInfo_ptr = new PsfOversamplingInfo(overPSFImage_malloc, nColumns_psf, nRows_psf, 1,
    										oversampleRegion);
    
    status = modelObj5d->AddOversampledPsfInfo(osampleInfo_ptr);
    TS_ASSERT_EQUALS(status, 0);

    int  overPsfPresent = modelObj5d->HasOversampledPSF();
    TS_ASSERT_EQUALS(overPsfPresent, true);
    
    delete modelObj5d;
    delete osampleInfo_ptr;
  }
};
  
