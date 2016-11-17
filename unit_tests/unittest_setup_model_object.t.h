// Unit tests for code in setup_model_object.cpp

// Evidence for the idea that we're now doing integration tests is the large
// number of other modules that have to be compiled along with model_object.cpp
// to get these tests to work...

// See run_unittest_setup_model_object.sh for how to compile & run this


// Things to test:


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
  ModelObject *modelObj4c;
  ModelObject *modelObj5a;
  ModelObject *modelObj5b;
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


  // Note that setUp() gets called prior to *each* individual test function!
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

    smallMaskImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallMaskImage[0] = 0.0;
    smallMaskImage[1] = 1.0;
    smallMaskImage[2] = 0.0;
    smallMaskImage[3] = 0.0;


    modelObj1 = new ModelObject();    
    // Initialize modelObj1; set up internal FunctionObjects vector (Exp + FlatSky) inside
    status = ReadConfigFile(filename1, true, functionList1, parameterList1, 
  								paramLimits1, FunctionBlockIndices1, paramLimitsExist1, userConfigOptions1);

    
  }

  void tearDown()
  {
    free(smallDataImage);
    free(smallErrorImage);
    free(smallWeightImage);
    free(smallMaskImage);
  }
  
  
  void testSetupMakeimage_simple( void )
  {
  ModelObject *theModel = nullptr;
  OptionsBase *optionsPtr;
  vector<int> nColumnsRowsVect;
  
  optionsPtr = new MakeimageOptions();
  nColumnsRowsVect.push_back(2);
  nColumnsRowsVect.push_back(2);
  
  theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL);
  
  long nDataVals_true = 4;
  long nDataVals = theModel->GetNDataValues();
  TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
  }
  
};
