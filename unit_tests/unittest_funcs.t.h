// Unit tests for image function classes
//
// "FUNCTION-SPECIFIC" denotes code which needs to be changed when testing
// a different function object (e.g., Sersic vs Exponential, etc.).

// See run_unittest_funcs.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>
#include <vector>
using namespace std;

// test stuff (not official image functions)
#include "function_objects/function_object.h"
#include "function_objects/func_gauss_extraparams.h"
// official image functions
#include "function_objects_1d/func1d_exp_test.h"
#include "function_objects/func_flatsky.h"
#include "function_objects/func_exp.h"
#include "function_objects/func_gaussian.h"
#include "function_objects/func_moffat.h"
#include "function_objects/func_sersic.h"
#include "function_objects/func_king.h"
#include "function_objects/func_king2.h"
#include "function_objects/func_pointsource.h"
#include "function_objects/func_broken-exp.h"
#include "function_objects/func_broken-exp2d.h"
#include "function_objects/func_edge-on-disk.h"
#include "function_objects/func_ferrersbar3d.h"
#include "function_objects/func_double-broken-exp.h"
#include "function_objects/func_nuker.h"
//#include "function_objects/func_spline-profile.h"

const double  DELTA = 1.0e-9;
const double  DELTA_e9 = 1.0e-9;
const double  DELTA_e7 = 1.0e-7;

const double PI = 3.14159265358979;


// Testing temporary 1D function (exponential with linear input and output,
// plus SetExtraParams testing)
class TestExp1DTest : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc1a;
  FunctionObject  *thisFunc1b;
  FunctionObject  *thisFunc2;
  FunctionObject  *thisFunc3;

public:
  void setUp()
  {
    bool  subsampleFlag = false;
    thisFunc1a = new Exponential1D_test();  // general testing (no numerical output)
    thisFunc1a->SetSubsampling(subsampleFlag);
    thisFunc1b = new Exponential1D_test();  // general testing (case without good extra params)
    thisFunc1b->SetSubsampling(subsampleFlag);
    thisFunc2 = new Exponential1D_test();
    thisFunc2->SetSubsampling(subsampleFlag);
    thisFunc3 = new Exponential1D_test();   // for testing numerical output
    thisFunc3->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc1a;
    delete thisFunc1b;
    delete thisFunc2;
    delete thisFunc3;
  }

  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 2;
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("h");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc1a->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc1a->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
    // check for (non-existent) parameter units
    TS_ASSERT_EQUALS(thisFunc1a->HasParameterUnits(), false)
    vector<string>  resultStrings;
    string  correctStr = "";
    thisFunc1a->GetParameterUnits(resultStrings);
    for (int i = 0; i < correctNParams; i++)
      TS_ASSERT_EQUALS(resultStrings[i], correctStr);
  }
  
  void testHasExtraParams( void )
  {
    bool  returnVal = thisFunc1a->HasExtraParams();
    TS_ASSERT_EQUALS( returnVal, true );
  }
  
  void testSetExtraParams1_GoodNameAndValue( void )
  {
    map<string, string> theMap;
    string  keyword = "floor";
    string  value = "0";
    theMap[keyword] = value;
    
    int  status = thisFunc1a->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 1 );

    bool  returnVal = thisFunc1a->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, true );
  }

  void testSetExtraParams1_EmptyMap( void )
  {
    map<string, string> theMap;

    int  status = thisFunc1b->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, -1 );

    bool  returnVal = thisFunc1b->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );
  }

  void testSetExtraParams1_BadName( void )
  {
    map<string, string> theMap;
    string  keyword = "interp";
    string  value = "100";
    theMap[keyword] = value;
    
    int  status = thisFunc1b->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 0 );

    bool  returnVal = thisFunc1b->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );
  }

  // non-numeric value for parameter expecting a number
  void testSetExtraParams1_BadValue( void )
  {
    map<string, string> theMap;
    string  keyword = "floor";
    string  value = "bob";
    theMap[keyword] = value;
    
    int  status = thisFunc1b->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, -3 );

    bool  returnVal = thisFunc1b->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );
  }

  // testing basic numerical output without ExtraParams complications
  void testCalculations( void )
  {
    // centered at x0 = 0
    double  x0 = 0.0;
    // FUNCTION-SPECIFIC:
    // test setup: exponential with I_0 = 1, h = 10,
    double  params[2] = {1.0, 10.0};
    
    thisFunc2->Setup(params, 0, x0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(0.0), 1.0, DELTA );
    // r = 1 value
    double  rEqualsOneValue = 1.0*exp(-1.0/10.0);
    TS_ASSERT_DELTA( thisFunc2->GetValue(1.0), rEqualsOneValue, DELTA );
  }

  // testing numerical output when we set ExtraParams
  void testCalculations_with_extra_params( void )
  {
    // centered at x0 = 0
    int  status;
    double  x0 = 0.0;
    // FUNCTION-SPECIFIC:
    // test setup: exponential with I_0 = 1, h = 10,
    double  params[2] = {1.0, 10.0};
    map<string, string> theMap;
    string  keyword = "floor";

    double  rEqualsOneValue = 1.0*exp(-1.0/10.0);
    
    thisFunc3->Setup(params, 0, x0);
    
    // Pre-ExtraParams:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc3->GetValue(0.0), 1.0, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc3->GetValue(1.0), rEqualsOneValue, DELTA );

    // Set ExtraParams = 0
    theMap[keyword] = string("0");
    status = thisFunc3->SetExtraParams(theMap);
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc3->GetValue(0.0), 1.0, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc3->GetValue(1.0), rEqualsOneValue, DELTA );

    // Set ExtraParams = 100
    theMap[keyword] = string("100");
    status = thisFunc3->SetExtraParams(theMap);
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc3->GetValue(0.0), 101.0, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc3->GetValue(1.0), 100.0 + rEqualsOneValue, DELTA );
  }
};



class TestExponential : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new Exponential();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 4;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("h");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames ); 
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular exponential with I_0 = 1, h = 10,
    double  params[4] = {90.0, 0.0, 1.0, 10.0};
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    double  rEqualsOneValue = 1.0*exp(-1.0/10.0);
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 11.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );

  }

  void testLabels( void )
  {
    string  result;
    
    string  correct = "";
    result = thisFunc->GetLabel();
    TS_ASSERT_EQUALS(result, correct);
    
    string  input = "this is a label";
    thisFunc->SetLabel(input);
    result = thisFunc->GetLabel();
    TS_ASSERT_EQUALS(result, input);
  }
  
  void testUnitNames( void )
  {
    int  nParams = 4;
    vector<string>  resultStrings;
    
    TS_ASSERT_EQUALS(thisFunc->HasParameterUnits(), true)

    vector<string>  correctStrings = {"deg (CCW from +y axis)", "", "counts/pixel", "pixels"};
    thisFunc->GetParameterUnits(resultStrings);
    for (int i = 0; i < nParams; i++)
      TS_ASSERT_EQUALS(resultStrings[i], correctStrings[i]);
  }
  
  void testIsBackground( void )
  {
    bool result = thisFunc->IsBackground();
    TS_ASSERT_EQUALS(result, false);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, true);
  }

  void testTotalFlux_calcs( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular or elliptical exponential with I_0 = 1, h = 10,
    double  I0 = 1.0;
    double  h = 10.0;
    double  ell1 = 0.0;
    double  ell2 = 0.25;
    double  params0[4] = {90.0, ell2, 0.0, h};
    double  params1[4] = {90.0, ell1, I0, h};
    double  params2[4] = {90.0, ell2, I0, h};


    // FUNCTION-SPECIFIC:
    thisFunc->Setup(params0, 0, x0, y0);
    double  correctCircFlux = 0.0;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA );
    
    thisFunc->Setup(params1, 0, x0, y0);
    correctCircFlux = 2.0*PI*I0*h*h;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA );
    
    thisFunc->Setup(params2, 0, x0, y0);
    double  correctEllFlux = (1 - ell2) * correctCircFlux;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctEllFlux, DELTA );
  }
};



class TestSersic : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new Sersic();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 5;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("n");
    correctParamNames.push_back("I_e");
    correctParamNames.push_back("r_e");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular Sersic with n = 2, I_e = 1, r_e = 10,
    double  params[5] = {90.0, 0.0, 2.0, 1.0, 10.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 39.333062332325284, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), 12.315472433581958, DELTA );
    // r = r_e value
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), 1.0, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), 1.0, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), 1.0, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), 1.0, DELTA );

  }

  void testUnitNames( void )
  {
    int  nParams = 5;
    vector<string>  resultStrings;
    
    TS_ASSERT_EQUALS(thisFunc->HasParameterUnits(), true)

    vector<string>  correctStrings = {"deg (CCW from +y axis)", "", "", "counts/pixel", "pixels"};
    thisFunc->GetParameterUnits(resultStrings);
    for (int i = 0; i < nParams; i++)
      TS_ASSERT_EQUALS(resultStrings[i], correctStrings[i]);
  }
  
  void testIsBackground( void )
  {
    bool result = thisFunc->IsBackground();
    TS_ASSERT_EQUALS(result, false);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, true);
  }

  void testTotalFlux_calcs( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular or elliptical Sersic with 
    double  I_e = 1.0;
    double  r_e = 10.0;
    double  ell1 = 0.0;
    double  ell2 = 0.25;
    double  params0[5] = {90.0, ell1, 1.0, 0.0, r_e};
    double  params1[5] = {90.0, ell1, 1.0, I_e, r_e};
    double  params2[5] = {90.0, ell2, 1.0, I_e, r_e};
    double  params3[5] = {90.0, ell1, 4.0, I_e, r_e};
    double  params4[5] = {90.0, ell2, 4.0, I_e, r_e};
    double  params5[5] = {90.0, ell1, 0.5, I_e, r_e};
    double  params6[5] = {90.0, ell1, 10.0, I_e, r_e};


    // FUNCTION-SPECIFIC:
    thisFunc->Setup(params0, 0, x0, y0);
    double  correctCircFlux = 0.0;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA );
    
    thisFunc->Setup(params1, 0, x0, y0);
    correctCircFlux = 1194.839934018338;   // from astro_utils.LSersic, using non-exact b_n
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA_e9 );
    
    thisFunc->Setup(params2, 0, x0, y0);
    double  correctEllFlux = (1 - ell2) * correctCircFlux;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctEllFlux, DELTA_e9 );

    thisFunc->Setup(params3, 0, x0, y0);
    correctCircFlux = 2266.523319362499;   // from astro_utils.LSersic, using non-exact b_n
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA_e9 );
    
    thisFunc->Setup(params4, 0, x0, y0);
    correctEllFlux = (1 - ell2) * correctCircFlux;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctEllFlux, DELTA_e9 );

    thisFunc->Setup(params5, 0, x0, y0);
    correctCircFlux = 906.370725112331;   // from astro_utils.LSersic, using non-exact b_n
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA_e9 );

    // n = 10 case diverges slightly more from Python estimate (but not enough
    // to worry about)
    thisFunc->Setup(params6, 0, x0, y0);
    correctCircFlux = 3546.3105962151512;   // from astro_utils.LSersic, using non-exact b_n
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA_e7 );
  }
};



class TestGaussian : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new Gaussian();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 4;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("sigma");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular Gaussian with I_e = 1, sigma = 10,
    double  params[4] = {90.0, 0.0, 1.0, 10.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    double  rEqualsOneValue = 1.0*exp(-1.0/(2*10.0*10.0));
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = sigma value
    double  rEqualsSigmaValue = 1.0*exp(-(10.0*10.0)/(2*10.0*10.0));
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );

  }
  
  void testUnitNames( void )
  {
    int  nParams = 4;
    vector<string>  resultStrings;
    
    TS_ASSERT_EQUALS(thisFunc->HasParameterUnits(), true)

    vector<string>  correctStrings = {"deg (CCW from +y axis)", "", "counts/pixel", "pixels"};
    thisFunc->GetParameterUnits(resultStrings);
    for (int i = 0; i < nParams; i++)
      TS_ASSERT_EQUALS(resultStrings[i], correctStrings[i]);
  }
  
  void testIsBackground( void )
  {
    bool result = thisFunc->IsBackground();
    TS_ASSERT_EQUALS(result, false);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, true);
  }

  void testTotalFlux_calcs( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular or elliptical Gaussian with I_0 = 1, sigma = 10,
    double  I0 = 1.0;
    double  sigma = 10.0;
    double  ell1 = 0.0;
    double  ell2 = 0.25;
    double  params1[4] = {90.0, ell1, I0, sigma};
    double  params2[4] = {90.0, ell2, I0, sigma};


    // FUNCTION-SPECIFIC:
    thisFunc->Setup(params1, 0, x0, y0);
    double  correctCircFlux = 2.0*PI*I0*sigma*sigma;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA );
    
    thisFunc->Setup(params2, 0, x0, y0);
    double  correctEllFlux = (1 - ell2) * correctCircFlux;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctEllFlux, DELTA );
  }
};


class TestGaussianExtraParams : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  FunctionObject  *thisFunc2, *thisFunc3;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new GaussianExtraParams();
    thisFunc->SetSubsampling(subsampleFlag);
    thisFunc2 = new GaussianExtraParams();
    thisFunc2->SetSubsampling(subsampleFlag);
    thisFunc3 = new GaussianExtraParams();
    thisFunc3->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
    delete thisFunc2;
    delete thisFunc3;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 4;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("sigma");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular Gaussian with I_e = 1, sigma = 10,
    double  params[4] = {90.0, 0.0, 1.0, 10.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    double  rEqualsOneValue = 1.0*exp(-1.0/(2*10.0*10.0));
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = sigma value
    double  rEqualsSigmaValue = 1.0*exp(-(10.0*10.0)/(2*10.0*10.0));
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );

  }
  
  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, true);
  }

  void testTotalFlux_calcs( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular or elliptical Gaussian with I_0 = 1, sigma = 10,
    double  I0 = 1.0;
    double  sigma = 10.0;
    double  ell1 = 0.0;
    double  ell2 = 0.25;
    double  params1[4] = {90.0, ell1, I0, sigma};
    double  params2[4] = {90.0, ell2, I0, sigma};


    // FUNCTION-SPECIFIC:
    thisFunc->Setup(params1, 0, x0, y0);
    double  correctCircFlux = 2.0*PI*I0*sigma*sigma;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctCircFlux, DELTA );
    
    thisFunc->Setup(params2, 0, x0, y0);
    double  correctEllFlux = (1 - ell2) * correctCircFlux;
    TS_ASSERT_DELTA( thisFunc->TotalFlux(), correctEllFlux, DELTA );
  }
  

  // ExtraParams tests
  void testHasExtraParams( void )
  {
    bool  returnVal = thisFunc->HasExtraParams();
    TS_ASSERT_EQUALS( returnVal, true );
  }
  
  void testSetExtraParams1_GoodNameAndValue( void )
  {
    map<string, string> theMap;
    string  keyword = "floor";
    string  value = "0";
    theMap[keyword] = value;
    
    int  status = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 1 );

    bool  returnVal = thisFunc->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, true );
  }

  void testSetExtraParams1_EmptyMap( void )
  {
    map<string, string> theMap;

    int  status = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, -1 );
  }

  void testSetExtraParams1_BadName( void )
  {
    map<string, string> theMap;
    string  keyword = "interp";
    string  value = "100";
    theMap[keyword] = value;
    
    int  status = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 0 );
  }

  // non-numeric value for parameter expecting a number
  void testSetExtraParams1_BadValue( void )
  {
    map<string, string> theMap;
    string  keyword = "floor";
    string  value = "bob";
    theMap[keyword] = value;
    
    int  status = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, -3 );
  }
  // testing numerical output when we set ExtraParams
  void testCalculations_with_extra_params( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular Gaussian with I_e = 1, sigma = 10,
    double  params[4] = {90.0, 0.0, 1.0, 10.0};
    map<string, string> theMap;
    string  keyword = "floor";
    int  status;

    double  rEqualsOneValue = 1.0*exp(-1.0/(2*10.0*10.0));
    double  rEqualsSigmaValue = 1.0*exp(-(10.0*10.0)/(2*10.0*10.0));
    
    thisFunc2->Setup(params, 0, x0, y0);
    
    // Pre-ExtraParams:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc2->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    // r = sigma value
    TS_ASSERT_DELTA( thisFunc2->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc2->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );

    // Set ExtraParams = 0
    theMap[keyword] = string("0");
    status = thisFunc2->SetExtraParams(theMap);
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc2->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    // r = sigma value
    TS_ASSERT_DELTA( thisFunc2->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc2->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );

    // Set ExtraParams = 100
    theMap[keyword] = string("100");
    status = thisFunc2->SetExtraParams(theMap);
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(10.0, 10.0), 101.0, DELTA );
    // r = 1 value
    TS_ASSERT_DELTA( thisFunc2->GetValue(11.0, 10.0), 100 + rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc2->GetValue(9.0, 10.0), 100 + rEqualsOneValue, DELTA );
    // r = sigma value
    TS_ASSERT_DELTA( thisFunc2->GetValue(20.0, 10.0), 100 + rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc2->GetValue(10.0, 20.0), 100 + rEqualsSigmaValue, DELTA );
  }
};


class TestFlatSky : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new FlatSky();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 1;
    correctParamNames.push_back("I_sky");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: FlatSky with I_0 = 1.0
    double  params1[1] = {1.0};
    double  params2[1] = {-1.0};
    
    
    thisFunc->Setup(params1, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value (should be identical)
    double  rEqualsOneValue = 1.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    
    // Same, but with negative value
    thisFunc->Setup(params2, 0, x0, y0);
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), -1.0, DELTA );
    // r = 1 value (should be identical)
    rEqualsOneValue = -1.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
  }

  void testUnitNames( void )
  {
    int  nParams = 1;
    vector<string>  resultStrings;
    
    TS_ASSERT_EQUALS(thisFunc->HasParameterUnits(), true)

    vector<string>  correctStrings = {"counts/pixel"};
    thisFunc->GetParameterUnits(resultStrings);
    for (int i = 0; i < nParams; i++)
      TS_ASSERT_EQUALS(resultStrings[i], correctStrings[i]);
  }
  
  void testIsBackground( void )
  {
    bool result = thisFunc->IsBackground();
    TS_ASSERT_EQUALS(result, true);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestMoffat : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new Moffat();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 5;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("fwhm");
    correctParamNames.push_back("beta");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular Moffat with I_e = 1, fwhm = 10, beta = 3.0
    double  params1[5] = {90.0, 0.0, 1.0, 10.0, 3.0};
    double  params2[5] = {90.0, 0.0, 1.0, 10.0, 1.0};
    double  params3[5] = {90.0, 0.5, 1.0, 10.0, 1.0};
    double  rEqualsOneValue, rEqualsFWHMValue;
    
    
    thisFunc->Setup(params1, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    rEqualsOneValue = 0.96944697430705418;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = fwhm value
    rEqualsFWHMValue = 0.11784501202946467;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsFWHMValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsFWHMValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsFWHMValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsFWHMValue, DELTA );

    // Same, but with beta = 1.0
    thisFunc->Setup(params2, 0, x0, y0);
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // r = 1 value
    rEqualsOneValue = 0.96153846153846145;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = fwhm value
    rEqualsFWHMValue = 0.20000000000000001;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsFWHMValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsFWHMValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsFWHMValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsFWHMValue, DELTA );

    // Same, but with ellipticity = 0.5
    thisFunc->Setup(params3, 0, x0, y0);
    // a = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 1.0, DELTA );
    // a = 1 value
    rEqualsOneValue = 0.96153846153846145;
    // along major axis
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    // along minor axis: delta-y = +/-0.5
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.5), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.5), rEqualsOneValue, DELTA );
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestModifiedKing : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new ModifiedKing();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 6;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("r_c");
    correctParamNames.push_back("r_t");
    correctParamNames.push_back("alpha");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations_alpha2( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular ModifiedKing with I_0 = 100, r_c = 5, r_t = 9, alpha = 2
    double  params[6] = {90.0, 0.0, 100.0, 5.0, 9.0, 2.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 100.0, DELTA );
    // r = 1 value (calculated with Python function astro_funcs.ModifiedKing_intensity)
    double  rEqualsOneValue = 92.59162886983584;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = 5 value (calculated with Python function astro_funcs.ModifiedKing_intensity)
    double  rEqualsFiveValue = 18.53857147439379;
    TS_ASSERT_DELTA( thisFunc->GetValue(15.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(5.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 5.0), rEqualsFiveValue, DELTA );
    // r = 10 value
    double  rEqualsSigmaValue = 0.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );
  }

  void testCalculations_alpha1( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular ModifiedKing with I_0 = 100, r_c = 5, r_t = 9, alpha = 1
    double  params[6] = {90.0, 0.0, 100.0, 5.0, 9.0, 1.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 100.0, DELTA );
    // r = 1 value (calculated with Python function astro_funcs.ModifiedKing_intensity)
    double  rEqualsOneValue = 94.96676163342828;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = 5 value (calculated with Python function astro_funcs.ModifiedKing_intensity)
    double  rEqualsFiveValue = 34.567901234567906;
    TS_ASSERT_DELTA( thisFunc->GetValue(15.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(5.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 5.0), rEqualsFiveValue, DELTA );
    // r = 10 value
    double  rEqualsSigmaValue = 0.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );
  }

  void testCalculations_alpha3( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular ModifiedKing with I_0 = 100, r_c = 5, r_t = 9, alpha = 3
    double  params[6] = {90.0, 0.0, 100.0, 5.0, 9.0, 3.0};
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 100.0, DELTA );
    // r = 1 value (calculated with Python function astro_funcs.ModifiedKing_intensity)
    double  rEqualsOneValue = 90.14642898820081;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = 5 value (calculated with Python function astro_funcs.ModifiedKing_intensity)
    double  rEqualsFiveValue = 9.744462560365605;
    TS_ASSERT_DELTA( thisFunc->GetValue(15.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(5.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 5.0), rEqualsFiveValue, DELTA );
    // r = 10 value
    double  rEqualsSigmaValue = 0.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestModifiedKing2 : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new ModifiedKing2();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 6;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("r_c");
    correctParamNames.push_back("c");
    correctParamNames.push_back("alpha");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations_alpha2( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular ModifiedKing2 with I_0 = 100, r_c = 5, c = 9/5, alpha = 2
    // (same as ModifiedKing with I_0 = 100, r_c = 5, r_t = 9, alpha = 2)
    double  params[6] = {90.0, 0.0, 100.0, 5.0, 9.0/5.0, 2.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 100.0, DELTA );
    // r = 1 value (calculated with Python function astro_funcs.ModifiedKing2_intensity)
    double  rEqualsOneValue = 92.59162886983584;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = 5 value (calculated with Python function astro_funcs.ModifiedKing2_intensity)
    double  rEqualsFiveValue = 18.53857147439379;
    TS_ASSERT_DELTA( thisFunc->GetValue(15.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(5.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 5.0), rEqualsFiveValue, DELTA );
    // r = 10 value
    double  rEqualsSigmaValue = 0.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );
  }

  void testCalculations_alpha1( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular ModifiedKing2 with I_0 = 100, r_c = 5, c = 0/5, alpha = 1
    // (same as ModifiedKing with I_0 = 100, r_c = 5, r_t = 9, alpha = 1)
    double  params[6] = {90.0, 0.0, 100.0, 5.0, 9.0/5.0, 1.0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 100.0, DELTA );
    // r = 1 value (calculated with Python function astro_funcs.ModifiedKing2_intensity)
    double  rEqualsOneValue = 94.96676163342828;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = 5 value (calculated with Python function astro_funcs.ModifiedKing2_intensity)
    double  rEqualsFiveValue = 34.567901234567906;
    TS_ASSERT_DELTA( thisFunc->GetValue(15.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(5.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 5.0), rEqualsFiveValue, DELTA );
    // r = 10 value
    double  rEqualsSigmaValue = 0.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );
  }

  void testCalculations_alpha3( void )
  {
    // centered at x0,y0 = 10,10
    double  x0 = 10.0;
    double  y0 = 10.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular ModifiedKing2 with I_0 = 100, r_c = 5, c = 9/5, alpha = 3
    // (same as ModifiedKing with I_0 = 100, r_c = 5, r_t = 9, alpha = 3)
    double  params[6] = {90.0, 0.0, 100.0, 5.0, 9.0/5.0, 3.0};
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 10.0), 100.0, DELTA );
    // r = 1 value (calculated with Python function astro_funcs.ModifiedKing2_intensity)
    double  rEqualsOneValue = 90.14642898820081;
    TS_ASSERT_DELTA( thisFunc->GetValue(11.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(9.0, 10.0), rEqualsOneValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 9.0), rEqualsOneValue, DELTA );
    // r = 5 value (calculated with Python function astro_funcs.ModifiedKing2_intensity)
    double  rEqualsFiveValue = 9.744462560365605;
    TS_ASSERT_DELTA( thisFunc->GetValue(15.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(5.0, 10.0), rEqualsFiveValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 5.0), rEqualsFiveValue, DELTA );
    // r = 10 value
    double  rEqualsSigmaValue = 0.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(20.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 20.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 10.0), rEqualsSigmaValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(10.0, 0.0), rEqualsSigmaValue, DELTA );
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestPointSource : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc2;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new PointSource();
    thisFunc->SetSubsampling(subsampleFlag);
    thisFunc2 = new PointSource();
    thisFunc2->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
    delete thisFunc2;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 1;
    correctParamNames.push_back("I_tot");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
  }

  void testBasic_PointSourceRelated( void )
  {
    bool  answer;

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->IsPointSource(), true );
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, true);
  }

  void testHasExtraParams( void )
  {
    bool  returnVal = thisFunc->HasExtraParams();
    TS_ASSERT_EQUALS( returnVal, true );
    returnVal = thisFunc->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );
  }
  
  void testSetExtraParams1_GoodNameAndValues( void )
  {
    map<string, string> theMap;
    string  keyword = "method";
    string  value_bicubic = "bicubic";
    string  value_lanczos2 = "lanczos2";
    vector<string>  extraParamLines;
    
    theMap[keyword] = value_bicubic;
    int  status = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 1 );

    bool  returnVal = thisFunc->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, true );
    
    string  returnVal2 = thisFunc->GetInterpolationType();
    TS_ASSERT_EQUALS( returnVal2, value_bicubic );
    
    thisFunc->GetExtraParamsDescription(extraParamLines);
    TS_ASSERT_EQUALS( extraParamLines.size(), 3 );
    if (extraParamLines.size() == 3) {
      string  correctString1 = "OPTIONAL_PARAMS_START";
      string  correctString2 = "method   bicubic";
      string  correctString3 = "OPTIONAL_PARAMS_END";
      TS_ASSERT_EQUALS( extraParamLines[0], correctString1 );
      TS_ASSERT_EQUALS( extraParamLines[1], correctString2 );
      TS_ASSERT_EQUALS( extraParamLines[2], correctString3 );
    }


    theMap[keyword] = value_lanczos2;
    status = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 1 );

    returnVal = thisFunc->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, true );
    
    returnVal2 = thisFunc->GetInterpolationType();
    TS_ASSERT_EQUALS( returnVal2, value_lanczos2 );

    
    thisFunc->GetExtraParamsDescription(extraParamLines);
    TS_ASSERT_EQUALS( extraParamLines.size(), 3 );
    if (extraParamLines.size() == 3) {
      string  correctString1 = "OPTIONAL_PARAMS_START";
      string  correctString2b = "method   lanczos2";
      string  correctString3 = "OPTIONAL_PARAMS_END";
      TS_ASSERT_EQUALS( extraParamLines[0], correctString1 );
      TS_ASSERT_EQUALS( extraParamLines[1], correctString2b );
      TS_ASSERT_EQUALS( extraParamLines[2], correctString3 );
    }
  }

  void testSetExtraParams2_EmptyMap( void )
  {
    map<string, string> theMap;

    int  status = thisFunc2->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, -1 );

    bool  returnVal = thisFunc2->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );

    bool  returnVal2 = thisFunc2->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal2, false );
  }

  void testSetExtraParams2_BadName( void )
  {
    map<string, string> theMap;
    string  keyword = "interp";
    string  value = "100";
    theMap[keyword] = value;
    
    int  status = thisFunc2->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, 0 );

    bool  returnVal = thisFunc2->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );

    bool  returnVal2 = thisFunc2->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal2, false );
  }

  // non-numeric value for parameter expecting a number
  void testSetExtraParams2_BadValue( void )
  {
    map<string, string> theMap;
    string  keyword = "method";
    string  value = "bob";
    theMap[keyword] = value;
    
    int  status = thisFunc2->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( status, -3 );

    bool  returnVal = thisFunc2->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal, false );

    bool  returnVal2 = thisFunc2->ExtraParamsSet();
    TS_ASSERT_EQUALS( returnVal2, false );
  }
};


class TestNukerLaw : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new NukerLaw();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 7;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_b");
    correctParamNames.push_back("r_b");
    correctParamNames.push_back("alpha");
    correctParamNames.push_back("beta");
    correctParamNames.push_back("gamma");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }

  void testCalculations_circular( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular nuker-law with I_b = 100, r_b = 10
    double  params[7] = {90.0, 0.0, 100.0, 10.0, 10.0, 1.5, 0.5};
    double  rEqualsOneValue, rEqualsRbValue, rEquals2RbValue;
    
    // test
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 1 value
    rEqualsOneValue = 338.924527730594;
    TS_ASSERT_DELTA( thisFunc->GetValue(101.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(99.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 99.0), rEqualsOneValue, DELTA);
    // r = r_break value
    double  r_b = 10.0;
    rEqualsRbValue = 100.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 + r_b, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 - r_b, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 - r_b), rEqualsRbValue, DELTA);
    // r = 2*r_break value
    rEquals2RbValue = 37.889215669066;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 + 2*r_b, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 - 2*r_b, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 - 2*r_b), rEquals2RbValue, DELTA);
  }

  void testCalculations_elliptical( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: ell=0.5 broken-exp with nuker-law with I_b = 100, r_b = 10
    double  params[7] = {90.0, 0.5, 100.0, 10.0, 10.0, 1.5, 0.5};
    double  rEqualsOneValue, rEqualsRbValue, rEquals2RbValue;
    
    // test
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 1 value; account for ellipticity = 0.5 for y offsets
    rEqualsOneValue = 338.924527730594;
    TS_ASSERT_DELTA( thisFunc->GetValue(101.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(99.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 99.5), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.5), rEqualsOneValue, DELTA);
    // r = r_break value
    double  r_b = 10.0;
    rEqualsRbValue = 100.0;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 + r_b, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 - r_b, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 - 0.5*r_b), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 + 0.5*r_b), rEqualsRbValue, DELTA);
    // r = 2*r_break value
    rEquals2RbValue = 37.889215669066;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 + 2*r_b, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 - 2*r_b, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 - r_b), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 + r_b), rEquals2RbValue, DELTA);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestBrokenExponential : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new BrokenExponential();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 7;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("h1");
    correctParamNames.push_back("h2");
    correctParamNames.push_back("r_break");
    correctParamNames.push_back("alpha");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }

  void testCalculations_circular( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular broken-exp with I_0 = 100, h1 = 50, h2 = 10, r_break = 50
    double  params[7] = {90.0, 0.0, 100.0, 50.0, 10.0, 50.0, 10.0};
    double  rEqualsOneValue, rEqualsRbValue, rEquals2RbValue;
    
    // test
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), 100.0, DELTA );
    // r = 1 value
    rEqualsOneValue = 98.019867330675524;
    TS_ASSERT_DELTA( thisFunc->GetValue(101.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(99.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 99.0), rEqualsOneValue, DELTA);
    // r = r_break value
    rEqualsRbValue = 36.584512991317212;
    TS_ASSERT_DELTA( thisFunc->GetValue(150.0, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(50.0, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 50.0), rEqualsRbValue, DELTA);
    // r = 2*r_break value
    rEquals2RbValue = 0.24787521766663584;
    TS_ASSERT_DELTA( thisFunc->GetValue(200.0, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 0.0), rEquals2RbValue, DELTA);
  }

  void testCalculations_elliptical( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: ell=0.5 broken-exp with I_0 = 100, h1 = 50, h2 = 10, r_break = 50
    double  params[7] = {90.0, 0.5, 100.0, 50.0, 10.0, 50.0, 10.0};
    double  rEqualsOneValue, rEqualsRbValue, rEquals2RbValue;
    
    // test
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), 100.0, DELTA );
    // r = 1 value; account for ellipticity = 0.5 for y offsets
    rEqualsOneValue = 98.019867330675524;
    TS_ASSERT_DELTA( thisFunc->GetValue(101.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(99.0, 100.0), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 99.5), rEqualsOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.5), rEqualsOneValue, DELTA);
    // r = r_break value
    rEqualsRbValue = 36.584512991317212;
    TS_ASSERT_DELTA( thisFunc->GetValue(150.0, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(50.0, 100.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 75.0), rEqualsRbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 125.0), rEqualsRbValue, DELTA);
    // r = 2*r_break value
    rEquals2RbValue = 0.24787521766663584;
    TS_ASSERT_DELTA( thisFunc->GetValue(200.0, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 100.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 50.0), rEquals2RbValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 150.0), rEquals2RbValue, DELTA);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestBrokenExponential2D : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new BrokenExponential2D();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 7;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("h1");
    correctParamNames.push_back("h2");
    correctParamNames.push_back("r_break");
    correctParamNames.push_back("alpha");
    correctParamNames.push_back("h_z");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT_EQUALS(paramNames, correctParamNames );
    
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestDoubleBrokenExponential : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new DoubleBrokenExponential();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 10;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("h1");
    correctParamNames.push_back("h2");
    correctParamNames.push_back("h3");
    correctParamNames.push_back("r_break1");
    correctParamNames.push_back("r_break2");
    correctParamNames.push_back("alpha1");
    correctParamNames.push_back("alpha2");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT_EQUALS(paramNames, correctParamNames );
    
  }

// # Sample profile for tests
// # I0 = 100, h1 = 20, h2 = 10, h3 = 5, r_brk1 = 10, r_brk2 = 20
// # alpha = 0.1, 0.1
// # alpha = 1, 1
// # alpha = 100, 100
// 
// # p01 = [100.0, 20.0, 10.0, 5.0, 10.0, 20.0, 0.1, 0.1]
// # I_01 = [62.42459371, 62.00155144, 43.84819263, 29.89701942,19.75152889, 12.66637718, 0.62207092])
// # 
// # p1 = [100.0, 20.0, 10.0, 5.0, 10.0, 20.0, 1.0, 1.0]
// I1 = (100.0, 99.501224163558689, 77.85410746224612, 58.58686690528004, 36.750989624929915,
//  20.81878008682871, 0.055308562573433168)

// # p100 = [100.0, 20.0, 10.0, 5.0, 10.0, 20.0, 100.0, 100.0]
// # I_100 = (100.0, 99.501247919268238, 77.880078307140494, 60.632048862625965,
// # 17.377394345044511, 8.2028121352973606, 3.059023205018258e-05)
// r = np.array([0.0, 0.1, 5.0, 10.0, 15.0, 20.0, 50.0])

  void testCalculations_circular( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: circular doubly broken-exp with I_0 = 100, h1,h2,3 = 20,10,5, 
    //								r_brk1,r_brk2 = 10,20, alpha1,alpha2 = 1,1
    double  params[10] = {90.0, 0.0, 100.0, 20.0,10.0,5.0, 10.0,20.0, 1.0,1.0};
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), 100.0, DELTA );
    // r = 0.1 value
    double rEqualsZeroPointOneValue = 99.501224163558689;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.1, 100.0), rEqualsZeroPointOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(99.9, 100.0), rEqualsZeroPointOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 99.9), rEqualsZeroPointOneValue, DELTA);
    // r = 10 value
    double rEqualsTen = 58.58686690528004;
    TS_ASSERT_DELTA( thisFunc->GetValue(110.0, 100.0), rEqualsTen, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(90.0, 100.0), rEqualsTen, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 90.0), rEqualsTen, DELTA);
    // r = 50 value
    double rEquals50Value = 0.055308562573433168;
    TS_ASSERT_DELTA( thisFunc->GetValue(150.0, 100.0), rEquals50Value, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(50.0, 100.0), rEquals50Value, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 50.0), rEquals50Value, DELTA);
    // r = 500 value
    double rEquals500Value = 4.5319905994148505e-41;
    TS_ASSERT_EQUALS( thisFunc->GetValue(600.0, 100.0), rEquals500Value);
  }

  void testCalculations_elliptical( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: ell=0.5 broken-exp with I_0 = 100, h1,h2,3 = 20,10,5, 
    //								r_brk1,r_brk2 = 10,20, alpha1,alpha2 = 1,1
    double  params[10] = {90.0, 0.5, 100.0, 20.0,10.0,5.0, 10.0,20.0, 1.0,1.0};
    
    // test
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = 0 value
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), 100.0, DELTA );
    // r = 0.1 value; account for ellipticity = 0.5 for y offsets
    double rEqualsZeroPointOneValue = 99.501224163558689;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.1, 100.0), rEqualsZeroPointOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(99.9, 100.0), rEqualsZeroPointOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 99.95), rEqualsZeroPointOneValue, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.05), rEqualsZeroPointOneValue, DELTA);
    // r = 10 value; account for ellipticity = 0.5 for y offsets
    double rEqualsTen = 58.58686690528004;
    TS_ASSERT_DELTA( thisFunc->GetValue(110.0, 100.0), rEqualsTen, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(90.0, 100.0), rEqualsTen, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 95.0), rEqualsTen, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 105.0), rEqualsTen, DELTA);
    // r = 50 value; account for ellipticity = 0.5 for y offsets
    double rEquals50Value = 0.055308562573433168;
    TS_ASSERT_DELTA( thisFunc->GetValue(150.0, 100.0), rEquals50Value, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(50.0, 100.0), rEquals50Value, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 75.0), rEquals50Value, DELTA);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 125.0), rEquals50Value, DELTA);
    // r = 500 value; account for ellipticity = 0.5 for y offsets
    double rEquals500Value = 4.5319905994148505e-41;
    TS_ASSERT_EQUALS( thisFunc->GetValue(600.0, 100.0), rEquals500Value);
    TS_ASSERT_EQUALS( thisFunc->GetValue(100.0, 350.0), rEquals500Value);
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestEdgeOnDisk : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new EdgeOnDisk();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 5;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("L_0");
    correctParamNames.push_back("h");
    correctParamNames.push_back("n");
    correctParamNames.push_back("z_0");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }
  
  void testCalculations( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: EdgeOnDisk at PA = 90
    double  h = 10.0;
    double  L0 = 1.0;
    double  z0 = 1.0;
    double  n1 = 1.0;
    double  nexp = 1000.0;
    double  params[5] = {90.0, L0, h, n1, z0};
    
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // r = z = 0 value
    double  centralValue = 2.0 * h * L0;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), centralValue, DELTA );

    // vertical quasi-exponential limit: flux at r=0, z=10
    double  params2[5] = {90.0, L0, h, nexp, z0};
    double  offsetFlux = centralValue * pow(2.0, 2.0/nexp) * exp(-10.0/z0);
    thisFunc->Setup(params2, 0, x0, y0);
    double x = thisFunc->GetValue(100.0, 110.0);
    //printf("offsetFlux = %.18f, x = %.18f\n", offsetFlux, x);
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 110.0), offsetFlux, DELTA );

  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


class TestFerrersBar3D : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc, *thisFunc_subsampled;
  
public:
  void setUp()
  {
    // FUNCTION-SPECIFIC:
    bool  subsampleFlag = false;
    thisFunc = new FerrersBar3D();
    thisFunc->SetSubsampling(subsampleFlag);
  }
  
  void tearDown()
  {
    delete thisFunc;
  }


  // and now the actual tests
  void testBasic( void )
  {
    vector<string>  paramNames;
    vector<string>  correctParamNames;
    // FUNCTION-SPECIFIC:
    int  correctNParams = 8;
    correctParamNames.push_back("PA");
    correctParamNames.push_back("inc");
    correctParamNames.push_back("barPA");
    correctParamNames.push_back("J_0");
    correctParamNames.push_back("R_bar");
    correctParamNames.push_back("q");
    correctParamNames.push_back("q_z");
    correctParamNames.push_back("n");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames );
    
  }

  void testCalculations( void )
  {
    // centered at x0,y0 = 100,100
    double  x0 = 100.0;
    double  y0 = 100.0;
    // FUNCTION-SPECIFIC:
    // test setup: face-on FerrersBar3D
    // Face-on ellipsoid with a = b = c
    double  PA = 0.0;
    double  inc = 0.0;
    double  barPA = 0.0;
    double  J_0 = 1.0;
    double  R_bar = 100.0;
    double  q = 1.0;
    double  q_z = 1.0;
    double  n = 2.0;
    double  params[8] = {PA, inc, barPA, J_0, R_bar, q, q_z, n};
    
    thisFunc->Setup(params, 0, x0, y0);
    
    // FUNCTION-SPECIFIC:
    // central value (Python numerical integration gives 106.66666666666666)
    double  centralValue = 106.666666573;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), centralValue, DELTA );

	// value outside bar should = 0
	TS_ASSERT_DELTA( thisFunc->GetValue(0.0, 100.0), 0.0, DELTA );
	TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 0.0), 0.0, DELTA );
	TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 300.0), 0.0, DELTA );
	TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 201.0), 0.0, DELTA );
	TS_ASSERT_DELTA( thisFunc->GetValue(201.0, 100.0), 0.0, DELTA );
	
	// at r = R_bar/2 along major axis, integrated flux should be
	// (Python numerical integration gives 51.96152422430678)
	double  halfRadiusValue = 51.9615240239;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 + R_bar/2.0), halfRadiusValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 - R_bar/2.0, 100.0), halfRadiusValue, DELTA );
    
    // Now with b = a/2 (q = 0.5).
    double  params2[8] = {PA, inc, barPA, J_0, R_bar, 0.5, q_z, n};
    thisFunc->Setup(params2, 0, x0, y0);
    // Central flux unchanged
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), centralValue, DELTA );
    // Half-radius value from before is at R_bar/2 along major axis, R_bar/4 along minor axis
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0 - R_bar/2.0), halfRadiusValue, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 + R_bar/2.0, 100.0), 0.0, DELTA );
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0 + R_bar/4.0, 100.0), halfRadiusValue, DELTA );
  }

  void testCanCalculateTotalFlux( void )
  {
    bool result = thisFunc->CanCalculateTotalFlux();
    TS_ASSERT_EQUALS(result, false);
  }
};


// class TestSplineProfile : public CxxTest::TestSuite 
// {
// 
// public:
// 
//   // and now the actual tests
//   void testFuncCreation( void )
//   {
//     FunctionObject *theFunc;
//     string funcName, correctName;
//     
//     theFunc = new SplineProfile;
//     
//     TS_ASSERT_EQUALS( theFunc->HasExtraParams(), true );
//     TS_ASSERT_EQUALS( theFunc->CanCalculateTotalFlux(), false );
//     TS_ASSERT_EQUALS( theFunc->IsPointSource(), false );
// 
//     // MAX_SPLINE_POINTS + PA, ell, 
//     int  correctNParams = MAX_SPLINE_POINTS + 5;
//     TS_ASSERT_EQUALS( theFunc->GetNParams(), correctNParams );
//     
//     funcName = theFunc->GetShortName();
//     correctName = "SplineProfile";
//     TS_ASSERT_EQUALS( funcName, correctName );
// 
//     delete theFunc;
//   }
// 
//   void testExtraParams_goodInput( void )
//   {
//     FunctionObject *theFunc;
//     int  retVal;
//     map<string,string> extraParamsMap;
//     
//     extraParamsMap["N"] = "3";
//     
//     theFunc = new SplineProfile;
//     
//     TS_ASSERT_EQUALS( theFunc->HasExtraParams(), true );
//     TS_ASSERT_EQUALS( theFunc->ExtraParamsSet(), false );
//     int  correctNParams = MAX_SPLINE_POINTS + 3;
//     TS_ASSERT_EQUALS( theFunc->GetNParams(), correctNParams );
//     
//     retVal = theFunc->SetExtraParams(extraParamsMap);
//     TS_ASSERT_EQUALS( retVal, 1 );
//     TS_ASSERT_EQUALS( theFunc->ExtraParamsSet(), true );
//     correctNParams = 3 + 3;
//     TS_ASSERT_EQUALS( theFunc->GetNParams(), correctNParams );
//     
//     delete theFunc;
//   }
// 
//   void testExtraParams_badInput( void )
//   {
//     FunctionObject *theFunc;
//     int  retVal;
//     map<string,string> extraParamsMap;
//     
//     theFunc = new SplineProfile;
//     
//     TS_ASSERT_EQUALS( theFunc->HasExtraParams(), true );
//     TS_ASSERT_EQUALS( theFunc->ExtraParamsSet(), false );
//     int  correctNParams = 4 + 3;
//     TS_ASSERT_EQUALS( theFunc->GetNParams(), correctNParams );
//     
//     // empty map --> -1
//     retVal = theFunc->SetExtraParams(extraParamsMap);
//     TS_ASSERT_EQUALS( retVal, -1 );
//     
//     // map with invalid keyword --> 0
//     extraParamsMap["alpha"] = "3";
//     retVal = theFunc->SetExtraParams(extraParamsMap);
//     TS_ASSERT_EQUALS( retVal, 0 );
// 
//     // map with invalid value --> -3
//     extraParamsMap.erase("alpha");
//     extraParamsMap["N"] = "bob";
//     retVal = theFunc->SetExtraParams(extraParamsMap);
//     TS_ASSERT_EQUALS( retVal, -3 );
//     
//     delete theFunc;
//   }
//   
//   void testCalculations( void )
//   {
//     FunctionObject *theFunc;
// 
// //  * I_0 = params[0 + offsetIndex ]; -- position angle of major axis
// //  * I_0 = params[1 + offsetIndex ]; -- ellipticity of isophotes
// //  * I_0 = params[2 + offsetIndex ]; -- central surf. brightness
// //  * r_1 = params[3 + offsetIndex ];    -- radius of 2nd interpolation data point
// //  * I_1 = params[4 + offsetIndex ];    -- surf. brightness of 2nd interpolation data point
// //  * r_2 = params[5 + offsetIndex ];    -- radius of 3rd interpolation data point
// //  * I_2 = params[6 + offsetIndex ];    -- surf. brightness of 3rd interpolation data point
// //  * r_3 = params[7 + offsetIndex ];    -- radius of 4th interpolation data point
// //  * I_3 = params[8 + offsetIndex ];    -- surf. brightness of 4th interpolation data point
//     double  inputParams[MAX_SPLINE_POINTS + 3] = {100.0, 10.0, 50.0, 20.0, 0.0, 30.0, 0.0};
// 
//     theFunc = new SplineProfile;
//     theFunc->SetSubsampling(false);
//     theFunc->Setup(inputParams, 0, 0.0);
// 
//     // check that we recover interpolation reference values
//     TS_ASSERT_DELTA( theFunc->GetValue(0.0), 100.0, DELTA );
//     TS_ASSERT_DELTA( theFunc->GetValue(10.0), 50.0, DELTA );
//     TS_ASSERT_DELTA( theFunc->GetValue(20.0), 0.0, DELTA );
//     TS_ASSERT_DELTA( theFunc->GetValue(30.0), 0.0, DELTA );
// 
//     // check that we recover interpolated values
//     TS_ASSERT_DELTA( theFunc->GetValue(5.0), 76.25, DELTA );
//     TS_ASSERT_DELTA( theFunc->GetValue(15.0), 21.25, DELTA );
//     TS_ASSERT_DELTA( theFunc->GetValue(25.0), -5.0, DELTA );
//   }
// };
