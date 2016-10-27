// Unit tests for 2D function classes
//
// "FUNCTION-SPECIFIC" denotes code which needs to be changed when testing
// a different function object (e.g., Sersic vs Exponential, etc.).

// See run_unittest_funcs.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>
#include <vector>
using namespace std;

#include "function_objects/function_object.h"
// FUNCTION-SPECIFIC:
#include "function_objects_1d/func1d_exp_test.h"
#include "function_objects/func_flatsky.h"
#include "function_objects/func_exp.h"
#include "function_objects/func_gaussian.h"
#include "function_objects/func_moffat.h"
#include "function_objects/func_sersic.h"
#include "function_objects/func_king.h"
#include "function_objects/func_king2.h"
#include "function_objects/func_broken-exp.h"
#include "function_objects/func_broken-exp2d.h"
#include "function_objects/func_edge-on-disk.h"

const double  DELTA = 1.0e-9;
const double  DELTA_e9 = 1.0e-9;
const double  DELTA_e7 = 1.0e-7;

const double PI = 3.14159265358979;


class TestExp1DTest : public CxxTest::TestSuite 
{
  FunctionObject  *thisFunc;

public:
  void setUp()
  {
    bool  subsampleFlag = false;
    thisFunc = new Exponential1D_test();
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
    int  correctNParams = 2;
    correctParamNames.push_back("mu_0");
    correctParamNames.push_back("h");

    // check that we get right number of parameters
    TS_ASSERT_EQUALS( thisFunc->GetNParams(), correctNParams );

    // check that we get correct set of parameter names
    thisFunc->GetParameterNames(paramNames);
    TS_ASSERT( paramNames == correctParamNames ); 
  }
  
  void testSetExtraParams( void )
  {
    map<string, string> theMap;
    string  keyword = "interp";
    string  value = "cspline";
    theMap[keyword] = value;
    
    int  returnVal = thisFunc->SetExtraParams(theMap);
    TS_ASSERT_EQUALS( returnVal, 0 );
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
//    subsampleFlag = true;
//    thisFunc_subsampled = new Exponential(subsampleFlag);
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
//    subsampleFlag = true;
//    thisFunc_subsampled = new Sersic(subsampleFlag);
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
    TS_ASSERT( paramNames == correctParamNames );
    
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
    double  centralFlux = 2.0 * h * L0;
    TS_ASSERT_DELTA( thisFunc->GetValue(100.0, 100.0), centralFlux, DELTA );

    // vertical quasi-exponential limit: flux at r=0, z=10
    double  params2[5] = {90.0, L0, h, nexp, z0};
    double  offsetFlux = centralFlux * pow(2.0, 2.0/nexp) * exp(-10.0/z0);
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
