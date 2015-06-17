// Unit tests for 2D function classes
//
// "FUNCTION-SPECIFIC" denotes code which needs to be changed when testing
// a different function object (e.g., Sersic vs Exponential, etc.).


#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>
#include <vector>
using namespace std;

#include "function_objects/function_object.h"
// FUNCTION-SPECIFIC:
#include "function_objects/func_exp.h"
#include "function_objects/func_gaussian.h"
#include "function_objects/func_sersic.h"
#include "function_objects/func_edge-on-disk.h"

#define DELTA  1.0e-10


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
    double  params[5] = {90.0, 0.0, 1.0, 10.0};
    
    
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
};
