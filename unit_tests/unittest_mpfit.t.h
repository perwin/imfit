// Unit/integration tests for mpfit.cpp

// See run_unittest_cmlineparser.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <stdlib.h>
#include <string>
using namespace std;
#include "mpfit.h"
#include "model_object.h"


// For use when calling mpfit
int myfunc_mpfit_good( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );

int myfunc_mpfit_good( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *theModel )
{

  theModel->ComputeDeviates(deviates, params);
  return 0;
}


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  double  *paramVector;
  double  *paramVector_null;
  mp_par  *parameterInfo;
  ModelObject  *theModel;
  mp_result  theResult;
  mp_config  *mpConfig;
  int  nData, nParams;

  void setUp()
  {
    int  i;
    nData = 5;
    nParams = 3;
    
    // set up parameter vector with all values = 1.0
    paramVector = (double *) calloc((size_t)nParams, sizeof(double));
    for (i = 0; i < nParams; i++)
      paramVector[i] = 1.0;
    // define paramVector_null to be null
    paramVector_null = NULL;
    
    parameterInfo = (mp_par *) calloc((size_t)nParams, sizeof(mp_par));
    // set lower & upper limits = 0,100 for all parameters
    for (i = 0; i < nParams; i++) {
      parameterInfo[i].limited[0] = 1;
      parameterInfo[i].limits[0] = 0.0;
      parameterInfo[i].limited[1] = 1;
      parameterInfo[i].limits[1] = 100.0;
    }

    mpConfig = (mp_config *) calloc((size_t)nParams, sizeof(mp_config));
  }
  
  void tearDown()
  {
     free(paramVector);
     free(parameterInfo);
     free(mpConfig);
  }
  
  


//Tests for mpfit() -- main we just test to see if various input errors are caught;
// formal testing of the actual minimization process is too involved for unit tests.

  void testMPFit_noUserFunc( void )
  {
    int  status;
    int  expectedStatus = MP_ERR_FUNC;
    mp_func  badFunc = 0;
    
    status = mpfit(badFunc, nData, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
  }

  void testMPFit_noDataPoints( void )
  {
    int  status;
    int  expectedStatus = MP_ERR_NPOINTS;
    
    // case of nData <= 0
    status = mpfit(myfunc_mpfit_good, 0, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
    status = mpfit(myfunc_mpfit_good, -1000, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);

    // case of null/unallocated data vector
    status = mpfit(myfunc_mpfit_good, nData, nParams, paramVector_null, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
  }

  void testMPFit_noParams( void )
  {
    int  status;
    int  expectedStatus = MP_ERR_NFREE;
    
    // case of nParams <= 0
    status = mpfit(myfunc_mpfit_good, nData, 0, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
    status = mpfit(myfunc_mpfit_good, nData, -1000, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
  }

  void testMPFit_allParamsFixed( void )
  {
    int  status;
    int  expectedStatus = MP_ERR_NFREE;
    
    // we have parameters, but all are fixed
    for (int i = 0; i < nParams; i++)
      parameterInfo[i].fixed = 1;
    
    status = mpfit(myfunc_mpfit_good, nData, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
  }

  void testMPFit_inconsistentBounds( void )
  {
    int  status;
    mp_par  *badParameterInfo;
    int  expectedStatus = MP_ERR_BOUNDS;
    
    badParameterInfo = (mp_par *) calloc((size_t)nParams, sizeof(mp_par));
    // start off with OK bounds (lower,upper = 0,100)
    for (int i = 0; i < nParams; i++) {
      badParameterInfo[i].limited[0] = 1;
      badParameterInfo[i].limits[0] = 0.0;
      badParameterInfo[i].limited[1] = 1;
      badParameterInfo[i].limits[1] = 100.0;
    }
    // inconsistent: lower > upper
    badParameterInfo[0].limits[0] = 100.0;
    badParameterInfo[0].limits[1] = 0.0;
    status = mpfit(myfunc_mpfit_good, nData, nParams, paramVector, badParameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
    
    free(badParameterInfo);
  }

  void testMPFit_paramsOutsideBounds( void )
  {
    int  status;
    int  expectedStatus = MP_ERR_INITBOUNDS;
    
    // input parameter below limit
    paramVector[0] = -1.0;
    status = mpfit(myfunc_mpfit_good, nData, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
    // input parameter above limit
    paramVector[0] = 101.0;
    status = mpfit(myfunc_mpfit_good, nData, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
    TS_ASSERT_EQUALS(status, expectedStatus);
    
    // reset parameter values for other tests
    paramVector[0] = 0.0;
  }

  void testMPFit_enoughDataPoints( void )
  {
    int  status, nDataCurrent;
    int  expectedStatus = MP_ERR_DOF;

    // case 1: too few data points
    nDataCurrent = nParams - 1;
    status = mpfit(myfunc_mpfit_good, nDataCurrent, nParams, paramVector, parameterInfo,
    				mpConfig, theModel, &theResult);
  }

// Things we can't test:
//   Return value of MP_ERR_PARAM, since mpfit.cpp is set up to prevent this from
// ever happening (or else catches the error previously)



  // Tests for InterpretMpfitResult()
  void testInterpretMpfitResult_failures( void )
  {
    string  outputString, expectedString;
    int  mpfitResult;
    
    expectedString = "ERROR: General input parameter error!";
    InterpretMpfitResult(MP_ERR_INPUT, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: User function produced non-finite values!";
    InterpretMpfitResult(MP_ERR_NAN, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: No user function was supplied!";
    InterpretMpfitResult(MP_ERR_FUNC, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: No user data points were supplied!";
    InterpretMpfitResult(MP_ERR_NPOINTS, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: No free parameters!";
    InterpretMpfitResult(MP_ERR_NFREE, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: Memory allocation error!";
    InterpretMpfitResult(MP_ERR_MEMORY, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: Initial values inconsistent with constraints!";
    InterpretMpfitResult(MP_ERR_INITBOUNDS, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: Initial constraints inconsistent!";
    InterpretMpfitResult(MP_ERR_BOUNDS, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: General input parameter error!";
    InterpretMpfitResult(MP_ERR_PARAM, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "ERROR: Not enough degrees of freedom!";
    InterpretMpfitResult(MP_ERR_DOF, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);
  }
  
  void testInterpretMpfitResult_successes( void )
  {
    string  outputString, expectedString;
    int  mpfitResult;
    
    expectedString = "SUCCESS: Convergence in fit-statistic value.";
    InterpretMpfitResult(MP_OK_CHI, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "SUCCESS: Convergence in parameter value.";
    InterpretMpfitResult(MP_OK_PAR, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "SUCCESS: Convergence in fit-statistic and parameter value.";
    InterpretMpfitResult(MP_OK_BOTH, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "SUCCESS: Convergence in orthogonality.";
    InterpretMpfitResult(MP_OK_DIR, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);
  }

  void testInterpretMpfitResult_terminations( void )
  {
    string  outputString, expectedString;
    int  mpfitResult;
    
    expectedString = "Terminated: Maximum number of iterations reached";
    InterpretMpfitResult(MP_MAXITER, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "Terminated: ftol too small; no further improvement";
    InterpretMpfitResult(MP_FTOL, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "Terminated: xtol too small; no further improvement";
    InterpretMpfitResult(MP_XTOL, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);

    expectedString = "Terminated: gtol too small; no further improvement";
    InterpretMpfitResult(MP_GTOL, outputString);
    TS_ASSERT_EQUALS(outputString, expectedString);
  }
};

