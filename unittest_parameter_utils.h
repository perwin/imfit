// Unit tests for condensation/expansion of parameter sets
//
// $ cxxtestgen --error-printer -o test_runner.cpp unittest_parameter_utils.h
// $ g++ -o test_runner test_runner.cpp parameter_utils.cpp -I/usr/local/include  -I$CXXTEST
// ./test_runner


#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>
using namespace std;

#include "parameter_utils.h"
#include "param_struct.h"   // for mp_par structure

#define DELTA  1.0e-10


// Testing of CountFixedParams

class TestCountFixedParams : public CxxTest::TestSuite 
{
  int  nParamsTot;
  mp_par *parameterInfo;
  
public:
  void setUp()
  {
    nParamsTot = 4;
    parameterInfo = (mp_par *) calloc((size_t)nParamsTot, sizeof(mp_par));
    for (int i = 0; i < nParamsTot; i++) {
      parameterInfo[i].fixed = 0;
    }
  }
  
  void tearDown()
  {
    free(parameterInfo);
  }


  // and now the actual tests
  // Test handling of 1 fixed parameter
  void testCondenseParams1( void )
  {
    int  nFixed;
    
    parameterInfo[0].fixed = 1;
    nFixed = CountFixedParams(parameterInfo, nParamsTot);
    TS_ASSERT_EQUALS(nFixed, 1);

    parameterInfo[0].fixed = 0;
    parameterInfo[3].fixed = 1;
    nFixed = CountFixedParams(parameterInfo, nParamsTot);
    TS_ASSERT_EQUALS(nFixed, 1);
  }

  // Test handling of 2 fixed parameters
  void testCondenseParams2( void )
  {
    int  nFixed;
    
    parameterInfo[0].fixed = 1;
    parameterInfo[1].fixed = 1;
    parameterInfo[2].fixed = 0;
    parameterInfo[3].fixed = 0;
    nFixed = CountFixedParams(parameterInfo, nParamsTot);
    TS_ASSERT_EQUALS(nFixed, 2);

    parameterInfo[0].fixed = 0;
    parameterInfo[1].fixed = 0;
    parameterInfo[2].fixed = 1;
    parameterInfo[3].fixed = 1;
    nFixed = CountFixedParams(parameterInfo, nParamsTot);
    TS_ASSERT_EQUALS(nFixed, 2);

    parameterInfo[0].fixed = 1;
    parameterInfo[1].fixed = 0;
    parameterInfo[2].fixed = 0;
    parameterInfo[3].fixed = 1;
    nFixed = CountFixedParams(parameterInfo, nParamsTot);
    TS_ASSERT_EQUALS(nFixed, 2);
  }

};




// Testing of CondenseParamVector

class TestCondenseParamVector : public CxxTest::TestSuite 
{
  int  nParamsTot;
  double  *origParams;
  bool  *fixedParamFlags;
  
public:
  void setUp()
  {
    nParamsTot = 4;
    
    origParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
    origParams[0] = 0.0;
    origParams[1] = 4.0;
    origParams[2] = 20.0;
    origParams[3] = 25.0;
    
    fixedParamFlags = (bool *)malloc((size_t)(nParamsTot*sizeof(bool)));

  }
  
  void tearDown()
  {
    free(origParams);
    free(fixedParamFlags);
  }


  // and now the actual tests
  
  // Test handling of 0 fixed parameters
  void testCondenseParams0( void )
  {
    double *freeParams;
    int  nFree;

	nFree = 4;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );

    // first test: initial parameter is fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[0]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[1]);
    TS_ASSERT_EQUALS(freeParams[2], origParams[2]);
    TS_ASSERT_EQUALS(freeParams[3], origParams[3]);
  }
  
  // Test handling of 1 fixed parameter
  void testCondenseParams1( void )
  {
    double *freeParams;
    int  nFree;

	nFree = 3;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );

    // first test: initial parameter is fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[1]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[2]);
    TS_ASSERT_EQUALS(freeParams[2], origParams[3]);

    // second test: last parameter is fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = true;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[0]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[1]);
    TS_ASSERT_EQUALS(freeParams[2], origParams[2]);

    // third test: second parameter is fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = true;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[0]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[2]);
    TS_ASSERT_EQUALS(freeParams[2], origParams[3]);

    
    free(freeParams);
  }

  // test handling of 2 fixed parameters
  void testCondenseParams2( void )
  {
    double *freeParams;
    int  nFree;

	nFree = 2;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );

    // first test: first 2 parameters are fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = true;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[2]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[3]);

    // second test: last 2 parameters are fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = true;
    fixedParamFlags[3] = true;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[0]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[1]);

    // third test: first and last parameters are fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = true;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[1]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[2]);

    // fourth test: first and third parameters are fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = true;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    TS_ASSERT_EQUALS(freeParams[0], origParams[1]);
    TS_ASSERT_EQUALS(freeParams[1], origParams[3]);

    
    free(freeParams);
  }
  
};




// Testing of ExpandParamVector
// void ExpandParamVector( double originalParams[], double condensedInputParams[],
// 							double outputParams[], int nParamsTot, int nParamsFree,
// 							bool fixedPars[] );

class TestExpandParamVector : public CxxTest::TestSuite 
{
  int  nParamsTot, nParamsFree;
  double  *origParams;
  bool  *fixedParamFlags;
  
public:
  void setUp()
  {
    nParamsTot = 4;
    
    origParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
    origParams[0] = 1.0;
    origParams[1] = 2.0;
    origParams[2] = 20.0;
    origParams[3] = 30.0;
    
    fixedParamFlags = (bool *)malloc((size_t)(nParamsTot*sizeof(bool)));

  }
  
  void tearDown()
  {
    free(origParams);
    free(fixedParamFlags);
  }


  // and now the actual tests
  
  // Test handling of 0 fixed parameters
  void testExpandParams0( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 4;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // first test: condense & restore original parameter set; *no* fixed params
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    

    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 1 fixed parameter, v1
  void testExpandParams1a( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 3;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // first test: condense & restore original parameter set; first param is fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 1 fixed parameter, v2
  void testExpandParams1b( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 3;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // second test: condense & restore original parameter set; last param is fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = true;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 1 fixed parameter, v3
  void testExpandParams1c( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 3;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // third test: condense & restore original parameter set; last param is fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = true;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 2 fixed parameters, v1
  void testExpandParams2a( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 2;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // first test: condense & restore original parameter set; first param is fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = true;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 2 fixed parameters, v2
  void testExpandParams2b( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 2;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // first test: condense & restore original parameter set; first param is fixed
    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = true;
    fixedParamFlags[3] = true;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 2 fixed parameters, v3
  void testExpandParams2c( void )
  {
    double *freeParams;
    double *outputFullParams;
    int  nFree;

	nFree = 2;
    freeParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    // first test: condense & restore original parameter set; first param is fixed
    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = true;

    CondenseParamVector(origParams, freeParams, nParamsTot, nFree, fixedParamFlags);
    
    ExpandParamVector(origParams, freeParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    free(freeParams);
    free(outputFullParams);
  }

  // Test handling of 1 fixed parameter and new input params
  void testExpandParams_newInput1a( void )
  {
    double *inputParams;
    double *outputFullParams;
    int  nFree;
    double  x1 = 100.0;
    double  x2 = 200.0;
    double  x3 = 300.0;

	nFree = 3;
    inputParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;
    
    inputParams[0] = x1;
    inputParams[1] = x2;
    inputParams[2] = x3;
    
    ExpandParamVector(origParams, inputParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], x1);
    TS_ASSERT_EQUALS(outputFullParams[2], x2);
    TS_ASSERT_EQUALS(outputFullParams[3], x3);
    
    
    free(inputParams);
    free(outputFullParams);
  }

  // Test handling of 1 fixed parameter and new input params
  void testExpandParams_newInput1b( void )
  {
    double *inputParams;
    double *outputFullParams;
    int  nFree;
    double  x1 = 100.0;
    double  x2 = 200.0;
    double  x3 = 300.0;

	nFree = 3;
    inputParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = true;
    
    inputParams[0] = x1;
    inputParams[1] = x2;
    inputParams[2] = x3;
    
    ExpandParamVector(origParams, inputParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], x1);
    TS_ASSERT_EQUALS(outputFullParams[1], x2);
    TS_ASSERT_EQUALS(outputFullParams[2], x3);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    
    free(inputParams);
    free(outputFullParams);
  }

  // Test handling of 1 fixed parameter and new input params, v3
  void testExpandParams_newInput1c( void )
  {
    double *inputParams;
    double *outputFullParams;
    int  nFree;
    double  x1 = 100.0;
    double  x2 = 200.0;
    double  x3 = 300.0;

	nFree = 3;
    inputParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = true;
    fixedParamFlags[3] = false;
    
    inputParams[0] = x1;
    inputParams[1] = x2;
    inputParams[2] = x3;
    
    ExpandParamVector(origParams, inputParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], x1);
    TS_ASSERT_EQUALS(outputFullParams[1], x2);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], x3);
    
    
    free(inputParams);
    free(outputFullParams);
  }

  // Test handling of 2 fixed parameters and new input params
  void testExpandParams_newInput2a( void )
  {
    double *inputParams;
    double *outputFullParams;
    int  nFree;
    double  x1 = 100.0;
    double  x2 = 200.0;

	nFree = 2;
    inputParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    fixedParamFlags[0] = true;
    fixedParamFlags[1] = true;
    fixedParamFlags[2] = false;
    fixedParamFlags[3] = false;
    
    inputParams[0] = x1;
    inputParams[1] = x2;
    
    ExpandParamVector(origParams, inputParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], origParams[1]);
    TS_ASSERT_EQUALS(outputFullParams[2], x1);
    TS_ASSERT_EQUALS(outputFullParams[3], x2);
    
    
    free(inputParams);
    free(outputFullParams);
  }

  // Test handling of 2 fixed parameters and new input params, v2
  void testExpandParams_newInput2b( void )
  {
    double *inputParams;
    double *outputFullParams;
    int  nFree;
    double  x1 = 100.0;
    double  x2 = 200.0;

	nFree = 2;
    inputParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    fixedParamFlags[0] = false;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = true;
    fixedParamFlags[3] = true;
    
    inputParams[0] = x1;
    inputParams[1] = x2;
    
    ExpandParamVector(origParams, inputParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], x1);
    TS_ASSERT_EQUALS(outputFullParams[1], x2);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], origParams[3]);
    
    
    free(inputParams);
    free(outputFullParams);
  }

  // Test handling of 2 fixed parameters and new input params, v3
  void testExpandParams_newInput2c( void )
  {
    double *inputParams;
    double *outputFullParams;
    int  nFree;
    double  x1 = 100.0;
    double  x2 = 200.0;

	nFree = 2;
    inputParams = (double *)calloc( (size_t)nFree, sizeof(double) );
    outputFullParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

    fixedParamFlags[0] = true;
    fixedParamFlags[1] = false;
    fixedParamFlags[2] = true;
    fixedParamFlags[3] = false;
    
    inputParams[0] = x1;
    inputParams[1] = x2;
    
    ExpandParamVector(origParams, inputParams, outputFullParams, nParamsTot, nFree, fixedParamFlags);
    TS_ASSERT_EQUALS(outputFullParams[0], origParams[0]);
    TS_ASSERT_EQUALS(outputFullParams[1], x1);
    TS_ASSERT_EQUALS(outputFullParams[2], origParams[2]);
    TS_ASSERT_EQUALS(outputFullParams[3], x2);
    
    
    free(inputParams);
    free(outputFullParams);
  }



};
