// Unit/integration tests for add_functions.cpp

// See run_unittest_add_functions.sh for how to compile & run the tests (assuming
// all of the necessary func_*.cpp files, along with model_object.cpp, have already 
// been compiled).

// Note that we don't try to test PrintAvailableFunctions() or ListFunctionParameters(),
// because these print to standard output (and can in principle be tested in the
// regression tests).

#include <cxxtest/TestSuite.h>

#include <string>
using namespace std;
#include "add_functions.h"


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  vector<string>  referenceFunctionNameList;
  int  nFunctionNames_ref;

  void setUp()
  {
    referenceFunctionNameList.clear();
    referenceFunctionNameList.push_back("BrokenExponential");
    referenceFunctionNameList.push_back("BrokenExponential2D");
    referenceFunctionNameList.push_back("BrokenExponentialDisk3D");
    referenceFunctionNameList.push_back("Core-Sersic");
    referenceFunctionNameList.push_back("EdgeOnDisk");
    referenceFunctionNameList.push_back("EdgeOnRing");
    referenceFunctionNameList.push_back("EdgeOnRing2Side");
    referenceFunctionNameList.push_back("Exponential");
    referenceFunctionNameList.push_back("ExponentialDisk3D");
    referenceFunctionNameList.push_back("Exponential_GenEllipse");
    referenceFunctionNameList.push_back("FlatSky");
    referenceFunctionNameList.push_back("Gaussian");
    referenceFunctionNameList.push_back("GaussianRing");
    referenceFunctionNameList.push_back("GaussianRing2Side");
    referenceFunctionNameList.push_back("GaussianRing3D");
    referenceFunctionNameList.push_back("ModifiedKing");
    referenceFunctionNameList.push_back("Moffat");
    referenceFunctionNameList.push_back("Sersic");
    referenceFunctionNameList.push_back("Sersic_GenEllipse");
    nFunctionNames_ref = referenceFunctionNameList.size();
  }

  // Tests for GetFunctionParameterNames()
  void testBadFunctionName( void )
  {
    string  badFname = "somefunction_that_doesntexist";
    vector<string>  paramNameList;
    int  result;
    
    result = GetFunctionParameterNames(badFname, paramNameList);
    TS_ASSERT_EQUALS(result, -1);
  }

  void testGoodFunctionName( void )
  {
    string  fname = "Exponential";
    vector<string>  paramNameList;
    vector<string>  correctParamNames;
    int  result;
    
    correctParamNames.push_back("PA");
    correctParamNames.push_back("ell");
    correctParamNames.push_back("I_0");
    correctParamNames.push_back("h");
    
    result = GetFunctionParameterNames(fname, paramNameList);
    TS_ASSERT_EQUALS(result, 0);
    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(paramNameList[i], correctParamNames[i]);
  }


  void testGetFunctionNameList( void )
  {
    vector<string>  fnameList;
    int  nFuncsReturned;
    
    GetFunctionNames(fnameList);
    
    nFuncsReturned = fnameList.size();
    TS_ASSERT_EQUALS(nFuncsReturned, nFunctionNames_ref);
    for (int i = 0; i < nFunctionNames_ref; i++)
      TS_ASSERT_EQUALS(fnameList[i], referenceFunctionNameList[i]);
  }
};

