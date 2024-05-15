// Unit/integration tests for add_functions.cpp

// See run_unittest_add_functions.sh for how to compile & run the tests

// Note that we don't try to test PrintAvailableFunctions() or ListFunctionParameters(),
// because these print to standard output (and can in principle be tested in the
// regression tests).

#include <cxxtest/TestSuite.h>

#include <string>
using namespace std;
#include "add_functions.h"
#include "model_object.h"
#include "config_file_parser.h"

// TODO: fix testAddFunctionsToModel_optionalParams so that it uses an
// actual optional-parameter function (not just pretending that FlatSky
// takes optional parameters). PointSource is a possibility (though we
// don't actually make any use of its optional parameter at the moment...)

#define SIMPLE_CONFIG_FILE "tests/imfit_reference/config_imfit_flatsky.dat"
#define SIMPLE_CONFIG_FILE2 "tests/makeimage_reference/config_imfit_pointsource-extra.dat"


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
    referenceFunctionNameList.push_back("DattathriPeanut3D");
    referenceFunctionNameList.push_back("EdgeOnDisk");
    referenceFunctionNameList.push_back("EdgeOnRing");
    referenceFunctionNameList.push_back("EdgeOnRing2Side");
    referenceFunctionNameList.push_back("Exponential");
    referenceFunctionNameList.push_back("ExponentialDisk3D");
    referenceFunctionNameList.push_back("Exponential_GenEllipse");
    referenceFunctionNameList.push_back("FerrersBar2D");
    referenceFunctionNameList.push_back("FerrersBar3D");
    referenceFunctionNameList.push_back("FlatBar");
    referenceFunctionNameList.push_back("FlatSky");
    referenceFunctionNameList.push_back("Gaussian");
    referenceFunctionNameList.push_back("GaussianRing");
    referenceFunctionNameList.push_back("GaussianRing2Side");
    referenceFunctionNameList.push_back("GaussianRing3D");
    referenceFunctionNameList.push_back("GaussianRingAz");
    referenceFunctionNameList.push_back("ModifiedKing");
    referenceFunctionNameList.push_back("ModifiedKing2");
    referenceFunctionNameList.push_back("Moffat");
    referenceFunctionNameList.push_back("PointSource");
    referenceFunctionNameList.push_back("PointSourceRot");
    referenceFunctionNameList.push_back("Sersic");
    referenceFunctionNameList.push_back("Sersic_GenEllipse");
    referenceFunctionNameList.push_back("TiltedSkyPlane");
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



  void testAddFunctionsToModel( void )
  {
    ModelObject *modelObj;
    vector<string>  fnameList;
    vector<string>  flabelList;
    vector<int> funcSetIndices;
    vector<double> parameterList;
    vector<mp_par>  paramLimits;
    bool  paramLimitsExist;
    configOptions  userConfigOptions;
    string filename = SIMPLE_CONFIG_FILE;
    int  status, nInputFuncs, nOutputFuncs;

    status = ReadConfigFile(filename, true, fnameList, flabelList, parameterList, paramLimits,
  							funcSetIndices, paramLimitsExist, userConfigOptions);

    modelObj = new ModelObject();
    status = AddFunctions(modelObj, fnameList, flabelList, funcSetIndices, false, -1);
    TS_ASSERT_EQUALS(status, 0);

    vector<string> outputFuncNames;
    modelObj->GetFunctionNames(outputFuncNames);
    nInputFuncs = fnameList.size();
    nOutputFuncs = outputFuncNames.size();
    TS_ASSERT_EQUALS(nOutputFuncs, nInputFuncs);
    if ((nInputFuncs > 0) && (nInputFuncs == nOutputFuncs)) {
      for (int i = 0; i < nInputFuncs; i++)
        TS_ASSERT_EQUALS(fnameList[i], outputFuncNames[i]);
    }
  }
  
  void testAddFunctionsToModel_with_labels( void )
  {
    ModelObject *modelObj;
    vector<string>  fnameList;
    vector<string>  flabelList;
    vector<int> funcSetIndices;
    vector<double> parameterList;
    vector<mp_par>  paramLimits;
    bool  paramLimitsExist;
    configOptions  userConfigOptions;
    string filename = SIMPLE_CONFIG_FILE;
    int  status, nInputFuncs, nOutputFuncs;

    status = ReadConfigFile(filename, true, fnameList, flabelList, parameterList, paramLimits,
  							funcSetIndices, paramLimitsExist, userConfigOptions);

    modelObj = new ModelObject();
    status = AddFunctions(modelObj, fnameList, flabelList, funcSetIndices, false, -1);
    TS_ASSERT_EQUALS(status, 0);

    vector<string> outputFuncNames;
    modelObj->GetFunctionNames(outputFuncNames);
    nInputFuncs = fnameList.size();
    nOutputFuncs = outputFuncNames.size();
    TS_ASSERT_EQUALS(nOutputFuncs, nInputFuncs);
    if ((nInputFuncs > 0) && (nInputFuncs == nOutputFuncs)) {
      for (int i = 0; i < nInputFuncs; i++)
        TS_ASSERT_EQUALS(fnameList[i], outputFuncNames[i]);
    }

    vector<string> outputFuncLabels;
    modelObj->GetFunctionLabels(outputFuncLabels);
    nInputFuncs = fnameList.size();
    nOutputFuncs = outputFuncNames.size();
    if ((nInputFuncs > 0) && (nInputFuncs == nOutputFuncs)) {
      for (int i = 0; i < nInputFuncs; i++)
        TS_ASSERT_EQUALS(flabelList[i], outputFuncLabels[i]);
    }
  }
  
    void testAddFunctionsToModel_optionalParams( void )
  {
    ModelObject *modelObj;
    vector<string>  fnameList;
    vector<string>  flabelList;
    vector<int> funcSetIndices;
    vector<double> parameterList;
    vector<mp_par>  paramLimits;
    bool  paramLimitsExist;
    configOptions  userConfigOptions;
    string filename = SIMPLE_CONFIG_FILE2;
    int  status, nInputFuncs, nOutputFuncs;

    vector< map<string, string> > optionalParamsVect;
    map<string, string> optionalParamsMap;
    string  keyword = "method";
    string  value = "bicubic";
    optionalParamsMap[keyword] = value;
    optionalParamsVect.push_back(optionalParamsMap);

    status = ReadConfigFile(filename, true, fnameList, flabelList, parameterList, paramLimits,
  							funcSetIndices, paramLimitsExist, userConfigOptions);

    modelObj = new ModelObject();
    // need to add a PSF image, otherwise ModelObject will object when we try to add
    // a PointSource function
    int  nColumns_psf = 4;
    int  nRows_psf = 4;
    int  nPixels_psf = nColumns_psf*nRows_psf;
    double psfPixels[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    status = modelObj->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels, false);
    TS_ASSERT_EQUALS(status, 0);
    status = AddFunctions(modelObj, fnameList, flabelList, funcSetIndices, false, -1, 
    						optionalParamsVect);
    TS_ASSERT_EQUALS(status, 0);

    vector<string> outputFuncNames;
    modelObj->GetFunctionNames(outputFuncNames);
    nInputFuncs = fnameList.size();
    nOutputFuncs = outputFuncNames.size();
    TS_ASSERT_EQUALS(nOutputFuncs, nInputFuncs);
    if ((nInputFuncs > 0) && (nInputFuncs == nOutputFuncs)) {
      for (int i = 0; i < nInputFuncs; i++)
        TS_ASSERT_EQUALS(fnameList[i], outputFuncNames[i]);
    }
  }

};
