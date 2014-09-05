// Unit tests for code in model_object.cpp

// See run_unittest_model_object.sh for how to compile & run this


// Things to test:


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
using namespace std;
#include "function_objects/function_object.h"
#include "model_object.h"
#include "add_functions.h"
#include "config_file_parser.h"
#include "param_struct.h"


#define CONFIG_FILE "tests/config_imfit_poisson.dat"


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  ModelObject *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  functionSetIndices;
  bool  paramLimitsExist;
  mp_par  *parameterInfo;
  int  status;
  configOptions  userConfigOptions;

  void setUp()
  {
    int  status;
    string  filename = CONFIG_FILE;
    
    theModel = new ModelObject();
    paramLimitsExist = false;
    
    status = ReadConfigFile(filename, true, functionList, parameterList, 
  								paramLimits, functionSetIndices, paramLimitsExist, userConfigOptions);
    status = AddFunctions(theModel, functionList, functionSetIndices, true);

  }

  // Tests for StripLeadingDashes function
  
  void testGetParamHeader( void )
  {
    std::string  headerLine;
    headerLine = theModel->GetParamHeader();
    printf("%s\n", headerLine.c_str());

  }
  
  

};
