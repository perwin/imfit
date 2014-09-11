// Unit tests for code in model_object.cpp

// See run_unittest_model_object.sh for how to compile & run this


// Things to test:


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
using namespace std;
#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object.h"
#include "add_functions.h"
#include "config_file_parser.h"
#include "param_struct.h"


#define CONFIG_FILE "tests/config_imfit_poisson.dat"


// Reference things
const string  headerLine_correct = "# X0_1		Y0_1		PA_1	ell_1	I_0_1	h_1	I_sky_2	";


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  ModelObject *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  FunctionBlockIndices;
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
    
    // Set up internal FunctionObjects vector (Exp + FlatSky) inside theModel
    status = ReadConfigFile(filename, true, functionList, parameterList, 
  								paramLimits, FunctionBlockIndices, paramLimitsExist, userConfigOptions);
    status = AddFunctions(theModel, functionList, FunctionBlockIndices, true);

  }

  
  void testGetParamHeader( void )
  {
    std::string  headerLine;
    headerLine = theModel->GetParamHeader();
    TS_ASSERT_EQUALS(headerLine, headerLine_correct);

  }
  
  void testSetAndGetStatisticType( void )
  {
    int  whichStat;
    bool cashStatUsed = false;
    
    // default should be chi^2
    whichStat = theModel->WhichFitStatistic();
    TS_ASSERT_EQUALS(whichStat, FITSTAT_CHISQUARE);
    
    theModel->UseCashStatistic();
    whichStat = theModel->WhichFitStatistic();
    TS_ASSERT_EQUALS(whichStat, FITSTAT_CASH);
    cashStatUsed = theModel->UsingCashStatistic();
    TS_ASSERT_EQUALS(cashStatUsed, true);

    // the following will generate a couple of warnings from within UseCashStatistic,
    // which is OK
    theModel->UseModifiedCashStatistic();
    whichStat = theModel->WhichFitStatistic();
    TS_ASSERT_EQUALS(whichStat, FITSTAT_MODCASH);
    cashStatUsed = theModel->UsingCashStatistic();
    TS_ASSERT_EQUALS(cashStatUsed, true);
  }
  

};
