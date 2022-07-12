// Unit tests for code in config_file_parser.cpp

// See run_unittest_configfileparser.sh for how to compile and run these tests.

// older compilation notes:
// $ cxxtestgen.py --error-printer -o test_runner.cpp unittest_config_parser.h
// $ g++ -o test_runner test_runner.cpp utilities.cpp config_file_parser.cpp -I/usr/local/include

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

#include "utilities_pub.h"
#include "config_file_parser.h"
#include "param_struct.h"


// Contents of good parameter file
// # Config file for makeimage, meant to be used to generate a 512x512 image.
// 
// NCOLS   512
// NROWS   512
// 
// # Generalized-ellipse Sersic object
// X0    256.0
// Y0    256.0
// FUNCTION   Sersic_GenEllipse
// PA    30.0
// ell    0.5
// c0    -1.0
// n      2.5
// I_e  100.0
// r_e   50.0

const string  TEST_CONFIGFILE_GOOD("./tests/config_makeimage_sersictest512.dat");
const string  TEST_CONFIGFILE_WITH_LABELS_GOOD("./tests/config_makeimage_sersictest512_with_labels.dat");
const string  TEST_CONFIGFILE_BAD1("./tests/makeimage_reference/config_makeimage_sersictest512_bad1.dat");
const string  TEST_CONFIGFILE_BAD2("./tests/makeimage_reference/config_makeimage_sersictest512_bad2.dat");
const string  TEST_CONFIGFILE_BAD3("./tests/makeimage_reference/config_makeimage_sersictest512_bad3.dat");
const string  TEST_CONFIGFILE_BAD4("./tests/imfit_reference/config_imfit_badparamline.dat");  // parameter line with name only
const string  TEST_CONFIGFILE_BADLIMIT4("./tests/config_imfit_sersictest512_badlimits4.dat");  // parameter line with only one limit

// new stuff
const string TEST_CONFIGFILE_OPTIONAL0("./tests/config_makeimage-optional_params0.dat");
// single image function, one optional params
const string TEST_CONFIGFILE_OPTIONAL1("./tests/config_makeimage-optional_params1.dat");
// single image function, two optional params
const string TEST_CONFIGFILE_OPTIONAL2("./tests/config_makeimage-optional_params2.dat");
// three image functions: first has two optional params, second none, third one
const string TEST_CONFIGFILE_OPTIONAL3("./tests/config_makeimage-optional_params3.dat");
const string TEST_CONFIGFILE_OPTIONAL_BAD("./tests/config_makeimage-optional_params_bad1.dat");


class NewTestSuite : public CxxTest::TestSuite 
{
public:

  // Tests for ValidParameterLine()
  void testValidParameterLine_Good( void )
  {
    string  goodLine1 = "X0   32.0    1,64\n";
    string  goodLine2 = "X0   32.0\n";
    string  goodLine3 = "X0   32.0   # this is a comment\n";
    string  goodLine4 = "X0   32.0   30,35  # this is a comment\n";
    string  goodLine5 = "X0   32.0   fixed  # this is a comment\n";
    bool  result1, result2, result3, result4, result5;
    
    result1 = ValidParameterLine(goodLine1);
    TS_ASSERT_EQUALS(result1, true);
    result2 = ValidParameterLine(goodLine2);
    TS_ASSERT_EQUALS(result2, true);
    result3 = ValidParameterLine(goodLine3);
    TS_ASSERT_EQUALS(result3, true);
    result4 = ValidParameterLine(goodLine4);
    TS_ASSERT_EQUALS(result4, true);
    result5 = ValidParameterLine(goodLine5);
    TS_ASSERT_EQUALS(result5, true);
  }

  void testValidParameterLine_Bad( void )
  {
    string  badLine1 = "X0\n";
    string  badLine1b = "X0	\n";
    string  badLine1c = "X0   \n";
    string  badLine2 = "X0   fixed\n";
    string  badLine3 = "X0   33,36\n";
    string  badLine4a = "X0   32.0   30\n";
    string  badLine4b = "X0   32.0   30,  \n";
    string  badLine4c = "X0   32.0   x\n";
    string  badLine4d = "X0   32.0   30,x\n";
    bool  result1, result1b, result1c, result2, result3;
    bool  result4a, result4b, result4c, result4d;
    
    result1 = ValidParameterLine(badLine1);
    TS_ASSERT_EQUALS(result1, false);
    result1b = ValidParameterLine(badLine1b);
    TS_ASSERT_EQUALS(result1b, false);
    result1c = ValidParameterLine(badLine1c);
    TS_ASSERT_EQUALS(result1c, false);
    result2 = ValidParameterLine(badLine2);
    TS_ASSERT_EQUALS(result2, false);
    result3 = ValidParameterLine(badLine3);
    TS_ASSERT_EQUALS(result3, false);
    
    result4a = ValidParameterLine(badLine4a);
    TS_ASSERT_EQUALS(result4a, false);
    result4b = ValidParameterLine(badLine4b);
    TS_ASSERT_EQUALS(result4b, false);
    result4c = ValidParameterLine(badLine4c);
    TS_ASSERT_EQUALS(result4c, false);
    result4d = ValidParameterLine(badLine4d);
    TS_ASSERT_EQUALS(result4d, false);
  }

  void testValidOptionalParameterLine_Good( void )
  {
    string  goodLine1 = "mode  alpha\n";
    string  goodLine2 = "mode  alpha   # this is a comment\n";
    bool  result1, result2;
    
    result1 = ValidParameterLine(goodLine1, true);
    TS_ASSERT_EQUALS(result1, true);
    result2 = ValidParameterLine(goodLine2, true);
    TS_ASSERT_EQUALS(result2, true);
  }

  void testValidOptionalParameterLine_Bad( void )
  {
    string  badLine1 = "mode\n";
    string  badLine2 = "mode  # this is a comment\n";
    bool  result1, result2;
    
    result1 = ValidParameterLine(badLine1, true);
    TS_ASSERT_EQUALS(result1, false);
    result2 = ValidParameterLine(badLine2, true);
    TS_ASSERT_EQUALS(result2, false);
  }



  // Test for capture of function label
  void testGetFunctionLabel( void )
  {
    string  noLabelLine1 = "FUNCTION Sersic";
    string  noLabelLine2 = "FUNCTION Sersic   # this is not a label";
    string  labelLine1 = "FUNCTION Sersic # LABEL inner bulge";
    string  labelLine2 = "FUNCTION Sersic  LABEL inner bulge";
    string  result1, result2, result3, result4;
    string  correctNoLabel = "";
    string  correctLabel = "inner bulge";
    
    result1 = GetFunctionLabel(noLabelLine1);
    TS_ASSERT_EQUALS(result1, correctNoLabel);
    result2 = GetFunctionLabel(noLabelLine2);
    TS_ASSERT_EQUALS(result2, correctNoLabel);
    result3 = GetFunctionLabel(labelLine1);
    TS_ASSERT_EQUALS(result3, correctLabel);
    result4 = GetFunctionLabel(labelLine2);
    TS_ASSERT_EQUALS(result4, correctLabel);
  }



  // Tests for AddFunctionNameAndLabel
  void testAddFunctionNameAndLabel( void )
  {
    string  noLabelLine1 = "FUNCTION Sersic";
    string  noLabelLine2 = "FUNCTION Sersic   # this is not a label";
    string  labelLine1 = "FUNCTION Sersic # LABEL inner bulge";
    string  labelLine2 = "FUNCTION Sersic  LABEL inner bulge";
    
    vector<string>  funcNames;
    vector<string>  funcLabels;
    string  correctFunc = "Sersic";
    string  correctNoLabel = "";
    string  correctLabel = "inner bulge";
    
    funcNames.clear();
    funcLabels.clear();
    AddFunctionNameAndLabel(noLabelLine1, funcNames, funcLabels);
    TS_ASSERT_EQUALS(funcNames[0], correctFunc);
    TS_ASSERT_EQUALS(funcLabels[0], correctNoLabel);

    funcNames.clear();
    funcLabels.clear();
    AddFunctionNameAndLabel(noLabelLine2, funcNames, funcLabels);
    TS_ASSERT_EQUALS(funcNames[0], correctFunc);
    TS_ASSERT_EQUALS(funcLabels[0], correctNoLabel);

    funcNames.clear();
    funcLabels.clear();
    AddFunctionNameAndLabel(labelLine1, funcNames, funcLabels);
    TS_ASSERT_EQUALS(funcNames[0], correctFunc);
    TS_ASSERT_EQUALS(funcLabels[0], correctLabel);

    funcNames.clear();
    funcLabels.clear();
    AddFunctionNameAndLabel(labelLine2, funcNames, funcLabels);
    TS_ASSERT_EQUALS(funcNames[0], correctFunc);
    TS_ASSERT_EQUALS(funcLabels[0], correctLabel);
  }



  // Tests for VetConfigFile()
  void testVetConfigfile_Good( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_GOOD.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_GOOD.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = 2;
    output = VetConfigFile(inputLines, origLineNumbers, false, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
  }

  void testVetConfigfile_Good_with_labels( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_WITH_LABELS_GOOD.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_WITH_LABELS_GOOD.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = 2;
    output = VetConfigFile(inputLines, origLineNumbers, false, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
  }

  // optional params: empty optional-params specification
  void testVetConfigfile_Good_OptionalParams0( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_OPTIONAL0.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_OPTIONAL0.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = 2;
    output = VetConfigFile(inputLines, origLineNumbers, false, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
  }


  // optional params: valid format for param-name + value
  void testVetConfigfile_Good_OptionalParams1( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_OPTIONAL1.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_OPTIONAL1.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = 2;
    output = VetConfigFile(inputLines, origLineNumbers, false, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
  }


  // Test bad config file (no function section)
  void testVetConfigfile_Bad1( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD1.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD1.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_NOFUNCSECTION;
    output = VetConfigFile(inputLines, origLineNumbers, false, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
  }


  // Test bad config file (no functions listed)
  void testVetConfigfile_Bad2( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD2.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD2.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_NOFUNCTIONS;
    output = VetConfigFile(inputLines, origLineNumbers, false, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);
  }


  // Test bad config file (2D mode, but no Y0 to get with first X0)
  void testVetConfigfile_Bad3( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD3.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD3.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_INCOMPLETEXY;
    int  trueBadLineNumber = 8;  // X0 line in TEST_CONFIGFILE_BAD3
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);    
    TS_ASSERT_EQUALS(possibleBadLineNumber, trueBadLineNumber);
  }


  // Test bad config file (2D mode, mangled parameter lines)
  void testVetConfigfile_Bad4( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD4.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD4.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_BADPARAMLINE;
    int  trueBadLineNumber = 12;  // PA line in TEST_CONFIGFILE_BAD4
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);    
    TS_ASSERT_EQUALS(possibleBadLineNumber, trueBadLineNumber);
  }


  // Test bad config file (2D mode, mangled optional parameter lines)
  void testVetConfigfile_BadOptionalParam( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_OPTIONAL_BAD.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_OPTIONAL_BAD.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_BADPARAMLINE;
    int  trueBadLineNumber = 13;  // "mode" line in TEST_CONFIGFILE_OPTIONAL_BAD
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);    
    TS_ASSERT_EQUALS(possibleBadLineNumber, trueBadLineNumber);
  }


  // Test bad config file (2D mode, parameter line with only one limit)
  void testVetConfigfile_Bad_OneLimit( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    int  possibleBadLineNumber = -1;
    int  k = 0;
    vector<int>  origLineNumbers;
    
    inputFileStream.open(TEST_CONFIGFILE_BADLIMIT4.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BADLIMIT4.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      k++;
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
        inputLines.push_back(inputLine);
        origLineNumbers.push_back(k);
      }
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_BADPARAMLINE;
    int  trueBadLineNumber = 9;  // First PA line in TEST_CONFIGFILE_BADLIMIT4
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);    
    TS_ASSERT_EQUALS(possibleBadLineNumber, trueBadLineNumber);
  }



  // Tests for ParseFunctionSection()
  void testParseFunctionSection( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<mp_par>  paramLimits1;
    vector<int>  functionSetIndices1;
    bool  paramLimitsExist1;
	int  status;
	vector<string> inputLines;
	inputLines.push_back("X0    200.0\n");
	inputLines.push_back("Y0    200.0\n");
	inputLines.push_back("FUNCTION   Gaussian\n");
	inputLines.push_back("PA    10.0\n");
	inputLines.push_back("ell    0.1\n");
	inputLines.push_back("I_0    1000.0\n");
	inputLines.push_back("sigma      5.0\n");
	inputLines.push_back("X0    256.0\n");
	inputLines.push_back("Y0    256.0\n");
	inputLines.push_back("FUNCTION   Sersic_GenEllipse\n");
	inputLines.push_back("PA    30.0\n");
	inputLines.push_back("ell    0.5\n");
	inputLines.push_back("c0    -1.0\n");
	inputLines.push_back("n      2.5\n");
	inputLines.push_back("I_e  100.0\n");
	inputLines.push_back("r_e   50.0\n");
	inputLines.push_back("FUNCTION   FlatSky\n");
	inputLines.push_back("I_0  1.0\n");

	vector<int> originalLineNos = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
									15, 16, 17, 18};

	status = ParseFunctionSection(inputLines, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, functionSetIndices1, 
    							paramLimitsExist1, originalLineNos);
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 3);
    TS_ASSERT_EQUALS(functionList1[0], "Gaussian");
    TS_ASSERT_EQUALS(functionList1[1], "Sersic_GenEllipse");
    TS_ASSERT_EQUALS(functionList1[2], "FlatSky");
    
    int nParams = 15;
    double correctParamVals[15] = {200.0,200.0, 10.0, 0.1, 1000.0,5.0,
    								256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0, 1.0};
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)functionSetIndices1.size(), 2);
    TS_ASSERT_EQUALS(functionSetIndices1[0], 0);
    TS_ASSERT_EQUALS(functionSetIndices1[1], 1);
    TS_ASSERT_EQUALS(paramLimitsExist1, false);
  }

  void testParseFunctionSection_bad( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<mp_par>  paramLimits1;
    vector<int>  functionSetIndices1;
    bool  paramLimitsExist1;
	int  status;
	vector<string> inputLines1, inputLines2;
	
	// Missing "Y0" line
	inputLines1.push_back("X0    200.0\n");
	inputLines1.push_back("FUNCTION   Gaussian\n");
	inputLines1.push_back("PA    10.0\n");
	inputLines1.push_back("ell    0.1\n");
	inputLines1.push_back("I_0    1000.0\n");
	inputLines1.push_back("sigma      5.0\n");

	vector<int> originalLineNos1 = {1, 2, 3, 4, 5, 6};

	status = ParseFunctionSection(inputLines1, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, functionSetIndices1, 
    							paramLimitsExist1, originalLineNos1);
    TS_ASSERT_EQUALS(status, -1);

	// Bad parameter limits
	inputLines2.push_back("X0    200.0\n");
	inputLines2.push_back("Y0    200.0\n");
	inputLines2.push_back("FUNCTION   Gaussian\n");
	inputLines2.push_back("PA    10.0     20.0,1.0\n");
	inputLines2.push_back("ell    0.1\n");
	inputLines2.push_back("I_0    1000.0\n");
	inputLines2.push_back("sigma      5.0\n");

	vector<int> originalLineNos2 = {1, 2, 3, 4, 5, 6, 7};

	status = ParseFunctionSection(inputLines2, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, functionSetIndices1, 
    							paramLimitsExist1, originalLineNos2, true);
    TS_ASSERT_EQUALS(status, -1);
  }


// Full version, for use by e.g. imfit -- reads parameter limits as well
//    configFileName = C++ string with name of configuration file
//    mode2D = true for 2D functions (imfit, makeimage), false for 1D (imfit1d)
//    functionNameList = output, will contain vector of C++ strings containing functions
//                   names from config file
//    parameterList = output, will contain vector of parameter values
//    parameterLimits = output, will contain vector of mp_par structures (specifying
//                   possible limits on parameter values)
//    fsetStartIndices = output, will contain vector of integers specifying
//                   which functions mark start of new function set
// int ReadConfigFile( const string& configFileName, const bool mode2D, vector<string>& functionNameList,
//                     vector<string>& functionLabels, vector<double>& parameterList, 
//                     vector<mp_par>& parameterLimits, vector<int>& fsetStartIndices, 
//                     bool& parameterLimitsFound, configOptions& configFileOptions, 
//                     vector< map<string, string> >& optionalParamsVect=EMPTY_MAP_VECTOR_CONFIGPARSER );


  // Tests for ReadConfigFile() -- imfit version
  void testReadConfigFile_imfit( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<mp_par>  paramLimits1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
    bool  paramLimitsExist1;
	int  status;
	
    status = ReadConfigFile(TEST_CONFIGFILE_GOOD, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, FunctionSetIndices1, 
    							paramLimitsExist1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
    int nParams = 8;
    double correctParamVals[8] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0};
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    TS_ASSERT_EQUALS(paramLimitsExist1, false);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

    // do it again, to make sure the input/output vectors are cleared, then appended
    // to, rather than just appended to without being cleared first
    status = ReadConfigFile(TEST_CONFIGFILE_GOOD, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, FunctionSetIndices1, 
    							paramLimitsExist1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    TS_ASSERT_EQUALS(paramLimitsExist1, false);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

    // do it again, this time passing in an extraParamsMap input parameter,
    // which should come back empty (since config file has no extra
    // params)
    vector< map<string, string> >  optionalParamsVect;
    status = ReadConfigFile(TEST_CONFIGFILE_GOOD, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, FunctionSetIndices1, 
    							paramLimitsExist1, userConfigOptions1, optionalParamsVect);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    TS_ASSERT_EQUALS(paramLimitsExist1, false);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

    TS_ASSERT_EQUALS((int)optionalParamsVect.size(), 1);
  }

  // Same as previous, but now with function labels
  void testReadConfigFile_imfit_with_labels( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<mp_par>  paramLimits1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
    bool  paramLimitsExist1;
	bool  status;
	
    status = ReadConfigFile(TEST_CONFIGFILE_WITH_LABELS_GOOD, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, FunctionSetIndices1, 
    							paramLimitsExist1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    TS_ASSERT_EQUALS(functionLabels[0], "main component");
    
    int nParams = 8;
    double correctParamVals[8] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0};
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    TS_ASSERT_EQUALS(paramLimitsExist1, false);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

    // do it again, to make sure the input/output vectors are cleared, then appended
    // to, rather than just appended to without being cleared first
    status = ReadConfigFile(TEST_CONFIGFILE_WITH_LABELS_GOOD, true, functionList1, functionLabels,
    							parameterList1, paramLimits1, FunctionSetIndices1, 
    							paramLimitsExist1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    TS_ASSERT_EQUALS(functionLabels[0], "main component");
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    TS_ASSERT_EQUALS(paramLimitsExist1, false);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

  }

  // Tests for ReadConfigFile() -- makeimage version
  void testReadConfigFile_makeimage( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
	bool  status;
	
    status = ReadConfigFile(TEST_CONFIGFILE_GOOD, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
    int nParams = 8;
    double correctParamVals[8] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0};
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

    // do it again, to make sure the input/output vectors are cleared, then appended
    // to, rather than just appended to without being cleared first
    status = ReadConfigFile(TEST_CONFIGFILE_GOOD, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++)
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");
  }

  void testReadConfigFile_makeimage_with_labels( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
	bool  status;
	
    status = ReadConfigFile(TEST_CONFIGFILE_WITH_LABELS_GOOD, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    TS_ASSERT_EQUALS(functionLabels[0], "main component");
    
    int nParams = 8;
    double correctParamVals[8] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0};
    
    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++) {
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    }
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");

    // do it again, to make sure the input/output vectors are cleared, then appended
    // to, rather than just appended to without being cleared first
    status = ReadConfigFile(TEST_CONFIGFILE_WITH_LABELS_GOOD, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    TS_ASSERT_EQUALS(functionLabels[0], "main component");

    TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
    for (int i = 0; i < nParams; i++)
      TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
    
    TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
    TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
    TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
    TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");
  }

  void testReadConfigFile_makeimage_OptionalParams( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
    vector< map<string, string> >  optionalParamsVect;
    bool  status;

    status = ReadConfigFile(TEST_CONFIGFILE_OPTIONAL1, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1, 
    						optionalParamsVect);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    if (functionList1.size() >= 1) {
      TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
      int nParams = 8;
      double correctParamVals[8] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0};
    
      TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
      for (int i = 0; i < nParams; i++) {
        TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
      }
    
      TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
      TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
      TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
      TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
      TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
      TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
      TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");
    }
   
    // check that we have a vector with correct number of elements
    TS_ASSERT_EQUALS((int)optionalParamsVect.size(), 1);
    // extract the first (and only) map in the vector, check that it has a key named "mode"
    map<string, string> firstMap = optionalParamsVect[0];
    TS_ASSERT_EQUALS((int)firstMap.size(), 1);
    map<string, string>::iterator key_iterator;
    key_iterator = firstMap.find("mode");
    TS_ASSERT_DIFFERS(key_iterator, firstMap.end());
    // check that the value indexed by "mode" is "alpha"
    if (optionalParamsVect.size() >= 1) {
      map<string, string> firstMap = optionalParamsVect[0];
      string resultString = firstMap["mode"];
      TS_ASSERT_EQUALS(resultString, "alpha");
    }    
  }

  void testReadConfigFile_makeimage_OptionalParams_multiparams( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
    vector< map<string, string> >  optionalParamsVect;
    bool  status;

    status = ReadConfigFile(TEST_CONFIGFILE_OPTIONAL2, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1, 
    						optionalParamsVect);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 1);
    if (functionList1.size() >= 1) {
      TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
      int nParams = 8;
      double correctParamVals[8] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0};
    
      TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
      for (int i = 0; i < nParams; i++) {
        TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
      }
    
      TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
      TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
      TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
      TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
      TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
      TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
      TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");
    }
   
    // check that we have a vector with correct number of elements
    TS_ASSERT_EQUALS((int)optionalParamsVect.size(), 1);
    // check that optional-params map has correct number of entries
    map<string, string> optParamsMap = optionalParamsVect[0];
    TS_ASSERT_EQUALS((int)optParamsMap.size(), 2);

    // Check that it has a key named "mode"
    map<string, string>::iterator key_iterator;
    key_iterator = optParamsMap.find("mode");
    TS_ASSERT_DIFFERS(key_iterator, optParamsMap.end());
    // check that the value indexed by "mode" is "alpha"
    if (optionalParamsVect.size() >= 1) {
      string resultString = optParamsMap["mode"];
      TS_ASSERT_EQUALS(resultString, "alpha");
    }    
    // extract the second entry in the map, check that it has a key named "degree"
    key_iterator = optParamsMap.find("degree");
    TS_ASSERT_DIFFERS(key_iterator, optParamsMap.end());
    // check that the value indexed by "degree" is "12"
    if (optionalParamsVect.size() >= 1) {
      string resultString = optParamsMap["degree"];
      TS_ASSERT_EQUALS(resultString, "12");
    }    
  }

  void testReadConfigFile_makeimage_OptionalParams_multiparams_multfunc( void )
  {
    vector<string>  functionList1;
    vector<string>  functionLabels;
    vector<double>  parameterList1;
    vector<int>  FunctionSetIndices1;
    configOptions  userConfigOptions1;
    vector< map<string, string> >  optionalParamsVect;
    map<string, string> optParamsMap;
    map<string, string>::iterator key_iterator;
    string  resultString;
    bool  status;

    status = ReadConfigFile(TEST_CONFIGFILE_OPTIONAL3, true, functionList1, functionLabels,
    						parameterList1, FunctionSetIndices1, userConfigOptions1, 
    						optionalParamsVect);
  
    TS_ASSERT_EQUALS(status, 0);
    
    TS_ASSERT_EQUALS((int)functionList1.size(), 3);
    if (functionList1.size() >= 1) {
      TS_ASSERT_EQUALS(functionList1[0], "Sersic_GenEllipse");
    
      int nParams = 16;
      double correctParamVals[16] = {256.0,256.0, 30.0,0.5,-1.0,2.5,100.0,50.0,
      			20.0,0.0,200.0,35.0, 30.0,0.5,100.0,10.0};
    
      TS_ASSERT_EQUALS((int)parameterList1.size(), nParams);
      for (int i = 0; i < nParams; i++) {
        TS_ASSERT_EQUALS(parameterList1[i], correctParamVals[i]);
      }
    
      TS_ASSERT_EQUALS((int)FunctionSetIndices1.size(), 1);
      TS_ASSERT_EQUALS(FunctionSetIndices1[0], 0);
    
      TS_ASSERT_EQUALS(userConfigOptions1.nOptions, 2);
      TS_ASSERT_EQUALS(userConfigOptions1.optionNames[0], "NCOLS");
      TS_ASSERT_EQUALS(userConfigOptions1.optionNames[1], "NROWS");
      TS_ASSERT_EQUALS(userConfigOptions1.optionValues[0], "512");
      TS_ASSERT_EQUALS(userConfigOptions1.optionValues[1], "512");
    }
   
    // check that we have a vector with correct number of elements
    TS_ASSERT_EQUALS((int)optionalParamsVect.size(), 3);
    
    // check each optional-parameter map
    
    // first function (two optional parameters)
    optParamsMap = optionalParamsVect[0];
    TS_ASSERT_EQUALS((int)optParamsMap.size(), 2);
    key_iterator = optParamsMap.find("mode");
    TS_ASSERT_DIFFERS(key_iterator, optParamsMap.end());
    resultString = optParamsMap["mode"];
    TS_ASSERT_EQUALS(resultString, "alpha");
    key_iterator = optParamsMap.find("degree");
    TS_ASSERT_DIFFERS(key_iterator, optParamsMap.end());
    // check that the value indexed by "degree" is "12"
    resultString = optParamsMap["degree"];
    TS_ASSERT_EQUALS(resultString, "12");

    // second function (no optional parameters specified)
    optParamsMap = optionalParamsVect[1];
    TS_ASSERT_EQUALS((int)optParamsMap.size(), 0);

    // third function (one optional parameter)
    optParamsMap = optionalParamsVect[2];
    TS_ASSERT_EQUALS((int)optParamsMap.size(), 1);
    key_iterator = optParamsMap.find("floor");
    TS_ASSERT_DIFFERS(key_iterator, optParamsMap.end());
    resultString = optParamsMap["floor"];
    TS_ASSERT_EQUALS(resultString, "1.5");
  }

};
