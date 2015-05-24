// Unit tests for code in config_file_parser.cpp

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

const string  TEST_CONFIGFILE_GOOD("tests/config_makeimage_sersictest512.dat");
const string  TEST_CONFIGFILE_BAD1("tests/config_makeimage_sersictest512_bad1.dat");
const string  TEST_CONFIGFILE_BAD2("tests/config_makeimage_sersictest512_bad2.dat");
const string  TEST_CONFIGFILE_BAD3("tests/config_makeimage_sersictest512_bad3.dat");
const string  TEST_CONFIGFILE_BAD4("tests/config_imfit_badparamline.dat");


class NewTestSuite : public CxxTest::TestSuite 
{
public:

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
    TS_ASSERT( output == correctOutput );
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT( output == correctOutput );
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
    TS_ASSERT( output == correctOutput );
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT( output == correctOutput );
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
    TS_ASSERT( output == correctOutput );    
    TS_ASSERT( possibleBadLineNumber == trueBadLineNumber );    
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
    int  trueBadLineNumber = 11;  // X0 line in TEST_CONFIGFILE_BAD3
    output = VetConfigFile(inputLines, origLineNumbers, true, &possibleBadLineNumber);
    TS_ASSERT_EQUALS(output, correctOutput);    
    TS_ASSERT_EQUALS(possibleBadLineNumber, trueBadLineNumber);    
  }


  // Tests for BadParameterLine()
  void testValidParameterLine_Good( void )
  {
    string  goodLine1 = "X0   32.0    1,64\n";
    string  goodLine2 = "X0   32.0\n";
    string  goodLine3 = "X0   32.0   # this is a comment\n";
    bool  result1, result2, result3;
    
    result1 = ValidParameterLine(goodLine1);
    TS_ASSERT_EQUALS(result1, true);
    result2 = ValidParameterLine(goodLine2);
    TS_ASSERT_EQUALS(result2, true);
    result3 = ValidParameterLine(goodLine3);
    TS_ASSERT_EQUALS(result3, true);
  }

  void testValidParameterLine_Bad( void )
  {
    string  badLine1 = "X0\n";
    string  badLine1b = "X0	\n";
    string  badLine1c = "X0   \n";
    string  badLine2 = "X0   fixed\n";
    string  badLine3 = "X0   33,36\n";
    bool  result1, result1b, result1c, result2, result3;
    
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
  }

};
