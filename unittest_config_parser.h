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
    
    inputFileStream.open(TEST_CONFIGFILE_GOOD.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_GOOD.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#'))
        inputLines.push_back(inputLine);
    }
    inputFileStream.close();

    correctOutput = 2;
    output = VetConfigFile(inputLines, false);
    TS_ASSERT( output == correctOutput );
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, true);
    TS_ASSERT( output == correctOutput );
  }

  // Test bad config file (no function section)
  void testVetConfigfile_Bad1( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD1.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD1.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#'))
        inputLines.push_back(inputLine);
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_NOFUNCSECTION;
    output = VetConfigFile(inputLines, false);
    TS_ASSERT( output == correctOutput );
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, true);
    TS_ASSERT( output == correctOutput );
  }

  // Test bad config file (no functions listed)
  void testVetConfigfile_Bad2( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD2.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD2.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#'))
        inputLines.push_back(inputLine);
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_NOFUNCTIONS;
    output = VetConfigFile(inputLines, false);
    TS_ASSERT( output == correctOutput );
    // check that this works in 2D mode also
    output = VetConfigFile(inputLines, true);
    TS_ASSERT( output == correctOutput );
  }

  // Test bad config file (2D mode, but no Y0 to get with first X0)
  void testVetConfigfile_Bad3( void )
  {
    ifstream  inputFileStream;
    string  inputLine;
    vector<string>  inputLines;
    int  output, correctOutput;
    
    inputFileStream.open(TEST_CONFIGFILE_BAD3.c_str());
    if( ! inputFileStream ) {
       cerr << "Error opening input stream for file " << TEST_CONFIGFILE_BAD3.c_str() << endl;
    }
    while ( getline(inputFileStream, inputLine) ) {
      // strip off leading & trailing spaces; turns a blank line with spaces/tabs
      // into an empty string
      TrimWhitespace(inputLine);
      // store non-empty, non-comment lines in a vector of strings
      if ((inputLine.size() > 0) && (inputLine[0] != '#'))
        inputLines.push_back(inputLine);
    }
    inputFileStream.close();

    correctOutput = CONFIG_FILE_ERROR_INCOMPLETEXY;
    output = VetConfigFile(inputLines, true);
    TS_ASSERT( output == correctOutput );    
  }

};
