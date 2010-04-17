// Experimental code for reading imfit parameter file
// Currently focused on getting the function names & associated parameters

// NEXT STEP:
//    Modify (or clone & modify) AddParameter so that it includes code
// from lines 235--255 to generate & store a new mp_par structure
//

/* ------------------------ Include Files (Header Files )--------------- */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "param_struct.h"
#include "utilities_pub.h"
#include "config_file_parser.h"

using namespace std;


/* ------------------- Function Prototypes ----------------------------- */
void AddParameter( string& currentLine, vector<double>& parameterList );
bool AddParameterAndLimit( string& currentLine, vector<double>& parameterList,
														vector<mp_par>& parameterLimits );
void AddFunctionName( string& currentLine, vector<string>& functionList );


/* ------------------------ Global Variables --------------------------- */

static string  fixedIndicatorString = "fixed";



/* ---------------- FUNCTION: AddParameter ----------------------------- */
// Parses a line, extracting the second element as a floating-point value and
// storing it in the parameterList vector.
void AddParameter( string& currentLine, vector<double>& parameterList ) {
  double  paramVal;
  vector<string>  stringPieces;
  
  ChopComment(currentLine);
  stringPieces.clear();
  SplitString(currentLine, stringPieces);
  // first piece is parameter name, which we ignore; second piece is initial value
  paramVal = strtod(stringPieces[1].c_str(), NULL);
  parameterList.push_back(paramVal);
}


/* ---------------- FUNCTION: AddParameterAndLimit --------------------- */
// Parses a line, extracting the second element as a floating-point value and
// storing it in the parameterList vector.  In addition, if parameter limits
// are present (3rd element of line), they are extracted and an mp_par structure
// is added to the parameterLimits vector.
// Returns true for the existence of a parameter limit; false if no limits were
// found
bool AddParameterAndLimit( string& currentLine, vector<double>& parameterList,
														vector<mp_par>& parameterLimits ) {
  double  paramVal;
  string  extraPiece;
  vector<string>  stringPieces, newPieces;
  mp_par  newParamLimit;
  bool  paramLimitsFound = false;
  
  ChopComment(currentLine);
  stringPieces.clear();
  SplitString(currentLine, stringPieces);
  // first piece is parameter name, which we ignore; second piece is initial value
  paramVal = strtod(stringPieces[1].c_str(), NULL);
  parameterList.push_back(paramVal);

  // OK, now we create a new mp_par structure and check for possible parameter limits
  bzero(&newParamLimit, sizeof(mp_par));
  if (stringPieces.size() > 2) {
    // parse and store parameter limits, if any
    paramLimitsFound = true;
    extraPiece = stringPieces[2];
    printf("Found a parameter limit: %s\n", extraPiece.c_str());
    if (extraPiece == fixedIndicatorString) {
      newParamLimit.fixed = 1;
    } else {
      if (extraPiece.find(',', 0) != string::npos) {
        newPieces.clear();
        SplitString(extraPiece, newPieces, ",");
        newParamLimit.limited[0] = 1;
        newParamLimit.limited[1] = 1;
        newParamLimit.limits[0] = strtod(newPieces[0].c_str(), NULL);
        newParamLimit.limits[1] = strtod(newPieces[1].c_str(), NULL);
      }
    }
  }
  parameterLimits.push_back(newParamLimit);
  
  return paramLimitsFound;
}


/* ---------------- FUNCTION: AddFunctionName -------------------------- */
// Parses a line, extracting the second element as a string and storing it in 
// the functionList vector.
void AddFunctionName( string& currentLine, vector<string>& functionList ) {
  vector<string>  stringPieces;
  
  ChopComment(currentLine);
  stringPieces.clear();
  SplitString(currentLine, stringPieces);
  // store the actual function name (remember that first token is "FUNCTION")
  functionList.push_back(stringPieces[1]);
}



/* ---------------- FUNCTION: ReadConfigFile --------------------------- */
// Limited version, for use by e.g. makeimage -- ignores parameter limits!
//    configFileName = C++ string with name of configuration file
//    mode2D = true for 2D functions (imfit, makeimage), false for 1D (imfit1d)
//    functionList = output, will contain vector of C++ strings containing functions
//                   names from config file
//    parameterList = output, will contain vector of parameter values
//    setStartFunctionNumber = output, will contain vector of integers specifying
//                   which functions mark start of new function set
int ReadConfigFile( string& configFileName, bool mode2D, vector<string>& functionList,
                     vector<double>& parameterList, vector<int>& setStartFunctionNumber )
{
  ifstream  inputFileStream;
  string  inputLine, currentLine;
  vector<string>  inputLines;
  vector<string>  stringPieces;
  int  functionNumber;
  int  i, nInputLines;
  
  inputFileStream.open(configFileName.c_str());
  if( ! inputFileStream ) {
     cerr << "Error opening input stream for file " << configFileName.c_str() << endl;
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

  nInputLines = inputLines.size();
  i = 0;
  functionNumber = 0;
  while (i < nInputLines) {
    if (inputLines[i].find("X0", 0) != string::npos) {
      setStartFunctionNumber.push_back(functionNumber);
      AddParameter(inputLines[i], parameterList);
      i++;
      if (mode2D) {
        // X0 line should always be followed by Y0 line in 2D mode
        if (inputLines[i].find("Y0", 0) == string::npos) {
          printf("*** WARNING: A 'Y0' line must follow each 'X0' line in the configuration file!\n");
          return -1;
        }
        AddParameter(inputLines[i], parameterList);
        i++;
      }
      continue;
    }
    if (inputLines[i].find("FUNCTION", 0) != string::npos) {
      AddFunctionName(inputLines[i], functionList);
      functionNumber++;
      i++;
      continue;
    }
    // OK, we only reach here if it's a regular (non-positional) parameter line
    AddParameter(inputLines[i], parameterList);
    i++;
  }
  
  return 0;
}



/* ---------------- FUNCTION: ReadConfigFile --------------------------- */
// Full version, for use by e.g. imfit -- reads parameter limits as well
//    configFileName = C++ string with name of configuration file
//    mode2D = true for 2D functions (imfit, makeimage), false for 1D (imfit1d)
//    functionList = output, will contain vector of C++ strings containing functions
//                   names from config file
//    parameterList = output, will contain vector of parameter values
//    parameterLimits = output, will contain vector of mp_par structures (specifying
//                   possible limits on parameter values)
//    setStartFunctionNumber = output, will contain vector of integers specifying
//                   which functions mark start of new function set
int ReadConfigFile( string& configFileName, bool mode2D, vector<string>& functionList,
                    vector<double>& parameterList, vector<mp_par>& parameterLimits,
                    vector<int>& setStartFunctionNumber, bool& parameterLimitsFound )
{
  ifstream  inputFileStream;
  string  inputLine, currentLine;
  vector<string>  inputLines;
  vector<string>  stringPieces;
  int  functionNumber, paramNumber;
  int  i, nInputLines;
  bool  pLimitFound;
  
  inputFileStream.open(configFileName.c_str());
  if( ! inputFileStream ) {
     cerr << "Error opening input stream for file " << configFileName.c_str() << endl;
  }
  while ( getline(inputFileStream, inputLine) ) {
    // strip off leading & trailing spaces; turns a blank line with spaces/tabs
    // into an empty string
    TrimWhitespace(inputLine);
    if ((inputLine.size() > 0) && (inputLine[0] != '#'))
      inputLines.push_back(inputLine);
  }
  inputFileStream.close();

  nInputLines = inputLines.size();
  i = 0;
  functionNumber = 0;
  paramNumber = 0;
  parameterLimitsFound = false;
  while (i < nInputLines) {
    if (inputLines[i].find("X0", 0) != string::npos) {
      //printf("X0 detected (i = %d)\n", i);
      setStartFunctionNumber.push_back(functionNumber);
      pLimitFound = AddParameterAndLimit(inputLines[i], parameterList, parameterLimits);
      paramNumber++;
      if (pLimitFound)
        parameterLimitsFound = true;
      i++;
      if (mode2D) {
        // X0 line should always be followed by Y0 line in 2D mode
        if (inputLines[i].find("Y0", 0) == string::npos) {
          printf("*** WARNING: A 'Y0' line must follow each 'X0' line in the configuration file!\n");
          return -1;
        }
        AddParameterAndLimit(inputLines[i], parameterList, parameterLimits);
        paramNumber++;
        i++;
      }
      continue;
    }
    if (inputLines[i].find("FUNCTION", 0) != string::npos) {
      printf("Function detected (i = %d)\n", i);
      AddFunctionName(inputLines[i], functionList);
      functionNumber++;
      i++;
      continue;
    }
    // OK, we only reach here if it's a regular (non-positional) parameter line
    //printf("Parameter detected (i = %d)\n", i);
    AddParameterAndLimit(inputLines[i], parameterList, parameterLimits);
    paramNumber++;
    i++;
  }
  
  return 0;
}

