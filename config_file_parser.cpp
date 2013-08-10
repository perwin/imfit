// Experimental code for reading imfit parameter file
// Currently focused on getting the function names & associated parameters

// NEXT STEP:
//    Modify (or clone & modify) AddParameter so that it includes code
// from lines 235--255 to generate & store a new mp_par structure
//

// Copyright 2010, 2011, 2012, 2013 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.



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
int AddParameterAndLimit( string& currentLine, vector<double>& parameterList,
														vector<mp_par>& parameterLimits, int origLineNumber );
void AddFunctionName( string& currentLine, vector<string>& functionList );
void ReportConfigError( int errorCode, int origLineNumber );


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
int AddParameterAndLimit( string& currentLine, vector<double>& parameterList,
														vector<mp_par>& parameterLimits, int origLineNumber ) {
  double  paramVal;
  double  lowerLimit, upperLimit;
  string  extraPiece;
  vector<string>  stringPieces, newPieces;
  mp_par  newParamLimit;
  int  paramLimitsFound = 0;
  
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
    paramLimitsFound = 1;
    extraPiece = stringPieces[2];
    //printf("Found a parameter limit: %s\n", extraPiece.c_str());
    if (extraPiece == fixedIndicatorString) {
      newParamLimit.fixed = 1;
    } else {
      if (extraPiece.find(',', 0) != string::npos) {
        newPieces.clear();
        SplitString(extraPiece, newPieces, ",");
        newParamLimit.limited[0] = 1;
        newParamLimit.limited[1] = 1;
        lowerLimit = strtod(newPieces[0].c_str(), NULL);
        upperLimit = strtod(newPieces[1].c_str(), NULL);
        if (lowerLimit >= upperLimit) {
          fprintf(stderr, "*** WARNING: first parameter limit for \"%s\" (%g) must be < second limit (%g)!\n",
          				stringPieces[0].c_str(), lowerLimit, upperLimit);
          if (lowerLimit == upperLimit) {
             fprintf(stderr, "    (To specify fixed parameters, use \"fixed\" keyword, not identical upper & lower parameter limits)\n");
          }
          fprintf(stderr, "    (Error on input line %d of configuration file)\n", origLineNumber);
          return -1;
        }
        if ((paramVal < lowerLimit) || (paramVal > upperLimit)) {
          fprintf(stderr, "*** WARNING: initial parameter value for \"%s\" (%g) must lie between the limits (%g,%g)!\n",
          				stringPieces[0].c_str(), paramVal, lowerLimit, upperLimit);
          fprintf(stderr, "    (Error on input line %d of configuration file)\n", origLineNumber);
          return -1;
        }
        newParamLimit.limits[0] = lowerLimit;
        newParamLimit.limits[1] = upperLimit;
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



/* ---------------- FUNCTION: VetConfigFile ---------------------------- */
// Utility function which checks a list of non-comment, non-empty lines from
// a config file for the following:
//    1. An X0 line (indicating the start of the function section)
//    2. At least one FUNCTION line after the X0 line
// Returns the index for the start of the function section on sucess, or
// error code on failure.  If an error can be traced to a particular line, then
// that line number is stored in *badLineNumber.
int VetConfigFile( vector<string>& inputLines, vector<int>& origLineNumbers, bool mode2D,
									int *badLineNumber )
{
  int  i, nInputLines;
  int  functionSectionStart = -1;
  bool  functionSectionFound = false;
  bool  allOK = false;
  bool  yValueOK = true;   // defaults to true for 1D case
  
  nInputLines = inputLines.size();
  // OK, locate the start of the function block (first line beginning with "X0")
  i = 0;
  while (i < nInputLines) {
    if (inputLines[i].find("X0", 0) != string::npos) {
      functionSectionStart = i;
      functionSectionFound = true;
      if (mode2D) {
        // for 2D (makeimage or imfit) mode, the next line must start with "Y0"
        if ( ((i + 1) == nInputLines) || (inputLines[i + 1].find("Y0", 0) == string::npos) ) {
          yValueOK = false;
          *badLineNumber = origLineNumbers[i];
        }
      }
      break;
    }
    i++;
  }
  
  if ((functionSectionFound) && (yValueOK)) {
    i = functionSectionStart;
    for (i = functionSectionStart; i < nInputLines; i++) {
      if (inputLines[i].find("FUNCTION", 0) != string::npos) {
        allOK = true;
        break;
      }
    }
  }

  if (allOK)
    return functionSectionStart;
  else {
    if (functionSectionFound) {
      if (yValueOK == false)
        return CONFIG_FILE_ERROR_INCOMPLETEXY;
      else
        return CONFIG_FILE_ERROR_NOFUNCTIONS;
    }
    else
      return CONFIG_FILE_ERROR_NOFUNCSECTION;
  }
}


/* ---------------- FUNCTION: ReportConfigError ------------------------ */
void ReportConfigError( int errorCode, int origLineNumber )
{
  switch (errorCode) {
    case CONFIG_FILE_ERROR_NOFUNCSECTION:
      fprintf(stderr, "\n*** ReadConfigFile: Unable to find start of function section in configuration file!");
      fprintf(stderr, " (no \"X0\" line found)\n");
      break;
    case CONFIG_FILE_ERROR_NOFUNCTIONS:
      fprintf(stderr, "\n*** ReadConfigFile: No actual functions found in configuration file!\n");
      break;
    case CONFIG_FILE_ERROR_INCOMPLETEXY:
      fprintf(stderr, "\n*** ReadConfigFile: Initial Y0 value (start of function section) missing in configuration file!\n");
      fprintf(stderr, "    (X0 specification on input line %d should be followed by Y0 specification on next line)\n", origLineNumber);
      break;
    default:
      fprintf(stderr, "\n*** ReadConfigFile: Unknown error code!\n");
      break;
  }
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
                     vector<double>& parameterList, vector<int>& setStartFunctionNumber,
                     configOptions& configFileOptions )
{
  ifstream  inputFileStream;
  string  inputLine, currentLine;
  vector<string>  inputLines;
  vector<string>  stringPieces;
  vector<int>  origLineNumbers;
  int  functionSectionStart, functionNumber;
  int  i, nInputLines;
  int  possibleBadLineNumber = -1;
  int  k = 0;
//  bool  functionSectionFound = false;
  
  inputFileStream.open(configFileName.c_str());
  if( ! inputFileStream ) {
     cerr << "Error opening input stream for file " << configFileName.c_str() << endl;
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

  nInputLines = inputLines.size();
  
  // OK, locate the start of the function block (first line beginning with "X0")
  functionSectionStart = VetConfigFile(inputLines, origLineNumbers, mode2D, &possibleBadLineNumber);
  if (functionSectionStart < 0) {
    ReportConfigError(functionSectionStart, possibleBadLineNumber);
    return -1;
  }

  // Parse the first (non-function-related) section here
  // We assume that each of lines has the form "CAPITAL_KEYWORD some_value"
  configFileOptions.nOptions = 0;
  configFileOptions.optionNames.clear();
  configFileOptions.optionValues.clear();
  for (i = 0; i < functionSectionStart; i++) {
    ChopComment(inputLines[i]);
    SplitString(inputLines[i], stringPieces);
    if (stringPieces.size() == 2) {
      configFileOptions.optionNames.push_back(stringPieces[0]);
      configFileOptions.optionValues.push_back(stringPieces[1]);
      configFileOptions.nOptions += 1;
    }
  }
  
  // OK, now parse the function section
  i = functionSectionStart;
  functionNumber = 0;
  while (i < nInputLines) {
    if (inputLines[i].find("X0", 0) != string::npos) {
      setStartFunctionNumber.push_back(functionNumber);
      AddParameter(inputLines[i], parameterList);
      i++;
      if (mode2D) {
        // X0 line should always be followed by Y0 line in 2D mode
        if (inputLines[i].find("Y0", 0) == string::npos) {
          fprintf(stderr, "*** WARNING: A 'Y0' line must follow each 'X0' line in the configuration file!\n");
          fprintf(stderr, "   (X0 specification on input line %d should be followed by Y0 specification on next line)\n",
          				origLineNumbers[i] - 1);
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
                    vector<int>& setStartFunctionNumber, bool& parameterLimitsFound,
                    configOptions& configFileOptions )
{
  ifstream  inputFileStream;
  string  inputLine, currentLine;
  vector<string>  inputLines;
  vector<string>  stringPieces;
  vector<int>  origLineNumbers;
  int  functionSectionStart, functionNumber, paramNumber;
  int  i, nInputLines;
  int  possibleBadLineNumber = -1;
  int  k = 0;
//  bool  functionSectionFound = false;
  int  pLimitFound;
  
  inputFileStream.open(configFileName.c_str());
  if( ! inputFileStream ) {
     cerr << "Error opening input stream for file " << configFileName.c_str() << endl;
  }
  while ( getline(inputFileStream, inputLine) ) {
    k++;
    // strip off leading & trailing spaces; turns a blank line with spaces/tabs
    // into an empty string
    TrimWhitespace(inputLine);
    if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
      inputLines.push_back(inputLine);
      origLineNumbers.push_back(k);
    }
  }
  inputFileStream.close();

  nInputLines = inputLines.size();

  // OK, locate the start of the function block (first line beginning with "X0")
  functionSectionStart = VetConfigFile(inputLines, origLineNumbers, mode2D, &possibleBadLineNumber);
  if (functionSectionStart < 0) {
    ReportConfigError(functionSectionStart, possibleBadLineNumber);
    return -1;
  }

  // Parse the first (non-function-related) section here
  // We assume that each of lines has the form "CAPITAL_KEYWORD some_value"
  configFileOptions.nOptions = 0;
  configFileOptions.optionNames.clear();
  configFileOptions.optionValues.clear();
  for (i = 0; i < functionSectionStart; i++) {
    ChopComment(inputLines[i]);
    SplitString(inputLines[i], stringPieces);
    if (stringPieces.size() == 2) {
      configFileOptions.optionNames.push_back(stringPieces[0]);
      configFileOptions.optionValues.push_back(stringPieces[1]);
      configFileOptions.nOptions += 1;
    }
  }
  
  // OK, now parse the function section
  i = functionSectionStart;
  functionNumber = 0;
  paramNumber = 0;
  parameterLimitsFound = false;
  while (i < nInputLines) {
    if (inputLines[i].find("X0", 0) != string::npos) {
      //printf("X0 detected (i = %d)\n", i);
      setStartFunctionNumber.push_back(functionNumber);
      pLimitFound = AddParameterAndLimit(inputLines[i], parameterList, parameterLimits,
      										origLineNumbers[i]);
      paramNumber++;
      if (pLimitFound < 0) {
        // Bad limit format or other problem -- bail out!
        return -1;
      }
      if (pLimitFound == 1)
        parameterLimitsFound = true;
      i++;
      if (mode2D) {
        // X0 line should always be followed by Y0 line in 2D mode
        if (inputLines[i].find("Y0", 0) == string::npos) {
          fprintf(stderr, "*** WARNING: A 'Y0' line must follow each 'X0' line in the configuration file!\n");
          fprintf(stderr, "   (X0 specification on input line %d should be followed by Y0 specification on next line)\n",
          				origLineNumbers[i] - 1);
          return -1;
        }
        pLimitFound = AddParameterAndLimit(inputLines[i], parameterList, parameterLimits,
        											origLineNumbers[i]);
        if (pLimitFound < 0) {
          // Bad limit format or other problem -- bail out!
          return -1;
        }
        if (pLimitFound == 1)
          parameterLimitsFound = true;
        paramNumber++;
        i++;
      }
      continue;
    }
    if (inputLines[i].find("FUNCTION", 0) != string::npos) {
      //printf("Function detected (i = %d)\n", i);
      AddFunctionName(inputLines[i], functionList);
      functionNumber++;
      i++;
      continue;
    }
    // OK, we only reach here if it's a regular (non-positional) parameter line
    //printf("Parameter detected (i = %d)\n", i);
    pLimitFound = AddParameterAndLimit(inputLines[i], parameterList, parameterLimits,
    										origLineNumbers[i]);
    if (pLimitFound < 0) {
      // Bad limit format or other problem -- bail out!
      return -1;
    }
    if (pLimitFound == 1)
      parameterLimitsFound = true;
    paramNumber++;
    i++;
  }
  
  return 0;
}

