/** @file
    \brief Public interfaces for code which parses imfit/makeimage config files

 */
// Currently focused on getting the function names & associated parameters

#ifndef _CONFIG_FILE_PARSER_H_
#define _CONFIG_FILE_PARSER_H_

#include <vector>
#include <string>

#include "param_struct.h"

using namespace std;

// Error codes returned by VetConfigFile
const int CONFIG_FILE_ERROR_NOFUNCSECTION = -1;
const int CONFIG_FILE_ERROR_NOFUNCTIONS   = -2;
const int CONFIG_FILE_ERROR_INCOMPLETEXY  = -3;
const int CONFIG_FILE_ERROR_BADPARAMLINE  = -4;

typedef struct {
  vector<string> optionNames;
  vector<string> optionValues;
  int  nOptions;
} configOptions;


/// \brief Utility function which does basic sanity-checking on config file
//
/// (This is only used inside ReadConfigFile, but is exposed in the header file so
/// we can do unit tests on it)
int VetConfigFile( vector<string>& inputLines, const vector<int>& origLineNumbers, 
					const bool mode2D, int *badLineNumber );

/// \brief Utility function: returns true if line has 2+ elements and 2nd is numeric
//
/// (This is only used inside VetConfigFile, but is exposed in the header file so
/// we can do unit tests on it)
bool ValidParameterLine( string& currentLine );

/// Function for use by makeimage
int ReadConfigFile( const string& configFileName, const bool mode2D, vector<string>& functionNameList,
                    vector<double>& parameterList, vector<int>& fblockStartIndices,
                     configOptions& configFileOptions );

/// Function for use by e.g. imfit: reads in parameters *and* parameter limits
int ReadConfigFile( const string& configFileName, const bool mode2D, vector<string>& functionNameList,
                    vector<double>& parameterList, vector<mp_par>& parameterLimits,
                    vector<int>& fblockStartIndices, bool& parameterLimitsFound,
                     configOptions& configFileOptions );


#endif  // _CONFIG_FILE_PARSER_H_
