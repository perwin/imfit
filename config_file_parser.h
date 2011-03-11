// Experimental code for reading imfit parameter file
// Currently focused on getting the function names & associated parameters

#ifndef _CONFIG_FILE_PARSER_H_
#define _CONFIG_FILE_PARSER_H_

#include <vector>
#include <string>

#include "param_struct.h"

using namespace std;

// Error codes returned by VetConfigFile
#define CONFIG_FILE_ERROR_NOFUNCSECTION  -1
#define CONFIG_FILE_ERROR_NOFUNCTIONS  -2
#define CONFIG_FILE_ERROR_INCOMPLETEXY  -3

typedef struct {
  vector<string> optionNames;
  vector<string> optionValues;
  int  nOptions;
} configOptions;


// Utility function (only used inside ReadConfigFile, but exposed here so
// we can do unit tests on it
int VetConfigFile( vector<string>& inputLines, vector<int>& origLineNumbers, bool mode2D,
									int *badLineNumber );

// First version is for use by e.g. makeimage: reads in parameters, but ignores
// parameter limits
int ReadConfigFile( string& configFileName, bool mode2D, vector<string>& functionList,
                    vector<double>& parameterList, vector<int>& setStartFunctionNumber,
                     configOptions& configFileOptions );

// This version is for use by e.g. imfit: reads in parameters *and* parameter limits
int ReadConfigFile( string& configFileName, bool mode2D, vector<string>& functionList,
                    vector<double>& parameterList, vector<mp_par>& parameterLimits,
                    vector<int>& setStartFunctionNumber, bool& parameterLimitsFound,
                     configOptions& configFileOptions );


#endif  // _CONFIG_FILE_PARSER_H_
