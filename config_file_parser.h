// Experimental code for reading imfit parameter file
// Currently focused on getting the function names & associated parameters

#ifndef _CONFIG_FILE_PARSER_H_
#define _CONFIG_FILE_PARSER_H_

#include <vector>
#include <string>

#include "param_struct.h"

using namespace std;


// First version is for use by e.g. makeimage: reads in parameters, but ignores
// parameter limits
int ReadConfigFile( string& configFileName, bool mode2D, vector<string>& functionList,
                    vector<double>& parameterList, vector<int>& setStartFunctionNumber );

// This version is for use by e.g. imfit: reads in parameters *and* parameter limits
int ReadConfigFile( string& configFileName, bool mode2D, vector<string>& functionList,
                    vector<double>& parameterList, vector<mp_par>& parameterLimits,
                    vector<int>& setStartFunctionNumber, bool& parameterLimitsFound );


#endif  // _CONFIG_FILE_PARSER_H_
