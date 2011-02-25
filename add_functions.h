/*   Public interfaces for function(s) which takes a list of user-specified
 * function objects and adds them to the ModelObject.
 */

#ifndef _ADD_FUNCTION_H_
#define _ADD_FUNCTION_H_

#include <string>
#include <vector>
#include "model_object.h"

using namespace std;


int AddFunctions( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices, bool subamplingFlag );

// Use the following to print out names of available functions/components
void PrintAvailableFunctions( );

// Use the following to print out a full list consisting of each function
// name ("FUNCTION <short-name>") followed by the ordered list of parameter
// names (suitable for copying and pasting into a config file for makeimage or imfit).
void ListFunctionParameters( );


#endif  // _ADD_FUNCTION_H_
