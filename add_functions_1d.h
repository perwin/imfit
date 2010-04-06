/*   Public interfaces for function(s) which takes a list of user-specified
 * function objects and adds them to the ModelObject.
 */

#ifndef _ADD_FUNCTION_H_
#define _ADD_FUNCTION_H_

#include <string>
#include <vector>
#include "model_object.h"

using namespace std;


int AddFunctions1d( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices );

// Use the following to print out names of available functions/components
void PrintAvailableFunctions( );


#endif  // _ADD_FUNCTION_H_
