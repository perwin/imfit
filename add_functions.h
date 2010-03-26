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


#endif  // _ADD_FUNCTION_H_
