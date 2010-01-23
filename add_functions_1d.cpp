/* FILE: add_functions.cpp ----------------------------------------------- */
/*
 * Function which takes a vector of strings listing function names and generates
 * the corresponding FunctionObjects, passing them to the input ModelObject
 *
 * This version works with "1-D" functions instead of 2D image-oriented functions.
 *
 */

#include <string>
#include <vector>

#include "model_object.h"

#include "function_object.h"
#include "func1d_exp.h"
#include "func1d_broken-exp.h"
#include "func1d_sersic.h"

using namespace std;


int AddFunctions( ModelObject *theModel, vector<string> &functionNameList )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  
  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    printf("Function: %s\n", currentName.c_str());
    if (currentName == "Exponential-1D") {
      thisFunctionObj = new Exponential1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Sersic-1D") {
      thisFunctionObj = new Sersic1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "BrokenExponential-1D") {
      thisFunctionObj = new BrokenExponential1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    // If we reach here, then something went wrong
    printf("*** AddFunctions: unidentified function name (\"%s\")\n", currentName.c_str());
    return - 1;
  }
  return 0;
}


/* END OF FILE: add_functions.cpp ---------------------------------------- */
