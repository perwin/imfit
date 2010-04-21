/* FILE: add_functions.cpp ----------------------------------------------- */
/*
 * Function which takes a vector of strings listing function names and generates
 * the corresponding FunctionObjects, passing them to the input ModelObject
 *
 * Places where you should insert/modify when adding a new function are indicated
 * by "CHANGE WHEN ADDING FUNCTION"
 *
 */

#include <string>
#include <vector>
#include <stdio.h>

#include "model_object.h"

// CHANGE WHEN ADDING FUNCTION -- add corresponding header file
#include "function_object.h"
#include "func_gaussian.h"
#include "func_sersic.h"
#include "func_exp.h"
#include "func_flat-exp.h"
#include "func_broken-exp.h"
#include "func_moffat.h"

using namespace std;


// CHANGE WHEN ADDING FUNCTION -- add function name to array, increment N_FUNCTIONS
const char  FUNCTION_NAMES[][30] = {"Exponential", "Sersic", "Gaussian", 
            "FlatExponential", "BrokenExponential", "Moffat"};
const int  N_FUNCTIONS = 6;


// CHANGE WHEN ADDING FUNCTION -- add corresponding
// 'if (currentName == "<function name>" {}' block
int AddFunctions( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices, bool subamplingFlag )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  
  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    printf("Function: %s\n", currentName.c_str());
    
    if (currentName == "Exponential") {
      thisFunctionObj = new Exponential(subamplingFlag);
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Sersic") {
      thisFunctionObj = new Sersic(subamplingFlag);
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Gaussian") {
      thisFunctionObj = new Gaussian(subamplingFlag);
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "FlatExponential") {
      thisFunctionObj = new FlatExponential(subamplingFlag);
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "BrokenExponential") {
      thisFunctionObj = new BrokenExponential(subamplingFlag);
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Moffat") {
      thisFunctionObj = new Moffat(subamplingFlag);
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    // If we reach here, then something went wrong
    printf("*** AddFunctions: unidentified function name (\"%s\")\n", currentName.c_str());
    return - 1;
  }
  
  // OK, we're done adding functions; now tell the model object to do some
  // final setup work
  // Tell model object about arrangement of functions into common-center sets
  theModel->DefineFunctionSets(functionSetIndices);
  
  // Tell model object to create vector of parameter labels
  theModel->PopulateParameterNames();
  return 0;
}



void PrintAvailableFunctions( )
{
  
  printf("\nAvailable function/components:\n");
  for (int i = 0; i < N_FUNCTIONS - 1; i++) {
    printf("%s, ", FUNCTION_NAMES[i]);
  }
  printf("%s.\n\n", FUNCTION_NAMES[N_FUNCTIONS - 1]);
    
}


/* END OF FILE: add_functions.cpp ---------------------------------------- */
