/* FILE: add_functions_1d.cpp -------------------------------------------- */
/*
 * Function which takes a vector of strings listing function names and generates
 * the corresponding FunctionObjects, passing them to the input ModelObject
 *
 * This version works with "1-D" functions instead of 2D image-oriented functions.
 *
 */

#include <string>
#include <vector>

#include "add_functions_1d.h"
#include "model_object.h"

#include "function_object.h"
#include "func1d_exp.h"
#include "func1d_gaussian.h"
#include "func1d_gaussian2side.h"
#include "func1d_moffat.h"
#include "func1d_sersic.h"
#include "func1d_core-sersic.h"
#include "func1d_broken-exp.h"
#include "func1d_delta.h"
#include "func1d_sech.h"
#include "func1d_sech2.h"
#include "func1d_vdksech.h"

using namespace std;


// CHANGE WHEN ADDING FUNCTION -- add function name to array, increment N_FUNCTIONS
const char  FUNCTION_NAMES[][30] = {"Exponential-1D", "Gaussian-1D", "Gaussian2Side-1D",
            "Moffat-1D", "Sersic-1D", "Core-Sersic-1D", "BrokenExponential-1D", 
            "Delta-1D", "Sech-1D", "Sech2-1D", "vdKSech-1D"};
const int  N_FUNCTIONS = 11;



int AddFunctions1d( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  
  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    printf("\tFunction: %s\n", currentName.c_str());
    
    if (currentName == "Exponential-1D") {
      thisFunctionObj = new Exponential1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Gaussian-1D") {
      thisFunctionObj = new Gaussian1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Gaussian2Side-1D") {
      thisFunctionObj = new Gaussian2Side1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Moffat-1D") {
      thisFunctionObj = new Moffat1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Sersic-1D") {
      thisFunctionObj = new Sersic1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Core-Sersic-1D") {
      thisFunctionObj = new CoreSersic1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "BrokenExponential-1D") {
      thisFunctionObj = new BrokenExponential1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Delta-1D") {
      thisFunctionObj = new Delta1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Sech-1D") {
      thisFunctionObj = new Sech1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Sech2-1D") {
      thisFunctionObj = new Sech21D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "vdKSech-1D") {
      thisFunctionObj = new vdKSech1D();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    // If we reach here, then something went wrong
    printf("*** AddFunctions1d: unidentified function name (\"%s\")\n", currentName.c_str());
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



void ListFunctionParameters( )
// Prints a list of function names, along with the ordered list of
// parameter names for each function (suitable for copying and pasting
// into a config file for makeimage or imfit).
{
  
  printf("\n<listing of function parameters not yet implemented!>\n\n");
}

/* END OF FILE: add_functions_1d.cpp ------------------------------------- */
