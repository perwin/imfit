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
#include <map>
#include <stdio.h>

#include "model_object.h"

// CHANGE WHEN ADDING FUNCTION -- add corresponding header file
#include "function_object.h"
#include "func_gaussian.h"
#include "func_sersic.h"
#include "func_exp.h"
#include "func_gen-exp.h"
#include "func_gen-sersic.h"
#include "func_flat-exp.h"
#include "func_broken-exp.h"
#include "func_broken-exp2d.h"
#include "func_moffat.h"
#include "func_flatsky.h"

using namespace std;


// CHANGE WHEN ADDING FUNCTION -- add function name to array, increment N_FUNCTIONS
const char  FUNCTION_NAMES[][30] = {"Exponential", "Exponential_GenEllipse", "Sersic", 
            "Sersic_GenEllipse", "Gaussian", "FlatExponential", "BrokenExponential", 
            "BrokenExponential2D", "Moffat", "FlatSky"};
const int  N_FUNCTIONS = 10;



// Code to create FunctionObject object factories
class factory
{
public:
    virtual FunctionObject* create() = 0;
};


template <class function_object_type>
class funcobj_factory : public factory
{
public:
   FunctionObject* create() { return new function_object_type(); }
};



void PopulateFactoryMap( map<string, factory*>& input_factory_map )
{
  string  classFuncName;

  // CHANGE WHEN ADDING FUNCTION -- add new pair of lines for new function-object class
  // Here we create the map of function-object names and factory objects
  Exponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Exponential>();
  
  Sersic::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Sersic>();
  
  GenSersic::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GenSersic>();
  
  GenExponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GenExponential>();
  
  BrokenExponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponential>();
  
  FlatExponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<FlatExponential>();
  
  BrokenExponential2D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponential2D>();
  
  Gaussian::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Gaussian>();

  Moffat::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Moffat>();

  FlatSky::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<FlatSky>();
}




int AddFunctions( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices, bool subsamplingFlag )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;

  PopulateFactoryMap(factory_map);

  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    printf("Function: %s\n", currentName.c_str());
    if (factory_map.count(currentName) < 1) {
      printf("*** AddFunctions: unidentified function name (\"%s\")\n", currentName.c_str());
      return - 1;
    }
    else {
      thisFunctionObj = factory_map[currentName]->create();
      thisFunctionObj->SetSubsampling(subsamplingFlag);
      theModel->AddFunction(thisFunctionObj);
    }
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
