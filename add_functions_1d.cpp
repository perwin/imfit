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
#include <map>
#include <stdio.h>

#include "add_functions_1d.h"
#include "model_object.h"

// CHANGE WHEN ADDING FUNCTION -- add corresponding header file
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

#include "func1d_n1543majmin_circbulge.h"
#include "func1d_n1543majmin.h"
#include "func1d_n1543majmin2.h"


using namespace std;


// CHANGE WHEN ADDING FUNCTION -- add function name to array, increment N_FUNCTIONS
// const char  FUNCTION_NAMES[][30] = {"Exponential-1D", "Gaussian-1D", "Gaussian2Side-1D",
//             "Moffat-1D", "Sersic-1D", "Core-Sersic-1D", "BrokenExponential-1D", 
//             "Delta-1D", "Sech-1D", "Sech2-1D", "vdKSech-1D", "n1543majmin-1D", "n1543majmin2-1D"};
// const int  N_FUNCTIONS = 12;


// Code to create FunctionObject object factories
// Abstract base class for FunctionObject factories
class factory
{
public:
    virtual FunctionObject* create() = 0;
};


// Template for derived FunctionObject factory classes
// (this implicitly sets up a whole set of derived classes, one for each
// FunctionOjbect class we substitute for the "function_object_type" placeholder)
template <class function_object_type>
class funcobj_factory : public factory
{
public:
   FunctionObject* create() { return new function_object_type(); }
};


// Miscellaneous function prototypes -- private to this module

void FreeFactories( map<string, factory*>& factory_map );



void PopulateFactoryMap( map<string, factory*>& input_factory_map )
{
  string  classFuncName;

  // CHANGE WHEN ADDING FUNCTION -- add new pair of lines for new function-object class
  // Here we create the map of function-object names (strings) and factory objects
  // (instances of the various template-specified factory subclasses)
  Exponential1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Exponential1D>();
  
  Gaussian1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Gaussian1D>();
  
  Gaussian2Side1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Gaussian2Side1D>();
  
  Moffat1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Moffat1D>();
  
  Sersic1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Sersic1D>();
  
  CoreSersic1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<CoreSersic1D>();
  
  BrokenExponential1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponential1D>();
  
  Delta1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Delta1D>();
  
  Sech1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Sech1D>();
  
  Sech21D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Sech21D>();
  
  vdKSech1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<vdKSech1D>();

  n1543MajMinCircBulge1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<n1543MajMinCircBulge1D>();

  n1543MajMin1D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<n1543MajMin1D>();

  n1543MajMin21D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<n1543MajMin21D>();
  
}



int AddFunctions1d( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;

  PopulateFactoryMap(factory_map);

  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    printf("\tFunction: %s\n", currentName.c_str());
    if (factory_map.count(currentName) < 1) {
      printf("*** AddFunctions: unidentified function name (\"%s\")\n", currentName.c_str());
      return - 1;
    }
    else {
      thisFunctionObj = factory_map[currentName]->create();
      theModel->AddFunction(thisFunctionObj);
    }
  }  
  // OK, we're done adding functions; now tell the model object to do some
  // final setup work
  // Tell model object about arrangement of functions into common-center sets
  theModel->DefineFunctionSets(functionSetIndices);
  
  // Tell model object to create vector of parameter labels
  theModel->PopulateParameterNames();

  // Avoid minor memory leak by freeing the individual funcobj_factory objects
  FreeFactories(factory_map);

  return 0;
}


// Function which frees the individual funcobj_factory objects inside the factory map
void FreeFactories( map<string, factory*>& factory_map )
{
  for (map<string, factory*>::iterator it = factory_map.begin(); it != factory_map.end(); ++it)
    delete it->second;
}



// void PrintAvailableFunctions( )
// {
//   
//   printf("\nAvailable function/components:\n");
//   for (int i = 0; i < N_FUNCTIONS - 1; i++) {
//     printf("%s, ", FUNCTION_NAMES[i]);
//   }
//   printf("%s.\n\n", FUNCTION_NAMES[N_FUNCTIONS - 1]);
//     
// }



void PrintAvailableFunctions( )
{
  string  currentName;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;

  PopulateFactoryMap(factory_map);

  // get list of keys (function names) and step through it
  map<string, factory*>::iterator  w;

  printf("\nAvailable function/components:\n\n");
  for (w = factory_map.begin(); w != factory_map.end(); w++) {
//    printf("%s, ", w->first.c_str());
    thisFunctionObj = w->second->create();
    currentName = thisFunctionObj->GetShortName();
    printf("%s\n", currentName.c_str());
    delete thisFunctionObj;
  }
  printf("\n\n");    

  // Avoid minor memory leak by freeing the individual funcobj_factory objects
  FreeFactories(factory_map);
}


void ListFunctionParameters( )
// Prints a list of function names, along with the ordered list of
// parameter names for each function (suitable for copying and pasting
// into a config file for makeimage or imfit).
{
  
  string  currentName;
  vector<string>  parameterNameList;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;

  PopulateFactoryMap(factory_map);

  // get list of keys (function names) and step through it
  map<string, factory*>::iterator  w;

  printf("\nAvailable function/components:\n");
  for (w = factory_map.begin(); w != factory_map.end(); w++) {
//    printf("%s, ", w->first.c_str());
    thisFunctionObj = w->second->create();
    currentName = thisFunctionObj->GetShortName();
    printf("\nFUNCTION %s\n", currentName.c_str());
    parameterNameList.clear();
    thisFunctionObj->GetParameterNames(parameterNameList);
    for (int i = 0; i < (int)parameterNameList.size(); i++)
      printf("%s\n", parameterNameList[i].c_str());
    delete thisFunctionObj;
  }
  printf("\n\n");
    
  // Avoid minor memory leak by freeing the individual funcobj_factory objects
  FreeFactories(factory_map);
}


/* END OF FILE: add_functions_1d.cpp ------------------------------------- */
