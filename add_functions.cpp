/* FILE: add_functions.cpp ----------------------------------------------- */
/*
 * Function which takes a vector of strings listing function names and generates
 * the corresponding FunctionObjects, passing them to the input ModelObject
 *
 * Places where you should insert/modify when adding a new function are indicated
 * by "CHANGE WHEN ADDING FUNCTION"
 *
 */

// Copyright 2010, 2011, 2012, 2013 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.


#ifndef USING_SCONS
#include "config.h"
#endif

#include <string>
#include <vector>
#include <map>
#include <stdio.h>

#include "model_object.h"

// CHANGE WHEN ADDING FUNCTION -- add corresponding header file
#include "function_objects/function_object.h"
#include "func_gaussian.h"
#include "func_sersic.h"
#include "func_exp.h"
#include "func_gen-exp.h"
#include "func_gen-sersic.h"
#include "func_core-sersic.h"
//#include "func_flat-exp.h"
#include "func_broken-exp.h"
#include "func_broken-exp2d.h"
#include "func_edge-on-ring.h"
#include "func_edge-on-ring2side.h"
#include "func_gaussian-ring.h"
#include "func_gaussian-ring2side.h"
#include "func_moffat.h"
#include "func_flatsky.h"
// modules requiring GSL:
#ifndef NO_GSL
#include "func_edge-on-disk.h"
#include "func_expdisk3d.h"
#include "func_brokenexpdisk3d.h"
#include "func_gaussianring3d.h"
#endif

// extra functions
#ifdef USE_EXTRA_FUNCS
//#include "func_expdisk3dv2.h"
#include "func_broken-exp-bar.h"
#include "func_brokenexpbar3d.h"
#include "func_boxytest3d.h"
#include "func_edge-on-disk_n4762.h"
#include "func_edge-on-disk_n4762v2.h"
#endif

using namespace std;


// CHANGE WHEN ADDING FUNCTION -- add function name to array, increment N_FUNCTIONS
// FUNCTION_NAMES is only used by PrintAvailableFunctions
// #ifndef NO_GSL
// const char  FUNCTION_NAMES[][30] = {"Exponential", "Exponential_GenEllipse", "Sersic", 
//             "Sersic_GenEllipse", "Core-Sersic", "Gaussian", "BrokenExponential", 
//             "BrokenExponentialBar", "BrokenExponential2D", "EdgeOnDisk", "Moffat", "FlatSky",
//             "EdgeOnDiskN4762", "EdgeOnDiskN4762v2", "EdgeOnRing", "EdgeOnRing2side",
//             "GaussianRing", "GaussianRing2Side", "ExponentialDisk3D", "GaussianRing3D"};
// const int  N_FUNCTIONS = 20;
// #else
// const char  FUNCTION_NAMES[][30] = {"Exponential", "Exponential_GenEllipse", "Sersic", 
//             "Sersic_GenEllipse", "Core-Sersic", "Gaussian", "BrokenExponential", 
//             "BrokenExponentialBar", "BrokenExponential2D", "Moffat", "FlatSky",
//             "EdgeOnDiskN4762", "EdgeOnDiskN4762v2", "EdgeOnRing", "EdgeOnRing2side",
//             "GaussianRing", "GaussianRing2Side"};
// const int  N_FUNCTIONS = 17;
// #endif


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
  Exponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Exponential>();
  
  Sersic::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Sersic>();
  
  GenSersic::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GenSersic>();
  
  CoreSersic::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<CoreSersic>();
  
  GenExponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GenExponential>();
  
  BrokenExponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponential>();
  
  BrokenExponential2D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponential2D>();
  
//   FlatExponential::GetClassShortName(classFuncName);
//   input_factory_map[classFuncName] = new funcobj_factory<FlatExponential>();
  
  EdgeOnRing::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnRing>();
  
  EdgeOnRing2Side::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnRing2Side>();
  
  GaussianRing::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GaussianRing>();
  
  GaussianRing2Side::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GaussianRing2Side>();
  
  Gaussian::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Gaussian>();

  Moffat::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<Moffat>();

  FlatSky::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<FlatSky>();

// functions requring GSL:
#ifndef NO_GSL 
  EdgeOnDisk::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnDisk>();

  ExponentialDisk3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<ExponentialDisk3D>();

  BrokenExponentialDisk3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponentialDisk3D>();

  GaussianRing3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GaussianRing3D>();
#endif

// extra functions
#ifdef USE_EXTRA_FUNCS
  // in testing/development
//  ExponentialDisk3Dv2::GetClassShortName(classFuncName);
//  input_factory_map[classFuncName] = new funcobj_factory<ExponentialDisk3Dv2>();

  BrokenExponentialBar::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponentialBar>();

  BrokenExpBar3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExpBar3D>();

  BoxyTest3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BoxyTest3D>();

  // weird extra stuff we may not keep
  EdgeOnDiskN4762::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnDiskN4762>();

  EdgeOnDiskN4762v2::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnDiskN4762v2>();

#endif
}




int AddFunctions( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices, bool subsamplingFlag )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;
  vector<string> factory_map_names;

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
  // OK, we're done adding functions; now tell the model object to do some final setup
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


/* END OF FILE: add_functions.cpp ---------------------------------------- */
