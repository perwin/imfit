/* FILE: add_functions.cpp ----------------------------------------------- */
/*
 * Function which takes a vector of strings listing function names and generates
 * the corresponding FunctionObjects, passing them to the input ModelObject
 *
 * Places where you should insert/modify when adding a new function are indicated
 * by "CHANGE WHEN ADDING FUNCTION"
 *
 */

// Copyright 2010--2017 by Peter Erwin.
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


#include <string>
#include <vector>
#include <map>
#include <stdio.h>

#include "model_object.h"
#include "add_functions.h"

// CHANGE WHEN ADDING FUNCTION -- add corresponding header file
#include "function_objects/function_object.h"
#include "func_gaussian.h"
#include "func_sersic.h"
#include "func_exp.h"
#include "func_gen-exp.h"
#include "func_gen-sersic.h"
#include "func_core-sersic.h"
#include "func_broken-exp.h"
#include "func_broken-exp2d.h"
#include "func_edge-on-ring.h"
#include "func_edge-on-ring2side.h"
#include "func_gaussian-ring.h"
#include "func_gaussian-ring2side.h"
#include "func_moffat.h"
#include "func_king.h"
#include "func_king2.h"
#include "func_flatsky.h"
// modules requiring GSL:
#ifndef NO_GSL
#include "func_edge-on-disk.h"
#include "func_expdisk3d.h"
#include "func_brokenexpdisk3d.h"
#include "func_gaussianring3d.h"
#include "func_ferrersbar3d.h"
#include "func_pointsource.h"
#endif

// extra functions
#ifdef USE_EXTRA_FUNCS
#include "func_broken-exp-bar.h"
#include "func_brokenexpbar3d.h"
#include "func_boxytest3d.h"
#include "func_edge-on-disk_n4762.h"
#include "func_edge-on-disk_n4762v2.h"
#include "func_logspiral.h"
#include "func_logspiral_gauss.h"
#include "func_nan.h"
#include "func_expdisk3d_trunc.h"
#include "func_triaxbar3d.h"
#include "func_triaxbar3d_sq.h"
#include "func_triaxbar3d_gengauss_sq.h"
#include "func_exp-higher-mom.h"
#include "func_double-broken-exp.h"
#endif

// extra functions useful for e.g. unit tests
#ifdef USE_TEST_FUNCS
#include "func_gauss_extraparams.h"
#endif

using namespace std;



// Code to create FunctionObject object factories
// Abstract base class for FunctionObject factories

// Note that we need to declare and then define a virtual destructor for this
// class to avoid annoying (Clang) compiler warnings due to the "delete it->second"
// line in FreeFactories() --
// "warning: delete called on 'factory' that is abstract but has non-virtual destructor"
// (see http://stackoverflow.com/questions/10024796/c-virtual-functions-but-no-virtual-destructors)
class factory
{
public:
    virtual FunctionObject* create() = 0;
    virtual ~factory() = 0;
};

factory::~factory() {};


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

  ModifiedKing::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<ModifiedKing>();

  ModifiedKing2::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<ModifiedKing2>();

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

  FerrersBar3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<FerrersBar3D>();

  PointSource::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<PointSource>();
#endif

// extra functions
#ifdef USE_EXTRA_FUNCS
  BrokenExponentialBar::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExponentialBar>();

  BrokenExpBar3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BrokenExpBar3D>();

  BoxyTest3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<BoxyTest3D>();

  ExpDisk3D_PerfectTrunc::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<ExpDisk3D_PerfectTrunc>();

  // weird extra stuff we may not keep
  EdgeOnDiskN4762::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnDiskN4762>();

  EdgeOnDiskN4762v2::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<EdgeOnDiskN4762v2>();

  LogSpiral::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<LogSpiral>();

  LogSpiralGauss::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<LogSpiralGauss>();

  NaNFunc::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<NaNFunc>();

  TriaxBar3D::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<TriaxBar3D>();

  TriaxBar3D_SuperQuadric::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<TriaxBar3D_SuperQuadric>();

  TriaxBar3D_GenGaussian_SuperQuadric::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<TriaxBar3D_GenGaussian_SuperQuadric>();

  ExponentialHigherMoment::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<ExponentialHigherMoment>();

  DoubleBrokenExponential::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<DoubleBrokenExponential>();
#endif

#ifdef USE_TEST_FUNCS
  GaussianExtraParams::GetClassShortName(classFuncName);
  input_factory_map[classFuncName] = new funcobj_factory<GaussianExtraParams>();
#endif

}



int AddFunctions( ModelObject *theModel, const vector<string> &functionNameList,
                  vector<int> &functionBlockIndices, const bool subsamplingFlag, 
                  const int verboseLevel, vector< map<string, string> > &extraParams )
{
  int  nFunctions = functionNameList.size();
  int  status;
  string  currentName;
  bool  extraParamsExist = false;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;

  PopulateFactoryMap(factory_map);

  if (extraParams.size() > 0)
    extraParamsExist = true;
  
  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    if (verboseLevel >= 0)
      printf("Function: %s\n", currentName.c_str());
    if (factory_map.count(currentName) < 1) {
      fprintf(stderr, "*** AddFunctions: unidentified function name (\"%s\")\n", currentName.c_str());
      return -1;
    }
    else {
      thisFunctionObj = factory_map[currentName]->create();
      thisFunctionObj->SetSubsampling(subsamplingFlag);
      if (extraParamsExist) {
        // specialize the function as requested by user (via config file)
        if (verboseLevel >= 0)
          printf("   Setting optional parameter(s) for %s...\n", currentName.c_str());
        status = thisFunctionObj->SetExtraParams(extraParams[i]);
        if (status < 0) {
          fprintf(stderr, "Error attempting to set extra/optional parameters for ");
          fprintf(stderr, "function \"%s\"\n", thisFunctionObj->GetShortName().c_str());
        }
      }
      status = theModel->AddFunction(thisFunctionObj);
      if (status < 0) {
        fprintf(stderr, "Error attempting to add function \"%s\"", 
        		thisFunctionObj->GetShortName().c_str());
        fprintf(stderr, " to ModelObject!\n");
        return status;
      }
    }
  }
  // OK, we're done adding functions; now tell the model object to do some final setup
  // Tell model object about arrangement of functions into common-center sets
  theModel->DefineFunctionBlocks(functionBlockIndices);
  
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


void PrintAvailableFunctions( )
{
  vector<string>  functionNames;

  GetFunctionNames(functionNames);
  printf("\nAvailable function/components:\n\n");
  for (int i = 0; i < (int)functionNames.size(); i++)
    printf("%s\n", functionNames[i].c_str());
  printf("\n\n");    
}


// Prints a list of function names, along with the ordered list of
// parameter names for each function (suitable for copying and pasting
// into a config file for makeimage or imfit).
void ListFunctionParameters( )
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



// Gets the ordered list of parameter names for the specified function and returns them
// in the input vector parameterNameList
//
// Returns -1 if functionName is not the name of a valid Function Object class
// (Based on code from André Luiz de Amorim.)
int GetFunctionParameterNames(string &functionName, vector<string> &parameterNameList)
{
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;
//  vector<string> factory_map_names;

  PopulateFactoryMap(factory_map);

  if (factory_map.count(functionName) < 1) {
    return - 1;
  }
  else {
    thisFunctionObj = factory_map[functionName]->create();
    thisFunctionObj->GetParameterNames(parameterNameList);
    delete thisFunctionObj;
  }

  // Avoid minor memory leak by freeing the individual funcobj_factory objects
  FreeFactories(factory_map);

  return 0;
}


// Gets the list of names of known functions and returns them
// in the input vector functionNameList
//
// (Based on code from André Luiz de Amorim.)
void GetFunctionNames( vector<string> &functionNameList )
{
  string  currentName;
  FunctionObject  *thisFunctionObj;
  map<string, factory*>  factory_map;

  PopulateFactoryMap(factory_map);

  // get list of keys (function names) and step through it
  map<string, factory*>::iterator  w;

  for (w = factory_map.begin(); w != factory_map.end(); w++) {
    thisFunctionObj = w->second->create();
    currentName = thisFunctionObj->GetShortName();
    functionNameList.push_back(currentName);
    delete thisFunctionObj;
  }

  // Avoid minor memory leak by freeing the individual funcobj_factory objects
  FreeFactories(factory_map);
}



/* END OF FILE: add_functions.cpp ---------------------------------------- */
