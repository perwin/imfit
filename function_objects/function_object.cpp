/* FILE: function_object.cpp ------------------------------------------- */
/* 
 *   This is the base class for the various function object classes.
 *   It really shouldn't be instantiated by itself.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *      So for a vector of x[], we expect the user to call Setup() once at
 *      the start, then loop through x[], calling GetValue() to compute
 *      the function results for each x.
 *
 *   MODIFICATION HISTORY:
 *     [v0.25]: 5 Dec: Minor cleanup.
 *     [v0.2]: 26---27 Nov: Modifications to accomodate derived classes
 * for working with 1-D fitting; also changed uses of char * to uses of
 * C++ string class.
 *     [v0.01]: 13 Nov 2009: Created (as modification of nonlinfit2's
 *   function_object class.
 */

// Copyright 2010--2022 by Peter Erwin.
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



/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <string.h>
#include <map>
#include <string>

#include "function_object.h"
#include "utilities_pub.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */


/* ---------------- CONSTRUCTOR ---------------------------------------- */

/// Basic constructor: sets functionName and shortFunctionName
FunctionObject::FunctionObject( )
{
  // These will get redefined by derived class's constructor
  functionName = "Base (undefined) function";
  shortFunctionName = "BaseFunction";
  extraParamsSet = false;
  
  // MULTIMFIT default values
  pixelScaling = 1.0;
  intensityScale = 1.0;
  imageRotation = 0.0;
}


/* ---------------- PUBLIC METHOD: SetSubsampling ---------------------- */
/// Turn pixel subsampling on or off (true = on, false = off).
void FunctionObject::SetSubsampling( bool subsampleFlag )
{
  doSubsampling = subsampleFlag;
}


/* ---------------- PUBLIC METHOD: SetZeroPoint ------------------------ */
/// Used to specify a magnitude zero point (for *1D* functions).
void FunctionObject::SetZeroPoint( double zeroPoint )
{
  ZP = zeroPoint;
}


/* ---------------- PUBLIC METHOD: SetLabel ---------------------------- */
/// Used to specify a string label for a particular function instance.
void FunctionObject::SetLabel( string &userLabel )
{
  label = userLabel;
}


/* ---------------- PUBLIC METHOD: SetImageParameters ------------------ */
/// Tell FunctionObject instance about updated image-description parameters
/// (for use with multimfit).
/// Note that pixScale = pixel size of image relative to reference image pixels;
/// internally, we will multiply input size parameters by 1/pixScale
void FunctionObject::SetImageParameters( double pixScale, double imageRot, double intensScale )
{
  pixelScaling = 1.0 / pixScale;
  imageRotation = imageRot;
  intensityScale = intensScale;
}


/* ---------------- PUBLIC METHOD: AdjustParametersForImage ------------ */
/// Convert input functionparameters using current set of image-description parameters
void FunctionObject::AdjustParametersForImage( const double inputFunctionParams[], 
										double adjustedFunctionParams[], int offsetIndex )
{
  for (int i = 0; i < nParams; i++)
    adjustedFunctionParams[i + offsetIndex] = inputFunctionParams[i + offsetIndex];
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */
/// Base method for 2D functions: pass current parameters into the function object,
/// storing them for when GetValue() is called, and pre-compute useful quantities.
/// The parameter array params contains *all* parameters for *all* components
/// in the overall model; offsetIndex is used to select the correct starting point
/// for *this* component's parameters.
void FunctionObject::Setup( double params[], int offsetIndex, double xc, double yc )
{
  ;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */
/// Base method for 1D functions: pass current parameters into the function object,
/// storing them for when GetValue() is called, and pre-compute useful quantities.
/// The parameter array params contains *all* parameters for *all* components
/// in the overall model; offsetIndex is used to select the correct starting point
/// for *this* component's parameters.
void FunctionObject::Setup( double params[], int offsetIndex, double xc )
{
  ;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
/// Base method for 2D functions: Compute and return actual function value at
/// pixel coordinates (x,y). This will be called once per pixel.
double FunctionObject::GetValue( double x, double y )
{
  return -1.0;   // dummy (bad) value indicating the user has accidentally
                 // accessed the base class, which should never happen...
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
/// Base method for 1D functions: Compute and return actual function value at
/// specified value of independent variable x.
double FunctionObject::GetValue( double x )
{
  return -1.0;   // dummy (bad) value indicating the user has accidentally
                 // accessed the base class, which should never happen...
}


/* ---------------- PUBLIC METHOD: GetDescription ---------------------- */

string FunctionObject::GetDescription( )
/// Return a string containing function name + short description.
{
  return functionName;
}


/* ---------------- PUBLIC METHOD: GetShortName ----------------------- */

string& FunctionObject::GetShortName( )
/// Return a string containing just the function name.
{
  return shortFunctionName;
}


/* ---------------- PUBLIC METHOD: GetLabel --------------------------- */

string& FunctionObject::GetLabel( )
/// Return a string containing just the function label (will be "" if not set).
{
  return label;
}


/* ---------------- PUBLIC METHOD: GetParameterNames ------------------- */
/// Add this function's parameter names to a vector of strings
void FunctionObject::GetParameterNames( vector<string> &paramNameList )
{
  for (int i = 0; i < nParams; i++)
    paramNameList.push_back(parameterLabels[i]);
}


/* ---------------- PUBLIC METHOD: GetParameterUnits ------------------- */
/// Add this function's parameter unit names (if they exist) to a vector of strings
void FunctionObject::GetParameterUnits( vector<string> &paramUnitList )
{
  for (int i = 0; i < nParams; i++)
    if (HasParameterUnits())
      paramUnitList.push_back(parameterUnits[i]);
    else
      paramUnitList.push_back("");
}


/* ---------------- PUBLIC METHOD: GetExtraParamsDescription ----------- */
/// Return lines describing the user-set extra/optional parameters, if they exist.
void FunctionObject::GetExtraParamsDescription( vector<string> &outputLines )
{
  string  newLine;
  
  outputLines.clear();
  if (extraParamsSet) {
    int  nExtraParams = inputExtraParams.size();
    outputLines.push_back("OPTIONAL_PARAMS_START");
	for (auto const& x : inputExtraParams) {
	  string  key = x.first;
	  string  val = x.second;
	  newLine = PrintToString("%s   %s", key.c_str(), val.c_str());
      outputLines.push_back(newLine);
    }
    outputLines.push_back("OPTIONAL_PARAMS_END");
  }
}


/* ---------------- PUBLIC METHOD: GetNParams -------------------------- */
/// Get number of parameters used by this function.
int FunctionObject::GetNParams( )
{
  return nParams;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

FunctionObject::~FunctionObject()
{
  ;   // nothing to do in base class
}



/* END OF FILE: function_object.cpp ------------------------------------ */
