/* FILE: func_pointsource.cpp ------------------------------------------ */
/* 
 *   This is the base class for the various function object classes.
 *   It really shouldn't be instantiated by itself.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x and y.
 *      GetValue() then completes the calculation, using the actual value
 *      of x and y, and returns the result.
 *      So for an image, we expect the user to call Setup() once at
 *      the start, then loop through the pixels of the image, calling
 *      GetValue() to compute the function results for each pixel coordinate
 *      (x,y).
 *
 * Sum_i I_norm,i = 1
 * 
 * I_tot = I_tot Sum_i I_norm,i
 *       = Sum_i (I_tot * I_norm,i)
 *
 *   MODIFICATION HISTORY:
 *     2 Oct 2017: Created (as modification of func_gaussian.cpp).
 */

// Copyright 2017 by Peter Erwin.
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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_pointsource.h"
#include "psf_interpolators.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 1;
const char  PARAM_LABELS[][20] = {"I_tot"};
const char  FUNCTION_NAME[] = "PointSource function";
const double PI = 3.14159265358979;

const char PointSource::className[] = "PointSource";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

PointSource::PointSource( )
{
  string  paramName;
  nParams = N_PARAMS;
  
  functionName = FUNCTION_NAME;
  shortFunctionName = className;   // defined in header file

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = false;
}

/* ---------------- DESTRUCTOR ----------------------------------------- */

PointSource::~PointSource( )
{
  if (interpolatorAllocated)
    delete psfInterpolator;
}


/* ---------------- PUBLIC METHOD: IsPointSource ----------------------- */

bool PointSource::IsPointSource( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: GetInterpolationType ---------------- */

string PointSource::GetInterpolationType( )
{
  return interpolationType;
}


// FIXME: remove this method?
/* ---------------- PUBLIC METHOD: AddPsfData -------------------------- */

void PointSource::AddPsfData( double *psfPixels, int nColumns_psf, int nRows_psf )
{
  psfInterpolator = new PsfInterpolator_bicubic(psfPixels, nColumns_psf, nRows_psf);
  interpolatorAllocated = true;
}


/* ---------------- PUBLIC METHOD: AddPsfInterpolator ------------------ */

void PointSource::AddPsfInterpolator( PsfInterpolator *theInterpolator )
{
  psfInterpolator = theInterpolator;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void PointSource::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  I_tot = params[0 + offsetIndex];
}


/* ---------------- PUBLIC METHOD: HasExtraParams ---------------------- */

bool PointSource::HasExtraParams( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: SetExtraParams ---------------------- */
// Returns -1 if map is empty, 0 if map is not empty but no valid parameter
// name is found. If map has valid parameter name, returns 1 if parameter
// value is OK, -3 if not.
int PointSource::SetExtraParams( map<string,string>& inputMap )
{
  // check for empty map
  if (inputMap.empty()) {
    printf("   PointSource::SetExtraParams: input map is empty!\n");
    return -1;
  }
  // only one possible parameter for this function, so no need to loop
  map<string,string>::iterator iter;
  for( iter = inputMap.begin(); iter != inputMap.end(); iter++) {
    if (iter->first == "method") {
      interpolationType = iter->second;
      // FIXME: do additional checks & setup here
      printf("   PointSource::SetExtraParams -- setting method = %s\n", 
       		interpolationType.c_str());
      extraParamsSet = true;
      return 1;
    }
  }
  return 0;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y).

double PointSource::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  normalizedIntensity;
  
  normalizedIntensity = psfInterpolator->GetValue(x_diff, y_diff);

  return I_tot * normalizedIntensity;
}



/* ---------------- PUBLIC METHOD: CanCalculateTotalFlux --------------- */

bool PointSource::CanCalculateTotalFlux( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: TotalFlux --------------------------- */

double PointSource::TotalFlux( )
{
  return I_tot;
}


/* END OF FILE: func_pointsource.cpp ----------------------------------- */
