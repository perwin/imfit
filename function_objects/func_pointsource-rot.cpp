/* FILE: func_pointsource_rot.cpp -------------------------------------- */
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
 *     11 Aug 2021: Created (as modification of func_pointsource.cpp).
 */

// Copyright 2021--2023 by Peter Erwin.
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

#include "func_pointsource-rot.h"
#include "psf_interpolators.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 2;
const char  PARAM_LABELS[][20] = {"PA", "I_tot"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "counts"};
const char  FUNCTION_NAME[] = "PointSourceRot function";
const double  DEG2RAD = 0.017453292519943295;
const double PI = 3.14159265358979;

const char PointSourceRot::className[] = "PointSourceRot";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

PointSourceRot::PointSourceRot( )
{

  nParams = N_PARAMS;
  
  functionName = FUNCTION_NAME;
  shortFunctionName = className;   // defined in header file

  // Set up vectors of parameter labels and units
  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }
  parameterUnitsExist = true;
  
  oversamplingScale = 1;
  doSubsampling = false;
}

/* ---------------- DESTRUCTOR ----------------------------------------- */

PointSourceRot::~PointSourceRot( )
{
  if (interpolatorAllocated)
    delete psfInterpolator;
}


/* ---------------- PUBLIC METHOD: IsPointSource ----------------------- */

bool PointSourceRot::IsPointSource( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: GetInterpolationType ---------------- */

string PointSourceRot::GetInterpolationType( )
{
  return interpolationType;
}



/* ---------------- PUBLIC METHOD: SetOversamplingScale ---------------- */

void PointSourceRot::SetOversamplingScale( int oversampleScale )
{
  oversamplingScale = oversampleScale;
}


// FIXME: remove this method?
/* ---------------- PUBLIC METHOD: AddPsfData -------------------------- */

void PointSourceRot::AddPsfData( double *psfPixels, int nColumns_psf, int nRows_psf )
{
  psfInterpolator = new PsfInterpolator_bicubic(psfPixels, nColumns_psf, nRows_psf);
  interpolatorAllocated = true;
}


/* ---------------- PUBLIC METHOD: AddPsfInterpolator ------------------ */

void PointSourceRot::AddPsfInterpolator( PsfInterpolator *theInterpolator )
{
  psfInterpolator = theInterpolator;
}


/* ---------------- PUBLIC METHOD: HasExtraParams ---------------------- */

bool PointSourceRot::HasExtraParams( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: SetExtraParams ---------------------- */
// Returns -1 if map is empty, 0 if map is not empty but no valid parameter
// name is found. If map has valid parameter name, returns 1 if parameter
// value is OK, -3 if not.
int PointSourceRot::SetExtraParams( map<string,string>& inputMap )
{
  // check for empty map
  if (inputMap.empty()) {
    printf("   PointSourceRot::SetExtraParams: input map is empty!\n");
    return -1;
  }
  // only one possible parameter for this function, so no need to loop
  map<string,string>::iterator iter;
  for( iter = inputMap.begin(); iter != inputMap.end(); iter++) {
    if (iter->first == "method") {
      if ((iter->second == "bicubic") || (iter->second == "Bicubic")) {
        interpolationType = "bicubic";
        break;
      }
      if ((iter->second == "lanczos2") || (iter->second == "Lanczos2")) {
        interpolationType = "lanczos2";
        break;
      }
      fprintf(stderr, "ERROR: unidentified interpolation type in PointSourceRot::SetExtraParams!\n");
      fprintf(stderr, "(\"%s\" is not a recognized interpolation type)\n",
      			iter->second.c_str());
      return -3;
    }
    else {
      fprintf(stderr, "ERROR: unrecognized extra-parameter name (\"%s\") ", iter->first.c_str());
      fprintf(stderr, " in PointSourceRot::SetExtraParams!\n");
      return 0;
    }
  }
  interpolationType = iter->second;
  extraParamsSet = true;
  printf("   PointSourceRot::SetExtraParams -- setting method = %s\n", 
       		interpolationType.c_str());
  return 1;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void PointSourceRot::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  I_tot = params[1 + offsetIndex];

  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y).
// Note that we multiply x_diff and y_diff by oversamplingScale so that this
// will work correctly when called from an OversampledRegion object.

double PointSourceRot::GetValue( double x, double y )
{
  double  x_diff = oversamplingScale*(x - x0);
  double  y_diff = oversamplingScale*(y - y0);
  double  xp, yp;
  double  normalizedIntensity;

  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp = -x_diff*sinPA + y_diff*cosPA;
  
  normalizedIntensity = psfInterpolator->GetValue(xp, yp);

  return I_tot * normalizedIntensity;
}



/* ---------------- PUBLIC METHOD: CanCalculateTotalFlux --------------- */

bool PointSourceRot::CanCalculateTotalFlux( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: TotalFlux --------------------------- */

double PointSourceRot::TotalFlux( )
{
  return I_tot;
}


/* END OF FILE: func_PointSourceRot_rot.cpp ------------------------------- */
