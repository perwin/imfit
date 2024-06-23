/* FILE: func_sersic.cpp ----------------------------------------------- */
/*
 *   Function object class for a Sersic function, with constant
 * ellipticity and position angle (pure elliptical, not generalized).
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
 *   NOTE: Currently, we assume input PA is in *degrees* [and then we
 * convert it to radians] relative to +x axis.
 *
 *   MODIFICATION HISTORY:
 *     [v0.4]  20--26 Mar 2010: Preliminary support for pixel subsampling.
 *     [v0.3]: 21 Jan 2010: Modified to treat x0,y0 as separate inputs.
 *     [v0.2]: 28 Nov 2009: Updated to new FunctionObject interface.
 *     [v0.1]: 19 Nov 2009: Created (as modification of func_exp.cpp.
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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>

#include "func_sersic.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 5;
const char  PARAM_LABELS[][20] = {"PA", "ell", "n", "I_e", "r_e"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "", "counts/pixel", "pixels"};
const char  FUNCTION_NAME[] = "Sersic function";
const double  DEG2RAD = 0.017453292519943295;
const double PI = 3.14159265358979;
const int  SUBSAMPLE_R = 10;

const char Sersic::className[] = "Sersic";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Sersic::Sersic( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;
  
  // Set up vectors of parameter labels and units
  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }
  parameterUnitsExist = true;
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: AdjustParametersForImage ------------ */
/// Rescale/adjust input function parameters using current set of image-description 
/// parameters
void Sersic::AdjustParametersForImage( const double inputFunctionParams[], 
										double adjustedFunctionParams[], int offsetIndex )
{
  // PA, ell, n, I_e, r_e
  adjustedFunctionParams[0 + offsetIndex] = inputFunctionParams[0 + offsetIndex] - imageRotation;
  adjustedFunctionParams[1 + offsetIndex] = inputFunctionParams[1 + offsetIndex];
  adjustedFunctionParams[2 + offsetIndex] = inputFunctionParams[2 + offsetIndex];
  adjustedFunctionParams[3 + offsetIndex] = intensityScale * inputFunctionParams[3 + offsetIndex];
  adjustedFunctionParams[4 + offsetIndex] = pixelScaling * inputFunctionParams[4 + offsetIndex];
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Sersic::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex] - imageRotation;
  ell = params[1 + offsetIndex];
  n = params[2 + offsetIndex];
  I_e = params[3 + offsetIndex] * intensityScale;
  r_e = params[4 + offsetIndex] * pixelScaling;

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  bn = Calculate_bn(n);
  invn = 1.0 / n;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a Sersic function at radius r,
// with the various parameters and derived values (n, b_n, r_e, etc.)
// pre-calculated by Setup().

double Sersic::CalculateIntensity( double r )
{
  double  intensity;
  
  intensity = I_e * exp( -bn * (pow((r/r_e), invn) - 1.0));
  return intensity;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double Sersic::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  
  nSubsamples = CalculateSubsamples(r);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-pixel
    double deltaSubpix = 1.0 / nSubsamples;
    double x_sub_start = x - 0.5 + 0.5*deltaSubpix;
    double y_sub_start = y - 0.5 + 0.5*deltaSubpix;
    double theSum = 0.0;
    for (int ii = 0; ii < nSubsamples; ii++) {
      double x_ii = x_sub_start + ii*deltaSubpix;
      for (int jj = 0; jj < nSubsamples; jj++) {
        double y_ii = y_sub_start + jj*deltaSubpix;
        x_diff = x_ii - x0;
        y_diff = y_ii - y0;
        xp = x_diff*cosPA + y_diff*sinPA;
        yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
        r = sqrt(xp*xp + yp_scaled*yp_scaled);
        theSum += CalculateIntensity(r);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// Sersic function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int Sersic::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Standard Sersic subsampling, based on Chien Peng's GALFIT algorithm
  if ((doSubsampling) && (r < 10.0)) {
    if ((r_e <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / r_e));
    else {
      if (r <= 4.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}


// Conditional compilation: if GSL library is *not* present, we don't have a
// Gamma function, so we can't calculate the total flux.

#ifdef NO_GSL
/* ---------------- PUBLIC METHOD: CanCalculateTotalFlux --------------- */

bool Sersic::CanCalculateTotalFlux( )
{
  return false;
}


/* ---------------- PUBLIC METHOD: TotalFlux --------------------------- */

double Sersic::TotalFlux( )
{
  return  -1.0;
}

#else
/* ---------------- PUBLIC METHOD: CanCalculateTotalFlux --------------- */

bool Sersic::CanCalculateTotalFlux( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: TotalFlux --------------------------- */

double Sersic::TotalFlux( )
{
  double  bn2n = pow(bn, 2.0*n);
  double  L = 2*PI*n*exp(bn)*I_e*(r_e*r_e);
  return  q*L*gsl_sf_gamma(2*n)/bn2n;
}
#endif

/* END OF FILE: func_sersic.cpp ---------------------------------------- */
