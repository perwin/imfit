/* FILE: func_gen-exp2.cpp --------------------------------------------- */
/* 
 *   Function object class for a generalized exponential function, with constant
 * (generalized) ellipticity and position angle.
 * This is an experimental version allowing interpolated values of c0.
 *
 *   Generalized ellipse is from Athanassoula et al. (1990), with parameterization
 * as in Peng et al. (2002).
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
 *     [v0.1]  17 July 2018: Created (as modification of func_gen-exp.cpp).
 */

// Copyright 2018 by Peter Erwin.
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

#include "func_gen-exp2.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 8;
const char  PARAM_LABELS[][20] = {"PA", "ell", "c01", "c02", "r1", "r2", "I_0", "h"};
const char FUNCTION_NAME[] = "Generalized-ellipse exponential function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char GenExponential2::className[] = "Exponential_GenEllipse2";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

GenExponential2::GenExponential2( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void GenExponential2::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  c0_1 = params[2 + offsetIndex ];
  c0_2 = params[3 + offsetIndex ];
  r1 = params[4 + offsetIndex ];
  r2 = params[5 + offsetIndex ];
  I_0 = params[6 + offsetIndex ];
  h = params[7 + offsetIndex ];
  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for an exponential function at radius r,
// with the various parameters and derived values (n, b_n, r_e, etc.)
// pre-calculated by Setup().

double GenExponential2::CalculateIntensity( double r )
{
  return I_0 * exp(-r/h);
}


/* ---------------- PRIVATE METHOD: CalculateRadius -------------------- */
// This function calculates the equivalent radius for a generalized ellipse,
// for a coordinate system where r=0 at deltaX = deltaY = 0.

double GenExponential2::CalculateRadius( double deltaX, double deltaY )
{
  double  r_ell, c0, ellExp, invEllExp, powerSum;
  double  xp, yp_scaled;
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio;
  // convert both to absolute values
  xp = fabs(deltaX*cosPA + deltaY*sinPA);
  yp_scaled = fabs((-deltaX*sinPA + deltaY*cosPA)/q);
  r_ell = sqrt(xp*xp + yp_scaled*yp_scaled);
  c0 = LinearInterp(r_ell, r1, r2, c0_1, c0_2);
  ellExp = c0 + 2.0;
  invEllExp = 1.0 / ellExp;

  powerSum = pow(xp, ellExp) + pow(yp_scaled, ellExp);
  return pow(powerSum, invEllExp);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double GenExponential2::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  r_ell, xp, yp_scaled;
  double  r, totalIntensity;
  double  c0_local, ellExp_local, invEllExp_local, powerSum;
  int  nSubsamples;
  
  // Compute radius for generalized ellipse
  // First, compute "radius" corresponding to perfect elliptical isophote
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
//   xp = fabs(x_diff*cosPA + y_diff*sinPA);
//   yp_scaled = fabs((-x_diff*sinPA + y_diff*cosPA)/q);
//   r_ell = sqrt(xp*xp + yp_scaled*yp_scaled);
//   // Now, compute corresponding c0 value, then compute corresponding "radius"
//   // for generalized ellipse
//   c0_local = LinearInterp(r_ell, r1, r2, c0_1, c0_2);
//   // generalized ellipse exponents
//   ellExp_local = c0_local + 2.0;
//   invEllExp_local = 1.0 / ellExp_local;
//   powerSum = pow(xp, ellExp_local) + pow(yp_scaled, ellExp_local);
//   r = pow(powerSum, invEllExp_local);
  
  r = CalculateRadius(x_diff, y_diff);
  
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
        r = CalculateRadius(x_diff, y_diff);
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
// exponential function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int GenExponential2::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Subsampling as for basic (pure-ellipse) exponential (func_exp.cpp)
  if ((doSubsampling) && (r < 10.0)) {
    if ((h <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / h));
    else {
      if (r <= 3.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}



/* END OF FILE: func_gen-exp2.cpp -------------------------------------- */
