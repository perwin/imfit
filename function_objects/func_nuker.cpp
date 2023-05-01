/* FILE: func_nuker.cpp ------------------------------------------------ */
/* 
 *
 *   Function object class for a Nuker-law function, with constant
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
 *     [v0.1]: 20 April 2023: Created (as modification of func_sersic.cpp.
 */

// Copyright 2023 by Peter Erwin.
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

#include "func_nuker.h"
// #include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 7;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_b", "r_b", "alpha", "beta", "gamma"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "counts/pixel",
								"pixels", "", "", ""};
const char  FUNCTION_NAME[] = "NukerLaw function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char NukerLaw::className[] = "NukerLaw";

// The following is the minimum allowable value of r, meant to
// avoid blowups due to the fact that the inner, power-law part of the
// Nuker-Law function becomes infinite at r = 0.
const double  R_MIN = 0.001;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

NukerLaw::NukerLaw( )
{
  
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


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void NukerLaw::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_b = params[2 + offsetIndex ];
  r_b = params[3 + offsetIndex ];
  alpha = params[4 + offsetIndex ];
  beta = params[5 + offsetIndex ];
  gamma = params[6 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  exponent = (beta - gamma)/alpha;
  Iprime = I_b * pow(2.0, exponent);
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a Nuker-law function at radius r,
// with the various parameters and derived values pre-calculated by Setup().

double NukerLaw::CalculateIntensity( double r )
{
  double  I1, I2;
  
  // kludge to handle cases when r is very close to zero:
  if (r < R_MIN)
    r = R_MIN;

  I1 = Iprime * pow(r_b/r, gamma);
  I2 = pow((1.0 + pow(r/r_b, alpha)), -exponent);
  return (I1*I2);
}



/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double NukerLaw::GetValue( double x, double y )
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
// Nuker-law function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int NukerLaw::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  if ((doSubsampling) && (r < 10.0)) {
    if (r <= 4.0)
      nSamples = 2 * SUBSAMPLE_R;
    else
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
  }
  return nSamples;
}



/* END OF FILE: func_nuker.cpp ----------------------------------------- */
