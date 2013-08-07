/* FILE: func_sersic.cpp ----------------------------------------------- */
/* VERSION 0.1
 *
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
 *     [v0.1]: 21 Sept 2012: Created (as modification of func_sersic.cpp.
 */

// Copyright 2012, 2013 by Peter Erwin.
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

#include "func_core-sersic.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 8;
const char  PARAM_LABELS[][20] = {"PA", "ell", "n", "I_b", "r_e", "r_b", "alpha", "gamma"};
const char  FUNCTION_NAME[] = "Core-Sersic function";
//const char SHORT_FUNCTION_NAME[] = "Core-Sersic";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char CoreSersic::className[] = "Core-Sersic";

const double  A0_M03 = 0.01945;
const double  A1_M03 = -0.8902;
const double  A2_M03 = 10.95;
const double  A3_M03 = -19.67;
const double  A4_M03 = 13.43;

// The following is the minimum allowable value of r, meant to
// avoid blowups due to the fact that the inner, power-law part of the
// Core-Sersic function becomes infinite at r = 0.
const double  R_MIN = 0.001;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

CoreSersic::CoreSersic( )
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

void CoreSersic::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  n = params[2 + offsetIndex ];
  I_b = params[3 + offsetIndex ];
  r_e = params[4 + offsetIndex ];
  r_b = params[5 + offsetIndex ];
  alpha = params[6 + offsetIndex ];
  gamma = params[7 + offsetIndex ];
// #ifdef DEBUG
//   printf("func_sersic: x0 = %g, y0 = %g, PA = %g, ell = %g, n = %g, r_e = %g, I_e = %g\n",
//           x0, y0, PA, ell, n, r_e, I_e);
// #endif

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  n2 = n*n;
  if (n > 0.36) {
    // Ciotti & Bertin (1999) approximation, good for all n > 0.36
    bn = 2*n - 0.333333333333333 + 0.009876543209876543/n
         + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
         - 7.1510122958919723e-05/(n2*n2);
  } else {
    // MacArthur et al. (2003) approximation for n < 0.36
    bn = A0_M03 + A1_M03*n + A2_M03*n2 + A3_M03*n2*n + A4_M03*n2*n2;
  }
  invn = 1.0 / n;
  Iprime = I_b * pow(2.0, -gamma/alpha) * exp( bn * pow( pow(2.0, 1.0/alpha) * r_b/r_e, (1.0/n) ));
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a Core-Sersic function at radius r,
// with the various parameters and derived values (n, b_n, r_e, etc.)
// pre-calculated by Setup().

double CoreSersic::CalculateIntensity( double r )
{
  double  powerlaw_part, exp_part, intensity;
  
  // kludge to handle cases when r is very close to zero:
  if (r < R_MIN)
    r = R_MIN;
  powerlaw_part = pow( 1.0 + pow(r_b/r, alpha), gamma/alpha );
  exp_part = exp( -bn * pow( ( pow(r, alpha) + pow(r_b, alpha) )/pow(r_e, alpha), 1.0/(alpha*n)) );
  intensity = Iprime * powerlaw_part * exp_part;
  return intensity;
}



/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double CoreSersic::GetValue( double x, double y )
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
    // start in center of leftmost/bottommost sub-picel
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
// Core-Sersic function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int CoreSersic::CalculateSubsamples( double r )
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



/* END OF FILE: func_core-sersic.cpp ----------------------------------- */
