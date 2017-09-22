/* FILE: func_multi-broken-exp.cpp ------------------------------------- */
/*
 *   Function object class for a double broken-exponential function, with constant
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
 *     [v0.1]  31 Aug 2017: Created as modification of func_broken-exp.cpp.
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

#include "func_multi-broken-exp.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 10;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_0", "h1", "h2", "h3",
								"r_break1", "r_break2", "alpha1", "alpha2"};
const char  FUNCTION_NAME[] = "Multi-Broken-Exponential function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char MultiBrokenExponential::className[] = "MultiBrokenExponential";


double CalculateScalingFactor( double h1, double h2, double h3, double r_brk1,
								double r_brk2, double alpha1, double alpha2 );


/* ---------------- CONSTRUCTOR ---------------------------------------- */

MultiBrokenExponential::MultiBrokenExponential( )
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

void MultiBrokenExponential::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_0 = params[2 + offsetIndex ];
  h1 = params[3 + offsetIndex ];
  h2 = params[4 + offsetIndex ];
  h3 = params[5 + offsetIndex ];
  r_b1 = params[6 + offsetIndex ];
  r_b2 = params[7 + offsetIndex ];
  alpha1 = params[8 + offsetIndex ];
  alpha2 = params[9 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  exponent2 = (1.0/alpha1)*(1.0/h1 - 1.0/h2);
  exponent3 = (1.0/alpha2)*(1.0/h2 - 1.0/h3);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  double  S = CalculateScalingFactor(h1, h2, h3, r_b1, r_b2, alpha1, alpha2);
  I_0_times_S = I_0 * S;
  delta_Rb1_scaled = r_b1/h2 - r_b1/h1;
  delta_Rb2_scaled = r_b2/h3 - r_b2/h2;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a broken-exponential function at radius r,
// with the various parameters and derived values (I_0*S, exponent, etc.)
// pre-calculated by Setup().
// NOTE: We assume that r >= 0, since GetValue() ensures that.

double MultiBrokenExponential::CalculateIntensity( double r )
{
  double  P1, P2, P3;
  double  I;
  
  // first piece
  P1 = I_0_times_S * exp(-r/h1);

  // second piece
  // check for possible overflow in exponentiation if r >> r_b1, and re-route around it:
  double scaledR1 = alpha1*(r - r_b1);
  if (scaledR1 > 100.0)
    P2 = exp(delta_Rb1_scaled - r/h2);
  else
	P2 = pow(1.0 + exp(scaledR1), exponent2);

  // third piece
  // check for possible overflow in exponentiation if r >> r_b2, and re-route around it:
  double scaledR2 = alpha2*(r - r_b2);
  if (scaledR2 > 100.0)
    P3 = exp(delta_Rb2_scaled - r/h3);
  else
	P3 = pow(1.0 + exp(scaledR2), exponent3);

  return P1*P2*P3;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double MultiBrokenExponential::GetValue( double x, double y )
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
// broken-exponential function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int MultiBrokenExponential::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // do subsampling of inner exponential only (same as for standard exponential)
  if ((doSubsampling) && (r < 10.0)) {
    if ((h1 <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / h1));
    else {
      if (r <= 3.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}




/* ----------------------------- OTHER FUNCTIONS -------------------------------- */

/* ---------------- CalculateScalingFactor -------------------------------------- */
double CalculateScalingFactor( double h1, double h2, double h3, double r_brk1,
								double r_brk2, double alpha1, double alpha2 )
{
  // n = 2
  double P2a = 1.0 + exp(-alpha1*r_brk1);
  double exp2 = (1.0/alpha1)*(1.0/h1 - 1.0/h2);
  double P2 = pow(P2a, exp2);

  // n = 3
  double P3a = 1.0 + exp(-alpha2*r_brk2);
  double exp3 = (1.0/alpha2)*(1.0/h2 - 1.0/h3);
  double P3 = pow(P3a, exp3);
	
  return 1.0/(P2*P3);
}


/* END OF FILE: func_multi-broken-exp.cpp ------------------------------ */
