/* FILE: func_broken-exp2d.cpp ----------------------------------------- */
/* VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces an approximate 2D edge-on broken-exponential (with radial
 * broken-exponential profile and vertical exponential profile).
 *
 *   TO-DO:
 *      [] Correct subsampling calcs to use z properly.
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
 *     [v0.1]  16 June 2010: Created as modification of func_broken-exp.cpp.
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


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_broken-exp2d.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 7;
const char  PARAM_LABELS[][20] = {"PA", "I_0", "h1", "h2", "r_break", "alpha", "h_z"};
const char  FUNCTION_NAME[] = "Broken-Exponential2D function";
//const char  SHORT_FUNCTION_NAME[] = "BrokenExponential2D";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char BrokenExponential2D::className[] = "BrokenExponential2D";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

BrokenExponential2D::BrokenExponential2D( )
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

void BrokenExponential2D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  I_0 = params[1 + offsetIndex ];
  h1 = params[2 + offsetIndex ];
  h2 = params[3 + offsetIndex ];
  r_b = params[4 + offsetIndex ];
  alpha = params[5 + offsetIndex ];
  h_z = params[6 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  exponent = (1.0/alpha) * (1.0/h1 - 1.0/h2);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  double  S = pow( (1.0 + exp(-alpha*r_b)), (-exponent) );
  I_0_times_S = I_0 * S;
  delta_Rb_scaled = r_b/h2 - r_b/h1;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a broken-exponential function 
// with vertical scale height h_z, at cylindrical radius r (along the major-axis of the
// component) and height z (perpendicular to the major-axis of the component).

double BrokenExponential2D::CalculateIntensity( double r, double z )
{
  double  I;
  double  rr = fabs(r);
  
  // check for possible overflow in exponentiation if rr >> r_b, and re-route around it:
  if ( alpha*(rr - r_b) > 100.0 ) {
    // Outer-exponential approximation:
    I = I_0_times_S * exp(delta_Rb_scaled - rr/h2);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    I = I_0_times_S * exp(-rr/h1) * pow( 1.0 + exp(alpha*(rr - r_b)), exponent );
  }
  return I * exp(-fabs(z)/h_z);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double BrokenExponential2D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, z_perp, r, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  z_perp = -x_diff*sinPA + y_diff*cosPA;   // "z" is y in the component reference frame
  r = sqrt(xp*xp + z_perp*z_perp);

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
        z_perp = -x_diff*sinPA + y_diff*cosPA;
        r = sqrt(xp*xp + z_perp*z_perp);
        theSum += CalculateIntensity(xp, z_perp);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(xp, z_perp);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// broken-exponential function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int BrokenExponential2D::CalculateSubsamples( double r )
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



/* END OF FILE: func_broken-exp2d.cpp ---------------------------------- */
