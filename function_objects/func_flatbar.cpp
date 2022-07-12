/* FILE: func_flat-bar.cpp --------------------------------------------- */
/*
 *   Function object class for a "flat" (vertically thin) bar with a broken
 * exponential major-axis profile, transitioning to single-exponential for
 * position angles > deltaPA_max.
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
 *     [v0.1]  11 Aug 2018: Created as modification of func_broken-exp.cpp.
 */

// Copyright 2018--2022 by Peter Erwin.
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
#include <tuple>

#include "func_flatbar.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 8;
const char  PARAM_LABELS[][20] = {"PA", "ell", "deltaPA_max", "I_0", "h1", "h2", 
									"r_break", "alpha"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "deg", "counts/pixel", 
								"pixels", "pixels", "pixels", "1/pixels"};
const char  FUNCTION_NAME[] = "FlatBar function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char FlatBar::className[] = "FlatBar";

const double  RAD2DEG = 180.0/M_PI;



/* ---------------- CONSTRUCTOR ---------------------------------------- */

FlatBar::FlatBar( )
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

void FlatBar::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  deltaPA_max = params[2 + offsetIndex];
  I_0 = params[3 + offsetIndex ];
  h1 = params[4 + offsetIndex ];
  h2 = params[5 + offsetIndex ];
  r_b = params[6 + offsetIndex ];
  alpha = params[7 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  deltaPA_max_rad = deltaPA_max * DEG2RAD;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a broken-exponential function at radius r.
// Since h2 and r_b will depend on position angle relative to the bar major
// axis, they have to be re-calculated in GetValue() and passed as parameters
// to this function, instead of simply being pre-calculated in Setup().
// NOTE: We assume that r >= 0, since GetValue() ensures that.

double FlatBar::CalculateIntensity( double r, double h2_adj, double r_b_adj )
{
  double  I;
  
  double  exponent = (1.0/alpha) * (1.0/h1 - 1.0/h2_adj);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  double  S = pow( (1.0 + exp(-alpha*r_b_adj)), (-exponent) );
  double  I_0_times_S = I_0 * S;
  double  delta_Rb_scaled = r_b_adj/h2_adj - r_b_adj/h1;

  // check for possible overflow in exponentiation if r >> r_b, and re-route around it:
  if ( alpha*(r - r_b_adj) > 100.0 ) {
    // Outer-exponential approximation:
    I = I_0_times_S * exp(delta_Rb_scaled - r/h2_adj);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    I = I_0_times_S * exp(-r/h1) * pow( 1.0 + exp(alpha*(r - r_b_adj)), exponent );
  }
  return I;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double FlatBar::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, yp_scaled, r, r_circ, totalIntensity;
  double  r_b_current, h2_current;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp = -x_diff*sinPA + y_diff*cosPA;
  yp_scaled = yp/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  r_circ = sqrt(xp*xp + yp*yp);
  
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
        yp = -x_diff*sinPA + y_diff*cosPA;
        yp_scaled = yp/q;
        r = sqrt(xp*xp + yp_scaled*yp_scaled);
        r_circ = sqrt(xp*xp + yp*yp);
        std::tie(r_b_current, h2_current) = GetAdjustedRbh2(xp, yp, r, r_circ);
        theSum += CalculateIntensity(r, h2_current, r_b_current);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else {
    double  deltaPA = fabs(atan(yp/xp));
    if (deltaPA > deltaPA_max_rad)
      //totalIntensity = I_0 * exp(-r/h1);
      h2_current = h1;
    else {
      std::tie(r_b_current, h2_current) = GetAdjustedRbh2(xp, yp, r, r_circ);
    }
    totalIntensity = CalculateIntensity(r, h2_current, r_b_current);
  }

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// broken-exponential function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int FlatBar::CalculateSubsamples( double r )
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


/* ---------------- PROTECTED METHOD: ParabolicInterpolation ---------------------- */
// Interpolate value of y from y = y1 (at x = x1) to y = y2 (at x = x2), evaluated
// at x.
// Normally, this means interpolating from y = h2 (at deltaPA = 0) to
// y = h1 (at deltaPA = deltaPA_max)
double FlatBar::ParabolicInterpolation( double x, double x1, double x2, double y1, double y2 )
{
  if (x < x1)
    return y1;
  if (x > x2)
    return y2;
  double  x21diff = x2 - x1;
  double  b = (y2 - y1) / (x21diff*x21diff);
  double  xDiff = x - x1;
  return y1 + b*xDiff*xDiff;
}



std::tuple<double, double>  FlatBar::GetAdjustedRbh2( double xp, double yp, 
								double r, double r_circ  )
{
  double  r_b_new, h2_new;
  // compute adjusted break radius
  if ((r_circ > 0.0) && (r > 0.0)) {
    r_b_new = (r/r_circ) * r_b;
    double  deltaPA = fabs(atan(yp/xp));
    h2_new = ParabolicInterpolation(deltaPA, 0.0, deltaPA_max_rad, h2, h1);
  }
  else {
    r_b_new = r_b;
    h2_new = h2;
  }
  return std::make_tuple(r_b_new, h2_new);
}


/* END OF FILE: func_flatbar.cpp --------------------------------------- */
