/* FILE: func_king2.cpp ------------------------------------------------ */
/* 
 *   Function object class for an King-model function, with constant
 * ellipticity and position angle (pure elliptical, not generalized).
 *   This is the variant where we have the concentration as a free
 * parameter, instead of the tidal/truncation radius r_t. For simplicity,
 * the code is identical to that in func_king.cpp, except that we
 * use c to derive r_t.
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
 *     [v0.1]: 20 Sept 2015: Created as modification of func_king.cpp
 */

// Copyright 2015--2022 by Peter Erwin.
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

#include "func_king2.h"

using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_0", "r_c", "c", "alpha"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "counts/pixel", 
								"pixels", "", ""};
const char  FUNCTION_NAME[] = "Modified King 2 function";
const double  DEG2RAD = 0.017453292519943295;
const double PI  =3.14159265358979;
const int  SUBSAMPLE_R = 10;

const char ModifiedKing2::className[] = "ModifiedKing2";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModifiedKing2::ModifiedKing2( )
{
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up vectors of parameter labels and units
  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void ModifiedKing2::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_0 = params[2 + offsetIndex ];
  r_c = params[3 + offsetIndex ];
  c = params[4 + offsetIndex ];
  alpha = params[5 + offsetIndex ];   // alpha = 2 for standard King model
  
  // extract r_t from c and r_c
  r_t = c * r_c;
  
  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  
  one_over_alpha = 1.0 / alpha;
  one_over_rc = 1.0 / r_c;
  constantTerm = 1.0 / pow(1.0 + (r_t/r_c)*(r_t/r_c), one_over_alpha);
  I_1 = I_0 * pow(1.0 - constantTerm, -alpha);
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
/// This function calculates the intensity for a modified King model function at radius r,
/// with the various parameters and derived values pre-calculated by Setup().
double ModifiedKing2::CalculateIntensity( double r )
{
  double  intensity = 0.0;
  if (r < r_t) {  // ensure we return 0 when r >= r_t
    double  r_over_rc = r * one_over_rc;
    double  variableTerm = 1.0 / pow(1.0 + r_over_rc*r_over_rc, one_over_alpha);
    double  secondPart = pow(variableTerm - constantTerm, alpha);
    intensity = I_1 * secondPart;
  }
  return intensity;
}



/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double ModifiedKing2::GetValue( double x, double y )
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
/// Function which determines the number of pixel subdivisions for sub-pixel integration,
/// given that the current pixel is a distance of r away from the center of the
/// modified King function.
///
/// SLIGHTLY KLUDGEY -- this is just the interpolation scheme for the Exponential function,
/// with r_c replacing h.
///
/// This function returns the number of x and y subdivisions; the total number of subpixels
/// will then be the return value *squared*.
int ModifiedKing2::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Standard exponential subsampling, based on Chien Peng's GALFIT algorithm
  if ((doSubsampling) && (r < 10.0)) {
    if ((r_c <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / r_c));
    else {
      if (r <= 3.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}


/* ---------------- PUBLIC METHOD: ParameterUnitsExist ----------------- */

bool ModifiedKing2::ParameterUnitsExist( )
{
  return true;
}


/* END OF FILE: func_king2.cpp ----------------------------------------- */
