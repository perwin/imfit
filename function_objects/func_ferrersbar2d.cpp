/* FILE: func_ferrersbar2d.cpp ----------------------------------------- */
/* 
 *   Function object class for a 2D version of a Ferrers ellipsoid, as used
 * by, e.g., GALFIT. (This is really a 2D cut through the equatorial plane
 * of a Ferrers ellipsoid, with intensity = local density.)
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
 *     [v0.1]  8 May 2019: Created (as modification of func_sersic.cpp).
 */

// Copyright 2019--2022 by Peter Erwin.
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

#include "func_ferrersbar2d.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "ell", "c0", "n", "I_0", "a_bar"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "", "", "counts/pixel", 
								"pixels"};
const char  FUNCTION_NAME[] = "2D version of Ferrers-ellipsoid function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char FerrersBar2D::className[] = "FerrersBar2D";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

FerrersBar2D::FerrersBar2D( )
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

void FerrersBar2D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  c0 = params[2 + offsetIndex ];
  n = params[3 + offsetIndex ];
  I_0 = params[4 + offsetIndex ];
  a_bar = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  // generalized ellipse exponents
  ellExp = c0 + 2.0;
  invEllExp = 1.0 / ellExp;
  // square of bar semi-major axis
  a_bar2 = a_bar*a_bar;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a 2D pseudo-Ferrers ellipsoid
// at radius r (where the intensity is the same as the actual 3D density
// in the plane of a proper 3D Ferrers ellipsoid).

double FerrersBar2D::CalculateIntensity( double r )
{
  double  m2;
  double  lumDensity;  // intensity, but numerically same as Ferrers ellipsoid lum.density
  
  // Calculate luminosity density for Ferrers ellipsoid in the z=0 plane
  m2 = r*r/a_bar2;
  if (m2 > 1.0)
  	lumDensity = 0.0;
  else
    lumDensity = I_0 * pow((1.0 - m2), n);
  return lumDensity;

}


/* ---------------- PRIVATE METHOD: CalculateRadius -------------------- */
// This function calculates the equivalent radius for a generalized ellipse,
// for a coordinate system where r=0 at deltaX = deltaY = 0.

double FerrersBar2D::CalculateRadius( double deltaX, double deltaY )
{
  double  xp, yp_scaled, powerSum;
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio;
  // convert both to absolute values
  xp = fabs(deltaX*cosPA + deltaY*sinPA);
  yp_scaled = fabs((-deltaX*sinPA + deltaY*cosPA)/q);
  // Compute radius for generalized ellipse
  powerSum = pow(xp, ellExp) + pow(yp_scaled, ellExp);
  return pow(powerSum, invEllExp);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double FerrersBar2D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  r, totalIntensity;
  int  nSubsamples;
  
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
// Sersic function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int FerrersBar2D::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Subsampling as for basic (pure-ellipse) Sersic (func_sersic.cpp)
  if ((doSubsampling) && (r < 10.0)) {
    if ((a_bar <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / a_bar));
    else {
      if (r <= 4.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}



/* END OF FILE: func_ferrersbar2d.cpp ---------------------------------- */
