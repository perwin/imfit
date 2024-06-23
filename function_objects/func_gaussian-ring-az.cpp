/* FILE: func_gaussian-ring-az.cpp ------------------------------------- */
/* 
 *   Class derived from FunctionObject; function_object.h) which produces an elliptical 
 * ring with a Gaussian profile.
 *
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
 *     [v0.1]  12 April 2019: Created as modification of func_gaussian-ring.cpp.
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

#include "func_gaussian-ring-az.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "ell", "A_maj", "A_min_rel", "R_ring", "sigma_r"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "counts/pixel", "",
								"pixels", "pixels"};
const char  FUNCTION_NAME[] = "Gaussian Ring with azimuthal variation function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char GaussianRingAz::className[] = "GaussianRingAz";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

GaussianRingAz::GaussianRingAz( )
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


/* ---------------- PUBLIC METHOD: AdjustParametersForImage ------------ */
/// Rescale/adjust input function parameters using current set of image-description 
/// parameters
void GaussianRingAz::AdjustParametersForImage( const double inputFunctionParams[], 
										double adjustedFunctionParams[], int offsetIndex )
{
  // PA, ell, A_maj, A_min_rel, R_ring, sigma_r
  adjustedFunctionParams[0 + offsetIndex] = inputFunctionParams[0 + offsetIndex]  - imageRotation;
  adjustedFunctionParams[1 + offsetIndex] = inputFunctionParams[1 + offsetIndex];
  adjustedFunctionParams[2 + offsetIndex] = intensityScale * inputFunctionParams[2 + offsetIndex];
  adjustedFunctionParams[3 + offsetIndex] = inputFunctionParams[3 + offsetIndex];
  adjustedFunctionParams[4 + offsetIndex] = pixelScaling * inputFunctionParams[4 + offsetIndex];
  adjustedFunctionParams[5 + offsetIndex] = pixelScaling * inputFunctionParams[5 + offsetIndex];
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void GaussianRingAz::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex]  - imageRotation;
  ell = params[1 + offsetIndex];
  A_maj = params[2 + offsetIndex] * intensityScale;
  A_min_rel = params[3 + offsetIndex];
  R_ring = params[4 + offsetIndex] * pixelScaling;   // major-axis radius of ring
  sigma_r = params[5 + offsetIndex] * pixelScaling;

  // pre-compute useful things for this round of invoking the function
  A_min = A_min_rel * A_maj;
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  A_mid = 0.5*(A_maj + A_min);
  deltaA = 0.5*(A_maj - A_min);

  twosigma_squared = 2.0 * sigma_r*sigma_r;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a 1-D Gaussian function with
// center = R_ring at position r.
double GaussianRingAz::CalculateIntensity( double r, double theta )
{
  double  A, I;
  double  r_diff = fabs(r - R_ring);
  
  A = A_mid + deltaA*cos(2*theta);
  I = A * exp(-(r_diff*r_diff)/twosigma_squared);
  return I;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// R is *non-negative*!
double GaussianRingAz::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r, totalIntensity;
  double  theta_image, theta_ellipse;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  theta_image = atan2(y_diff, x_diff);
  theta_ellipse = PA_rad - theta_image;
  
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
        theta_image = atan2(y_diff, x_diff);
        theta_ellipse = PA_rad - theta_image;
        theSum += CalculateIntensity(r, theta_ellipse);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r, theta_ellipse);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a (scaled) distance of r away from the center of the
// ring.
// For now, we don't do any subsampling (because it's easier and also because Gaussians
// usually aren't as prone to subsampling effects).
int GaussianRingAz::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  return nSamples;
}



/* END OF FILE: func_gaussian-ring-az.cpp ------------------------------ */
