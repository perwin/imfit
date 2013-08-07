/* FILE: func_edge-on-disk_n4762.cpp ----------------------------------- */
/* VERSION 0.2
 *
 *
 *   TO-DO:
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
 *     [v0.2]  25 Sept: Corrected calculation of r to be actual cylindrical
 * radius; changed xp and z_perp to R and z; updated subsample calculations
 * to incorporate z-dependence as well as r-dependence.
 *     [v0.1]  24 Sept 2010: Created as modification of func_edge-on-disk.cpp.
 */

// Copyright 2010, 2011, 2012 by Peter Erwin.
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

#include "func_edge-on-disk_n4762.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "I_0", "h2", "r_b", "alpha", "h_z"};
const char  FUNCTION_NAME[] = "Edge-on Disk (NGC 4762 variant) function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char EdgeOnDiskN4762::className[] = "EdgeOnDisk_n4762";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

EdgeOnDiskN4762::EdgeOnDiskN4762( )
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

void EdgeOnDiskN4762::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  I_0 = params[1 + offsetIndex ];
  h2 = params[2 + offsetIndex ];
  r_b = params[3 + offsetIndex ];
  alpha = params[4 + offsetIndex ];
  h_z = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  exponent = (1.0 / (alpha*h2));
  double  S = pow( (1.0 + exp(-alpha*r_b)), exponent );
  I_0_times_S = I_0 * S;
  delta_Rb_scaled = r_b/h2;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */

double EdgeOnDiskN4762::CalculateIntensity( double r, double z )
{
  double  verticalScaling, I_radial;
  double  rr = fabs(r);
  
  // check for possible overflow in exponentiation if r >> r_b, and re-route around it:
  if ( alpha*(rr - r_b) > 100.0 ) {
    // Outer-exponential approximation:
    I_radial = I_0_times_S * exp(delta_Rb_scaled - rr/h2);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    I_radial = I_0_times_S * pow( 1.0 + exp(alpha*(rr - r_b)), -exponent );
  }

  verticalScaling = exp(-fabs(z)/h_z);
  return I_radial * verticalScaling;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double EdgeOnDiskN4762::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  R, z, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame
  R = x_diff*cosPA + y_diff*sinPA;
  z = -x_diff*sinPA + y_diff*cosPA;   // "z" is y in the component reference frame

  nSubsamples = CalculateSubsamples(R, z);
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
        R = x_diff*cosPA + y_diff*sinPA;
        z = -x_diff*sinPA + y_diff*cosPA;
        theSum += CalculateIntensity(R, z);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(R, z);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// r=0 line and/or z away from the disk plane.
// Note that since this particular function object uses a radial profile which is *flat*
// in its inner part, we ignore the r-dependence for subsampling.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int EdgeOnDiskN4762::CalculateSubsamples( double r, double z )
{
  int  nSamples = 1;
  double  z_abs = fabs(z);
  
  // based on standard exponential-function subsampling
  if ( (doSubsampling) && (z_abs < 10.0) ) {
    if ((h_z <= 1.0) && (z_abs <= 1.0)) {
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / h_z));
    }
    else {
      if (z_abs <= 3.0)
        nSamples = 2 * SUBSAMPLE_R;
      else {
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / z_abs));
      }
    }
  }
  return nSamples;
}



/* END OF FILE: func_edge-on-disk_n4762.cpp ---------------------------- */
