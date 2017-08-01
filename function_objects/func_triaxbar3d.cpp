/* FILE: func_triaxbar3d.cpp ------------------------------------------- */
/* VERSION 0.1
 *
 *   Experimental function object class for a 3D triaxial ellipsoid,
 * seen at position angle PA and inclination inc.
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
 *     [v0.1]: 1 Aug 2017: Created (as modification of func_ferrersbar3d.cpp.
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
#include <gsl/gsl_errno.h>

#include "func_triaxbar3d.h"
#include "integrator.h"
#include "helper_funcs_3d.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 7;
const char  PARAM_LABELS[][20] = {"PA", "inc", "barPA", "J_0", "sigma", "q", "q_z"};
const char  FUNCTION_NAME[] = "TriaxBar3D function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const double  INTEGRATION_MULTIPLIER = 2;

const char TriaxBar3D::className[] = "TriaxBar3D";


/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensity_TriaxBar( double s, void *params );




/* ---------------- CONSTRUCTOR ---------------------------------------- */

TriaxBar3D::TriaxBar3D( )
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

  // Stuff related to GSL integration  
  gsl_set_error_handler_off();
  F.function = LuminosityDensity_TriaxBar;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void TriaxBar3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  barPA = params[2 + offsetIndex];
  J_0 = params[3 + offsetIndex ];
  sigma = params[4 + offsetIndex ];
  q = params[5 + offsetIndex ];
  q_z = params[6 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  // bar PA rotations are computed relative to +y axis; convert to +x-axis reference
  barPA_rad = (barPA + 90.0) * DEG2RAD;
  cosBarPA = cos(barPA_rad);
  sinBarPA = sin(barPA_rad);
  inc_rad = inclination * DEG2RAD;
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);
  
  twosigma_squared = 2.0 * sigma*sigma;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double TriaxBar3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity, error;
  double  integLimit;
  double  xyParameters[11];
  
  // Calculate x,y in component (projected sky) reference frame, corrected for
  // rotation of line of nodes
  xp = x_diff*cosPA + y_diff*sinPA;
  yp = -x_diff*sinPA + y_diff*cosPA;

  // Calculate (x,y,z)_start point in component's native xyz reference frame, 
  // corresponding to intersection of line-of-sight ray with projected sky frame
  x_d0 = xp;
  y_d0 = yp * cosInc;
  z_d0 = yp * sinInc;

  // Set up parameter vector for the integration (everything that stays unchanged
  // for this particular xp,yp location)
  xyParameters[0] = x_d0;
  xyParameters[1] = y_d0;
  xyParameters[2] = z_d0;
  xyParameters[3] = cosInc;
  xyParameters[4] = sinInc;
  xyParameters[5] = cosBarPA;
  xyParameters[6] = sinBarPA;
  xyParameters[7] = J_0;
  xyParameters[8] = q;
  xyParameters[9] = q_z;
  xyParameters[10] = twosigma_squared;
  F.params = xyParameters;

  // [] NOTE: ideally, we should compute the integration limits directly, given
  // the known orientation of the ellipsoid, as the inner and outer boundaries of
  // the ellipsoid along the current line of sight; this would be superior to our current
  // -number,+number integration...

  // integrate out to +/- integLimit, which is multiple of Gaussian sigma^2
  integLimit = INTEGRATION_MULTIPLIER * twosigma_squared;
  totalIntensity = Integrate(F, -integLimit, integLimit);

  return totalIntensity;
}



// Compute equivalent radius given position (xp, yp, zp) relative to center of
// object for triaxial ellipsoid with b/a = q and c/a = q_z
double CalculateTriaxRadiusSquared_simple( double xp, double yp, double zp, double q, 
											double q_z )
{
  double  xx, yy_scaled, zz_scaled;
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio;
  // convert both to absolute values
  xx = fabs(xp);
  yy_scaled = fabs(yp / q);
  zz_scaled = fabs(zp / q_z);
  // Compute radius for generalized ellipse
  return xx*xx + yy_scaled*yy_scaled + zz_scaled*zz_scaled;
}

// Compute equivalent radius given position (xp, yp, zp) relative to center of
// object for triaxial ellipsoid with generalized ellipsoidal contours
// double CalculateTriaxRadius_general( double xp, double yp, double zp, double q, double q_z )
// {
//   double  xx, yy_scaled, zz_scaled, powerSum;
//   // Calculate x,y in component reference frame, and scale y by 1/axis_ratio;
//   // convert both to absolute values
//   xx = fabs(xp);
//   yy_scaled = fabs(yp / q);
//   zz_scaled = fabs(zp / q_z);
//   // Compute radius for generalized ellipse
//   powerSum = pow(xx, ellExp) + pow(yy_scaled, ellExp) + pow(zz_scaled, ellExp);
//   return pow(powerSum, invEllExp);
// }






/* ----------------------------- OTHER FUNCTIONS -------------------------------- */


/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 * This version is for a triaxial ellipsoid, which is rotated by barPA degrees relative
 * to component line of nodes.
 */ 
double LuminosityDensity_TriaxBar( double s, void *params )
{
  double  x_bar, y_bar, z_bar, r2;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double cosInc = paramsVect[3];
  double sinInc = paramsVect[4];
  double cosBarPA = paramsVect[5];
  double sinBarPA = paramsVect[6];
  double J_0 = paramsVect[7];
  double  q = paramsVect[8];
  double  q_z = paramsVect[9];
  double  twosigma_squared = paramsVect[10];
  
  // Determine 3D Cartesian coordinates in bar's native frame of reference
  Compute3dObjectCoords(s, x_d0, y_d0, z_d0, sinInc, cosInc, cosBarPA, sinBarPA,
							x_bar, y_bar, z_bar);
  
  // Calculate luminosity density for Gaussian ellipsoid at x_bar, y_bar, z_bar
  r2 = CalculateTriaxRadiusSquared_simple(x_bar, y_bar, z_bar, q, q_z);
  return J_0 * exp(-r2/twosigma_squared);
}



/* END OF FILE: func_triaxbar3d.cpp ------------------------------------ */
