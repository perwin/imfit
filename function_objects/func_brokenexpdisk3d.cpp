/* FILE: func_brokenexp3d.cpp ------------------------------------------ */
/* VERSION 0.2
 *
 *   Experimental function object class for a 3D exponential disk (luminosity
 * density = radial broken exponential with scale lengths h1 and h2, break
 * radius r_b, and vertical vertical sech^(2/n) profile with scale heigh h_z, seen
 * at position angle PA and inclination inc.
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
 *     [v0.2]: 10 Aug 2013: Modified to use sech^(2/n) instead of simple exponential
 * as vertical profile.
 *     [v0.1]: 25 Oct 2012: Created (as modification of func_expdisk3d.cpp).
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
#include <gsl/gsl_errno.h>

#include "func_brokenexpdisk3d.h"
#include "integrator.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 9;
const char  PARAM_LABELS[][20] = {"PA", "inc", "J_0", "h1", "h2", "r_break", "alpha", "n", "z_0"};
const char  FUNCTION_NAME[] = "BrokenExponentialDisk3D function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const double  COSH_LIMIT = 100.0;

const double  INTEGRATION_MULTIPLIER = 20;

const char BrokenExponentialDisk3D::className[] = "BrokenExponentialDisk3D";


/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensityBED( double s, void *params );




/* ---------------- CONSTRUCTOR ---------------------------------------- */

BrokenExponentialDisk3D::BrokenExponentialDisk3D( )
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
  F.function = LuminosityDensityBED;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void BrokenExponentialDisk3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  J_0 = params[2 + offsetIndex ];
  h1 = params[3 + offsetIndex ];
  h2 = params[4 + offsetIndex ];
  r_b = params[5 + offsetIndex ];
  alpha = params[6 + offsetIndex ];
//  h_z = params[7 + offsetIndex ];
  n = params[7 + offsetIndex ];
  z_0 = params[8 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  inc_rad = inclination * DEG2RAD;
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);
  //tanInc = tan(inc_rad);
  
  // broken-exponential stuff
  exponent = (1.0/alpha) * (1.0/h1 - 1.0/h2);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  double  S = pow( (1.0 + exp(-alpha*r_b)), (-exponent) );
  J_0_times_S = J_0 * S;
  delta_Rb_scaled = r_b/h2 - r_b/h1;
  
  // vertical-profile stuff
  alphaVert = 2.0/n;
  scaledZ0 = alphaVert*z_0;
  two_to_alpha = pow(2.0, alphaVert);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double BrokenExponentialDisk3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity;
  double  integLimit;
  double  xyParameters[17];
  
  // Calculate x,y in component (projected sky) reference frame
  xp = x_diff*cosPA + y_diff*sinPA;
  yp = -x_diff*sinPA + y_diff*cosPA;

  // Calculate (x,y,z)_start in component's native xyz reference frame, corresponding to
  // intersection of line-of-sight ray with projected sky frame
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
  xyParameters[5] = J_0;
  xyParameters[6] = h1;
  xyParameters[7] = h2;
  xyParameters[8] = r_b;
  xyParameters[9] = alpha;
//  xyParameters[10] = h_z;
  xyParameters[10] = J_0_times_S;
  xyParameters[11] = delta_Rb_scaled;
  xyParameters[12] = exponent;
  xyParameters[13] = z_0;
  xyParameters[14] = scaledZ0;
  xyParameters[15] = two_to_alpha;
  xyParameters[16] = alphaVert;
  F.params = xyParameters;

  // integrate out to +/- integLimit, which is larger of (multiple of break radius)
  // and (multiple of h2)
  // (NOTE: it seems like it would be faster to precalculate integLimit in the
  // Setup() call above; for some reason doing it that way makes the whole thing
  // take ~ 4 times longer!)
  integLimit = fmax(INTEGRATION_MULTIPLIER * r_b, INTEGRATION_MULTIPLIER * h2);
//  printf("x,y = %f,%f; x_d0,y_d0,z_d0 = %f,%f,%f; cosInc,sinInc = %f,%f\n",
//  	x,y, x_d0,y_d0,z_d0, cosInc,sinInc);
  totalIntensity = Integrate(F, -integLimit, integLimit);

  return totalIntensity;
}







/* ----------------------------- OTHER FUNCTIONS -------------------------------- */


/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 */ 
double LuminosityDensityBED( double s, void *params )
{
  double  y_d, z_d, z, R, J_rad, lumDensity;
  double  verticalScaling, sech;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double cosInc = paramsVect[3];
  double sinInc = paramsVect[4];
  double J_0 = paramsVect[5];
  double h1 = paramsVect[6];
  double h2 = paramsVect[7];
  double r_b = paramsVect[8];
  double alpha = paramsVect[9];
//  double h_z = paramsVect[10];
  double J_0_times_S = paramsVect[10];
  double delta_Rb_scaled = paramsVect[11];
  double exponent = paramsVect[12];
  double z_0 = paramsVect[13];
  double scaledZ0 = paramsVect[14];
  double two_to_alpha = paramsVect[15];
  double alphaVert = paramsVect[16];
  
  // Given s and the pre-defined parameters, determine our 3D location (x_d,y_d,z_d)
  // [by construction, x_d = x_d0]
  y_d = y_d0 + s*sinInc;
  z_d = z_d0 - s*cosInc;
  
  // Convert 3D Cartesian coordinate to R,z coordinate
  R = sqrt(x_d0*x_d0 + y_d*y_d);
  z = fabs(z_d);

  // Calculate radial component
  // check for possible overflow in exponentiation if r >> r_b, and re-route around it:
  if ( alpha*(R - r_b) > 100.0 ) {
    // Outer-exponential approximation:
    J_rad = J_0_times_S * exp(delta_Rb_scaled - R/h2);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    J_rad = J_0_times_S * exp(-R/h1) * pow( 1.0 + exp(alpha*(R - r_b)), exponent );
  }

  // if combination of n*z/z_0 is large enough, switch to simple exponential
  if ((z/scaledZ0) > COSH_LIMIT)
    verticalScaling = two_to_alpha * exp(-z/z_0);
  else {
    sech = 1.0 / cosh(z/scaledZ0);
    verticalScaling = pow(sech, alphaVert);
  }

  lumDensity = J_rad * verticalScaling;
  return lumDensity;
}



/* END OF FILE: func_brokenexp3d.cpp ----------------------------------- */
