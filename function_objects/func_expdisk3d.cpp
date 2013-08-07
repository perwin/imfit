/* FILE: func_expdisk3d.cpp -------------------------------------------- */
/* VERSION 0.1
 *
 *   Experimental function object class for a 3D exponential disk (luminosity
 * density = radial exponential with scale length h and vertical sech^(2/n) profile
 * with scale heigh h_z), seen at position angle PA and inclination inc.
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
 *     [v0.2]{ 24 Aug: Correct implementation of line-of-sight integration (using
 *   sky plane passing through component center as s=0 start point for integration).
 *     [v0.1]: 18--19 Aug 2012: Created (as modification of func_exp.cpp.
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




// Working Python code:

// def ExpFunc( s, x_d0, y_d0, inc_rad, p ):
// 
// 	I0 = p[0]
// 	h = p[1]
// 	h_z = p[2]
// 	
// 	x_d = x_d0
// 	y_d = y_d0 + s*math.sin(inc_rad)
// 	z_d = -s*math.cos(inc_rad)
// 	
// 	R = math.sqrt(x_d**2 + y_d**2)
// 	z = abs(z_d)
// 	
// 	return I0 * math.exp(-R/h)*math.exp(-z/h_z)
// 
// 
// 
// def DoIntegration( x_d0, y_d0, inc_rad, someFunc, params ):
// 	
// 	I,err = integrate.quad(someFunc, -INT_LIMIT, INT_LIMIT, args=(x_d0,y_d0,inc_rad,params))



/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <gsl/gsl_errno.h>

#include "func_expdisk3d.h"
#include "integrator.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "inc", "J_0", "h", "n", "z_0"};
const char  FUNCTION_NAME[] = "ExponentialDisk3D function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const double  COSH_LIMIT = 100.0;

const double  INTEGRATION_MULTIPLIER = 20;

const char ExponentialDisk3D::className[] = "ExponentialDisk3D";


/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensity( double s, void *params );




/* ---------------- CONSTRUCTOR ---------------------------------------- */

ExponentialDisk3D::ExponentialDisk3D( )
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
  F.function = LuminosityDensity;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void ExponentialDisk3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  J_0 = params[2 + offsetIndex ];
  h = params[3 + offsetIndex ];
  n = params[4 + offsetIndex ];
  z_0 = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  inc_rad = inclination * DEG2RAD;
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);
  //tanInc = tan(inc_rad);
  
  alpha = 2.0/n;
  scaledZ0 = alpha*z_0;
  two_to_alpha = pow(2.0, alpha);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double ExponentialDisk3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity, error;
  double  integLimit;
  double  xyParameters[11];
  int  nSubsamples;
  int  nEvals;
  
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
  xyParameters[6] = h;
  xyParameters[7] = z_0;
  xyParameters[8] = scaledZ0;
  xyParameters[9] = two_to_alpha;
  xyParameters[10] = alpha;
  F.params = xyParameters;

  // integrate out to +/- integLimit, which is multiple of exp. scale length
  // (NOTE: it seems like it would be faster to precalculate integLimit in the
  // Setup() call above; for some reason doing it that way makes the whole thing
  // take ~ 4 times longer!)
  integLimit = INTEGRATION_MULTIPLIER * h;
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
double LuminosityDensity( double s, void *params )
{
  double  y_d, z_d, z, R, lumDensity;
  double  verticalScaling, sech;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double cosInc = paramsVect[3];
  double sinInc = paramsVect[4];
  double J_0 = paramsVect[5];
  double h = paramsVect[6];
  double z_0 = paramsVect[7];
  double scaledZ0 = paramsVect[8];
  double two_to_alpha = paramsVect[9];
  double alpha = paramsVect[10];
  
  // Given s and the pre-defined parameters, determine our 3D location (x_d,y_d,z_d)
  // [by construction, x_d = x_d0]
  y_d = y_d0 + s*sinInc;
  z_d = z_d0 - s*cosInc;
  
  // Convert 3D Cartesian coordinate to R,z coordinate
  R = sqrt(x_d0*x_d0 + y_d*y_d);
  z = fabs(z_d);

  // if combination of n*z/z_0 is large enough, switch to simple exponential
  if ((z/scaledZ0) > COSH_LIMIT)
    verticalScaling = two_to_alpha * exp(-z/z_0);
  else {
    sech = 1.0 / cosh(z/scaledZ0);
    verticalScaling = pow(sech, alpha);
  }

  lumDensity = J_0 * exp(-R/h) * verticalScaling;
  return lumDensity;
}



/* END OF FILE: func_expdisk3d.cpp ------------------------------------- */
