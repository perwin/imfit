/* FILE: func_ferrersbar3d.cpp ----------------------------------------- */
/* VERSION 0.1
 *
 *   Experimental function object class for a 3D Ferrers ellipsoid (aka "Ferrers bar"),
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
 *     [v0.1]: 16 Oct 2012: Created (as modification of func_exp3d.cpp.
 */

// Copyright 2012, 2013, 2017 by Peter Erwin.
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
#include <gsl/gsl_errno.h>

#include "func_ferrersbar3d.h"
#include "integrator.h"
#include "helper_funcs_3d.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 8;
const char  PARAM_LABELS[][20] = {"PA", "inc", "barPA", "J_0", "R_bar", "q", "q_z", "n"};
const char  FUNCTION_NAME[] = "FerrersBar3D function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

// Ferrers bar always has zero luminosity density outside its boundaries
const double  INTEGRATION_MULTIPLIER = 1;

const char FerrersBar3D::className[] = "FerrersBar3D";


/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensity_FerrersBar( double s, void *params );




/* ---------------- CONSTRUCTOR ---------------------------------------- */

FerrersBar3D::FerrersBar3D( )
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
  F.function = LuminosityDensity_FerrersBar;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void FerrersBar3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  barPA = params[2 + offsetIndex];
  J_0 = params[3 + offsetIndex ];
  R_bar = params[4 + offsetIndex ];
  q = params[5 + offsetIndex ];
  q_z = params[6 + offsetIndex ];
  n = params[7 + offsetIndex ];

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
  
  a2 = R_bar*R_bar;
  b2 = q*q*a2;
  c2 = q_z*q_z*a2;

  // for Ferrers bar, which has zero luminosity density outside its boundaries,
  // the maximum integLimit can be defined as the larger of R_bar*sin(inc) 
  // [imagine end-on bar inclined w.r.t. line of sight] and c [in the case
  // of a face-on Ferrers bar]
  
  // Possible approach: determine distance from sky plane to major-axis "spine"
  // of bar [fails for perfectly end-on bar!], then work out +/- integration
  // distance around that 
  
  integrationLimit = 1.01 * max((sqrt(b2)*sinInc), R_bar);
//  integrationLimit = 1.01*sqrt(b2);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double FerrersBar3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity, error;
  double  xyParameters[12];
  int  nSubsamples;
  int  nEvals;
  
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
  xyParameters[8] = a2;
  xyParameters[9] = b2;
  xyParameters[10] = c2;
  xyParameters[11] = n;
  F.params = xyParameters;

  // [] NOTE: ideally, we should compute the integration limits directly, given
  // the known orientation of the ellipsoid, as the inner and outer boundaries of
  // the ellipsoid along the current line of sight; this would be superior to our current
  // -number,+number integration...
  // integrate out to +/- integLimit
  totalIntensity = Integrate(F, -integrationLimit, integrationLimit);

  return totalIntensity;
}







/* ----------------------------- OTHER FUNCTIONS -------------------------------- */


/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 * This version is for a Ferrers ellipsoid, which is rotated by barPA degrees relative
 * to component line of nodes.
 */ 
double LuminosityDensity_FerrersBar( double s, void *params )
{
  double  y_d, z_d, x_bar, y_bar, z_bar, m2, lumDensity;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double cosInc = paramsVect[3];
  double sinInc = paramsVect[4];
  double cosBarPA = paramsVect[5];
  double sinBarPA = paramsVect[6];
  double J_0 = paramsVect[7];
  double a2 = paramsVect[8];
  double b2 = paramsVect[9];
  double c2 = paramsVect[10];
  double n = paramsVect[11];
  
  // Determine 3D Cartesian coordinates in bar's native frame of reference
//   Compute3dObjectCoords(s, x_d0, y_d0, z_d0, sinInc, cosInc, cosBarPA, sinBarPA,
// 							x_bar, y_bar, z_bar);
  std::tie(x_bar, y_bar, z_bar) = Compute3dObjectCoords(s, x_d0, y_d0, z_d0, sinInc, cosInc, 
  														cosBarPA, sinBarPA);

  // Calculate luminosity density for Ferrers ellipsoid at x_bar, y_bar, z_bar
  m2 = y_bar*y_bar/a2 + x_bar*x_bar/b2 + z_bar*z_bar/c2;
  if (m2 > 1.0)
  	lumDensity = 0.0;
  else
    lumDensity = J_0 * pow((1.0 - m2), n);
  return lumDensity;
}



/* END OF FILE: func_ferrersbar3d.cpp ---------------------------------- */
