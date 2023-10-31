/* FILE: func_bpbar3d.cpp ---------------------------------------------- */
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
 *     [v0.1]  5 Aug 2023: Created.
 */

// Copyright 2023 by Peter Erwin.
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

#include "func_bpbar3d.h"
#include "integrator.h"
#include "helper_funcs_3d.h"

using namespace std;


// Global (common) parameters
// PA
// inc
// barPA
// z_0 = vertical scale height for Part 1 and base vertical scale height for Part 2
// 
// 1. Flat outer part
// J_0
// ell_outer
// h1
// h2
// R_brk
// alpha
// 
// 2. BP bulge
// J_0
// ell_bp
// h_bp
// 
// z_bp_max [same as A_pea in Datthatri et al.]
// 	To ensure z_0 has same meaning for both components, the "peanut scaling"
// 	should be *mulitiplicative*
// 		z_h(x,y,z) = Z_max * [Gaussian(x,y,R_max,sigma_bp) + Gaussian(x,y,-R_max,sigma_bp)] * z_0
// R_max [same as R_pea in Datthatri et al.]
// sigma_bp [same as sigma_pea in Datthatri et al.]



/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 16;
const char  PARAM_LABELS[][20] = {"PA", "inclination", "barPA", "z_0", "ell_outer", 
									"J_0_outer", "h1", "h2", "r_break", "alpha", "ell_bp", 
									"J_0_bp", "h_bp", "r_bp_max", "z_bp_max", "sigma_bp"};
const char  FUNCTION_NAME[] = "BPBar3D function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char BPBar3D::className[] = "BPBar3D";

const double  INTEGRATION_MULTIPLIER = 5;
const double  RAD2DEG = 180.0/M_PI;


/* ---------------- Local Functions ------------------------------------ */

double CalculateIntensity_BrokenExp( double r, double h1, double h2_adj, double r_b_adj, double alpha );
double DoubleGaussian2D( double x, double y, double R_max, double twosigma_squared );
double LuminosityDensity_BPBar3D( double s, void *params );



/* ---------------- CONSTRUCTOR ---------------------------------------- */

BPBar3D::BPBar3D( )
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
  F.function = LuminosityDensity_BPBar3D;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void BPBar3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  double  ell_outer, ell_bp;
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  barPA = params[2 + offsetIndex];
  z_0 = params[3 + offsetIndex];
  ell_outer = params[4 + offsetIndex];
  J_0_outer = params[5 + offsetIndex ];
  h1 = params[6 + offsetIndex ];
  h2 = params[7 + offsetIndex ];
  r_b = params[8 + offsetIndex ];
  alpha = params[9 + offsetIndex ];
  ell_bp = params[10 + offsetIndex ];
  J_0_bp = params[11 + offsetIndex ];
  r_bp_max = params[12 + offsetIndex ];
  z_bp_max = params[13 + offsetIndex ];
  h_bp = params[14 + offsetIndex ];
  sigma_bp = params[15 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  double  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  // bar PA rotations are computed relative to +y axis; convert to +x-axis reference
  double  barPA_rad = (barPA + 90.0) * DEG2RAD;
  cosBarPA = cos(barPA_rad);
  sinBarPA = sin(barPA_rad);
  double  inc_rad = inclination * DEG2RAD;
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);
  q_outer = 1.0 - ell_outer;
  q_bp = 1.0 - ell_bp;
  twosigma_squared = 2.0 * sigma_bp*sigma_bp;

  integrationLimit = INTEGRATION_MULTIPLIER * r_b;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double BPBar3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity, error;
  double  xyParameters[20];
  int  nSubsamples;
  int  nEvals;
  
  // printf("BPBar3D::GetValue -- x,y = %f,%f\n", x, y);
  // Calculate xp,yp in component (projected sky) reference frame, corrected for
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
  xyParameters[7] = z_0;
  xyParameters[8] = q_outer;
  xyParameters[9] = J_0_outer;
  xyParameters[10] = h1;
  xyParameters[11] = h2;
  xyParameters[12] = r_b;
  xyParameters[13] = alpha;
  xyParameters[14] = q_bp;
  xyParameters[15] = J_0_bp;
  xyParameters[16] = r_bp_max;
  xyParameters[17] = z_bp_max;
  xyParameters[18] = h_bp;
  xyParameters[19] = twosigma_squared;
  F.params = xyParameters;

  // printf("Calling Integrate...\n");
  // integrate out to +/- integLimit
  totalIntensity = Integrate(F, -integrationLimit, integrationLimit);
  // printf("Done.\n");

  return totalIntensity;
}



/* ----------------------------- OTHER FUNCTIONS -------------------------------- */

double CalculateIntensity_BrokenExp( double r, double h1, double h2_adj, double r_b_adj, double alpha )
{
  double  I;
  
  double  exponent = (1.0/alpha) * (1.0/h1 - 1.0/h2_adj);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  double  S = pow( (1.0 + exp(-alpha*r_b_adj)), (-exponent) );
  double  delta_Rb_scaled = r_b_adj/h2_adj - r_b_adj/h1;

  // check for possible overflow in exponentiation if r >> r_b, and re-route around it:
  if ( alpha*(r - r_b_adj) > 100.0 ) {
    // Outer-exponential approximation:
    I = S * exp(delta_Rb_scaled - r/h2_adj);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    I = S * exp(-r/h1) * pow( 1.0 + exp(alpha*(r - r_b_adj)), exponent );
  }
  return I;
}


double DoubleGaussian2D( double x, double y, double R_max, double twosigma_squared )
{
  double  x_diff_right2 = (x - R_max)*(x - R_max);
  double  x_diff_left2 = (x + R_max)*(x + R_max);
  double  y_diff2 = y*y;
  // right_part = contribution from Gaussian centered in right half of bar (x > 0)
  double  right_part = exp( -x_diff_right2/twosigma_squared - y_diff2/twosigma_squared );
  // left_part = contribution from Gaussian centered in left half of bar (x < 0)
  double  left_part = exp( -x_diff_left2/twosigma_squared - y_diff2/twosigma_squared );
  return right_part + left_part;
}


/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 * This version is for a "BP bar", which is rotated by barPA degrees relative
 * to component line of nodes.
 */ 
double LuminosityDensity_BPBar3D( double s, void *params )
{
  double  y_d, z_d, x_bar, y_bar, z_bar;
  double  y_bar_scaled, R_outer;
  double  y_bp_scaled, R_bp, z_h_bp, bp_z_scaling;
  double  lumDensity_outer, lumDensity_bp;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double cosInc = paramsVect[3];
  double sinInc = paramsVect[4];
  double cosBarPA = paramsVect[5];
  double sinBarPA = paramsVect[6];
  double z_0 = paramsVect[7];
  double q_outer = paramsVect[8];
  double J_0_outer = paramsVect[9];
  double h1 = paramsVect[10];
  double h2 = paramsVect[11];
  double r_b = paramsVect[12];
  double alpha = paramsVect[13];
  double q_bp = paramsVect[14];
  double J_0_bp = paramsVect[15];
  double r_bp_max = paramsVect[16];
  double z_bp_max = paramsVect[17];
  double h_bp = paramsVect[18];
  double twosigma_squared = paramsVect[19];

  // Determine 3D Cartesian coordinates in bar's native frame of reference
  std::tie(x_bar, y_bar, z_bar) = Compute3dObjectCoords(s, x_d0, y_d0, z_d0, sinInc, 
  														cosInc, cosBarPA, sinBarPA);

  // * Outer (~FlatBar) component
  y_bar_scaled = y_bar/q_outer;
  R_outer = sqrt(x_bar*x_bar + y_bar_scaled*y_bar_scaled);
  lumDensity_outer = CalculateIntensity_BrokenExp(R_outer, h1, h2, r_b, alpha);
  lumDensity_outer = J_0_outer * lumDensity_outer * exp(-z_bar/z_0);

  // * BP component
  y_bp_scaled = y_bar/q_bp;
  // computer modified vertical scale height
  bp_z_scaling = z_bp_max * DoubleGaussian2D(x_bar, y_bp_scaled, r_bp_max, twosigma_squared);
  z_h_bp = bp_z_scaling * z_0;
  R_bp = sqrt(x_bar*x_bar + y_bp_scaled*y_bp_scaled);
//   lumDensity_bp = J_0_bp * exp(-R_bp/h_bp) * exp(-z_bar/z_0);
  lumDensity_bp = (J_0_bp / bp_z_scaling) * exp(-R_bp/h_bp) * exp(-z_bar/z_h_bp);
  
  return lumDensity_outer + lumDensity_bp;
}


/* END OF FILE: func_bpbar3d.cpp --------------------------------------- */

