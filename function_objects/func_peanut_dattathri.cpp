/* FILE: func_peanut_dattathri.cpp ------------------------------------- */
/* 
 *   Function object class for a 3D peanut-shaped/boxy bar model, as described
 * in Dattathri et al. (2024; MNRAS 530: 1195).
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
 *   NOTE: We assume input PA is in *degrees* [and then we convert it to radians] 
 * relative to +x axis.
 *
 */	

// Copyright 2023 by Shashank Datthatri.
// Copyright 2012--2016,2024 by Peter Erwin.
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

#include "func_peanut_dattathri.h"
#include "integrator.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 12;
const char  PARAM_LABELS[][20] = {"PA_lon", "inc", "barPA", "J0_bar", "R_bar_x", "q", "q_z", 
				"R_peanut", "A_peanut", "sigma_peanut", "c_bar_par", "c_bar_perp"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "deg", "deg", "counts/voxel", 
				"pixels", "", "", "pixels", "pixels", "pixels", "", ""};
const char  FUNCTION_NAME[] = "Peanut function (v2)";
const double  DEG2RAD = 0.017453292519943295;

const double  INTEGRATION_MULTIPLIER = 10;

const char DattathriPeanut3D::className[] = "DattathriPeanut3D";


/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensityDattathriPeanut3D( double s, void *params );




/* ---------------- CONSTRUCTOR ---------------------------------------- */

DattathriPeanut3D::DattathriPeanut3D( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up vectors of parameter labels and units
  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }
  parameterUnitsExist = true;

  // Stuff related to GSL integration  
  gsl_set_error_handler_off();
  F.function = LuminosityDensityDattathriPeanut3D;
  
  doSubsampling = false;
}



/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void DattathriPeanut3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inc = params[1 + offsetIndex];
  psi_bar = params[2 + offsetIndex];
  J0_bar = params[3 + offsetIndex];
  R_bar_x = params[4 + offsetIndex];
  q = params[5 + offsetIndex];
  q_z = params[6 + offsetIndex];
  R_peanut = params[7 + offsetIndex];
  A_peanut = params[8 + offsetIndex];
  sigma_peanut = params[9 + offsetIndex];
  c_bar_par = params[10 + offsetIndex];
  c_bar_perp = params[11 + offsetIndex];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  inc_rad = inc * DEG2RAD;
  double  barPA_rad = psi_bar * DEG2RAD;

  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);
  cosPsi = cos(barPA_rad);
  sinPsi = sin(barPA_rad);

  R_bar_y = q * R_bar_x;
  R_bar_z = q_z * R_bar_x;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double DattathriPeanut3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity;
  double  integLimit;
  double  xyParameters[16];
//   int  nSubsamples;
  
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
  xyParameters[3] = J0_bar;
  xyParameters[4] = R_bar_x;
  xyParameters[5] = R_bar_y;
  xyParameters[6] = R_bar_z;
  xyParameters[7] = R_peanut;
  xyParameters[8] = A_peanut;
  xyParameters[9] = sigma_peanut;
  xyParameters[10] = c_bar_par;
  xyParameters[11] = c_bar_perp;
  xyParameters[12] = cosInc;
  xyParameters[13] = sinInc;
  xyParameters[14] = cosPsi;
  xyParameters[15] = sinPsi;
  F.params = xyParameters;

  // integrate out to +/- integLimit, which is multiple of exp. scale length
  // (NOTE: it seems like it would be faster to precalculate integLimit in the
  // Setup() call above; for some reason doing it that way makes the whole thing
  // take ~ 4 times longer!)
  integLimit = INTEGRATION_MULTIPLIER * R_bar_x;
  totalIntensity = Integrate(F, -integLimit, integLimit);

  return totalIntensity;
}





/* ----------------------------- OTHER FUNCTIONS -------------------------------- */

/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 */ 

double LuminosityDensityDattathriPeanut3D( double s, void *params )
{
  double  y_d, z_d, z, R, lumDensity;
  double  verticalScaling, sech;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double J0_bar = paramsVect[3];
  double R_bar_x = paramsVect[4];
  double R_bar_y = paramsVect[5];
  double R_bar_z = paramsVect[6];
  double R_peanut = paramsVect[7];
  double A_peanut = paramsVect[8];
  double sigma_peanut = paramsVect[9];
  double c_bar_par = paramsVect[10];
  double c_bar_perp = paramsVect[11];
  double cosInc = paramsVect[12];
  double sinInc = paramsVect[13];
  double cosPsi = paramsVect[14];
  double sinPsi = paramsVect[15];
  
  // Given s and the pre-defined parameters, determine our 3D location (x_d,y_d,z_d)
  // [by construction, x_d = x_d0]
  y_d = y_d0 + s*sinInc;
  z_d = z_d0 - s*cosInc;
  
  double xbar = x_d0*cosPsi + y_d*sinPsi;
  double ybar = -x_d0*sinPsi + y_d*cosPsi;
  double zbar = z_d;

  // vertical scale height at this (xbar, ybar) location
  double z0 = A_peanut*exp(-( pow((xbar - R_peanut),2.0)/(2.0*pow(sigma_peanut,2.0)) + 
  				pow(ybar,2.0)/(2.*pow(sigma_peanut,2.0)) )) + 
  				A_peanut*exp(-( pow((xbar + R_peanut),2.0)/(2.0* pow(sigma_peanut,2.0)) + 
  				pow(ybar,2.0)/(2.0*pow(sigma_peanut,2.0)) )) + R_bar_z;

  // scaled 3D radius value, for use in computing density
  double Rs = pow( (pow( fabs(xbar)/R_bar_x, c_bar_perp ) + 
  				pow( fabs(ybar)/R_bar_y, c_bar_perp ) ), c_bar_par/c_bar_perp ) + 
  				pow( fabs(zbar)/z0, c_bar_par);
  Rs = pow(Rs, 1.0/c_bar_par);
  /*
  double extent=5.;
  double f_cutoff=sqrt( pow(xbar/(R_bar_x*extent),2.)+pow(ybar/(R_bar_y*extent),2.)+pow(zbar/(R_bar_z*extent),2.));
  if(f_cutoff<1.) f_cutoff=1-f_cutoff*f_cutoff;
  else f_cutoff=0;*/
  lumDensity = J0_bar * pow(cosh(-Rs),-2.0);
  return lumDensity;
}



/* END OF FILE: func_peanut_dattathri.cpp ------------------------------ */
