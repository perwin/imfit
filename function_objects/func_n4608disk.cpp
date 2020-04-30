/* FILE: func_n4608disk.cpp -------------------------------------------- */
/*
 *   Function object class for a component that combines an instance
 * of BrokenExponential and an instance of GaussianRingAz, for use in
 * modeling WFC3-UVIS images of NGC 4608.
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
 *     [v0.1]  30 Apr 2020: Created as modification of func_brokenexp.cpp.
 */

// Copyright 2010--2020 by Peter Erwin.
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

#include "func_n4608disk.h"
#include "func_broken-exp.h"
#include "func_gaussian-ring-az.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 13;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_0", "h1", "h2", "r_break", "alpha",
				"PA_ring", "ell_ring", "A_maj", "A_min", "R_ring", "sigma_r"};
const char  FUNCTION_NAME[] = "NGC4608 main-disk function (BrokenExponential + GaussianRingAz)";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char N4608Disk::className[] = "N4608Disk";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

N4608Disk::N4608Disk( )
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


/* ---------------- PUBLIC METHOD: SetSubsampling ---------------------- */
/// Turn pixel subsampling on or off (true = on, false = off).
void N4608Disk::SetSubsampling( bool subsampleFlag )
{
  doSubsampling = subsampleFlag;
  // specify subsampling for component functions
  funcBrokenExp.SetSubsampling(subsampleFlag);
  funcGaussianRingAz.SetSubsampling(subsampleFlag);
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void N4608Disk::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_0 = params[2 + offsetIndex ];
  h1 = params[3 + offsetIndex ];
  h2 = params[4 + offsetIndex ];
  r_b = params[5 + offsetIndex ];
  alpha = params[6 + offsetIndex ];
  // Initial versions of parameters for GaussianRingAz
  PA_ring = params[7 + offsetIndex];
  ell_ring = params[8 + offsetIndex ];
  A_maj_rel = params[9 + offsetIndex ];
  A_min_rel = params[10 + offsetIndex ];
  R_ring = params[11 + offsetIndex ];
  sigma_r = params[12 + offsetIndex ];

  // And now pass the initial params to the appropriate sub-functions

  // BrokenExp is simple, since its parameters start at offsetIndex and
  // don't need modifications
  funcBrokenExp.Setup(params, offsetIndex, x0, y0);
  
  // GaussianRingAz is a bit more complicated
  // compute absolute intensities for GaussianRingAz
  double  A_maj = A_maj_rel * I_0;
  double  A_min = A_min_rel * I_0;
  ringParams[0] = PA_ring;
  ringParams[1] = ell_ring;
  ringParams[2] = A_maj;
  ringParams[3] = A_min;
  ringParams[4] = R_ring;
  ringParams[5] = sigma_r;
  funcGaussianRingAz.Setup(ringParams, 0, x0, y0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double N4608Disk::GetValue( double x, double y )
{
  double  I_BrokenExp = funcBrokenExp.GetValue(x, y);
  double  I_GaussianRingAz = funcGaussianRingAz.GetValue(x, y);

  return I_BrokenExp + I_GaussianRingAz;
}


/* END OF FILE: func_n4608disk.cpp ------------------------------------- */
