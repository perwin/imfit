/* FILE: func_incflatsky.cpp ---------------------------------------------- */
/* 
 *   This is a derived class which provides for an inclined flat sky
 * background plane.  It has only three parameters, the central brightness,
 * and the spatial derivatives of its brightness in the x and y directions.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 24 April 2010: Created (as modification of func_flatsky.cpp).
 */

// Copyright 2010--2016 by Peter Erwin.
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

#include "func_incflatsky.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 3;
const char  PARAM_LABELS[][20] = {"I0", "dIdx", "dIdy"};
const char  FUNCTION_NAME[] = "Inclined flat sky background function";

const char InclinedFlatSky::className[] = "InclinedFlatSky";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

InclinedFlatSky::InclinedFlatSky( )
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

void InclinedFlatSky::Setup( double params[], int offsetIndex,
                             double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  I0 = params[0 + offsetIndex ];
  dIdx = params[1 + offsetIndex ];
  dIdy = params[2 + offsetIndex ];
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double InclinedFlatSky::GetValue( double x, double y )
{
  return I0 + (x-x0)*dIdx + (x-y0)*dIdy
}


/* ---------------- PUBLIC METHOD: IsBackground ------------------------ */

bool InclinedFlatSky::IsBackground( )
{
  return true;
}


/* END OF FILE: func_incflatsky.cpp --------------------------------------- */
