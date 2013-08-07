/* FILE: func_flatsky.cpp ---------------------------------------------- */
/* VERSION 0.1
 *
 *   This is a derived class which provides for a constant ("flat") sky
 * background.  It has only one parameter, and no dependence on pixel
 * position whatsoever.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 24 April 2010: Created (as modification of func_gaussian.cpp).
 */

// Copyright 2010, 2013 by Peter Erwin.
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

#include "func_flatsky.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 1;
const char  PARAM_LABELS[][20] = {"I_sky"};
const char  FUNCTION_NAME[] = "Flat sky background function";

const char FlatSky::className[] = "FlatSky";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

FlatSky::FlatSky( )
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

void FlatSky::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  I_sky = params[0 + offsetIndex ];
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double FlatSky::GetValue( double x, double y )
{
  return I_sky;
}



/* END OF FILE: func_flatsky.cpp --------------------------------------- */
