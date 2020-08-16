/* FILE: func_simple-checkerboard.cpp ---------------------------------- */
/* 
 *   This is a derived class which provides for a constant ("flat") sky
 * background.  It has only one parameter, and no dependence on pixel
 * position whatsoever.
 *
 */

// Copyright 2020by Peter Erwin.
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

#include "func_simple-checkerboard.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 1;
const char  PARAM_LABELS[][20] = {"I_pos"};
const char  FUNCTION_NAME[] = "Simple checkerboard background function";

const char SimpleCheckerboard::className[] = "SimpleCheckerboard";

const int  SPACING = 10;  // size of individual square


/* ---------------- CONSTRUCTOR ---------------------------------------- */

SimpleCheckerboard::SimpleCheckerboard( )
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

void SimpleCheckerboard::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  I_pos = params[0 + offsetIndex ];
  
  step = 2*SPACING;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double SimpleCheckerboard::GetValue( double x, double y )
{
  double  I;
  int  x_int = (int)floor(x - x0 + 0.5);
  int  y_int = (int)floor(y - y0 + 0.5);
  int  x_quotient = x_int / SPACING;
  int  y_quotient = y_int / SPACING;
  
  I = 0.0;
  if ((y_quotient % 2) == 0) {  // even (C-counting) row (row 0, 2, etc.)
    if ((x_quotient % 2) == 0)  // even (C-counting) column
      I = I_pos;
  } else {   // odd (C-counting) row
    if ((x_quotient % 2) == 1)  // odd (C-counting) column
      I = I_pos;
  }
  return I;
}


/* ---------------- PUBLIC METHOD: IsBackground ------------------------ */

bool SimpleCheckerboard::IsBackground( )
{
  return true;
}



/* END OF FILE: func_simple-checkerboard.cpp --------------------------- */
