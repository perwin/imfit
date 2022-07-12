/* FILE: func_tilted-sky-plane.cpp ---------------------------------------------- */
/* 
 *   This is a derived class which provides for a sky background modeled as an
 * inclined ("tilted") plane.
 *
 *   Inspired by the InclinedFlatSky function of Dan Prole (djampro).
 *
 */

// Copyright 2020--2022 by Peter Erwin.
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

#include "func_tilted-sky-plane.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 3;
const char  PARAM_LABELS[][20] = {"I_0", "m_x", "m_y"};
const char  PARAM_UNITS[][30] = {"counts/pixel", "", ""};
const char  FUNCTION_NAME[] = "Tilted sky-plane background function";

const char TiltedSkyPlane::className[] = "TiltedSkyPlane";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

TiltedSkyPlane::TiltedSkyPlane( )
{

  nParams = N_PARAMS;
  
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up vectors of parameter labels and units
  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }
  parameterUnitsExist = true;
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void TiltedSkyPlane::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  I_0 = params[0 + offsetIndex ];
  m_x = params[1 + offsetIndex ];   // slope in x-direction
  m_y = params[2 + offsetIndex ];   // slope in y-direction
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double TiltedSkyPlane::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;

  return I_0 + m_x*x_diff + m_y*y_diff;
}


/* ---------------- PUBLIC METHOD: IsBackground ------------------------ */

bool TiltedSkyPlane::IsBackground( )
{
  return true;
}



/* END OF FILE: func_tilted-sky-plane.cpp ------------------------------ */
