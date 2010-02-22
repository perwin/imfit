/* FILE: func_exp.cpp -------------------------------------------------- */
/* VERSION 0.3
 *
 *   Function object class for an exponential function, with constant
 * ellipticity and position angle (pure elliptical, not generalized).
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
 *     [v0.3]: 21 Jan 2010: Modified to treat x0,y0 as separate inputs.
 *     [v0.2]: 28 Nov 2009: Updated to new FunctionObject interface.
 *     [v0.1]: 18 Nov 2009: Created (as modification of func_gauss.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_exp.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 4;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_0", "h"};
const char FUNCTION_NAME[] = "Exponential function";
const double  DEG2RAD = 0.017453292519943295;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Exponential::Exponential( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Exponential::Setup( double params[], int offsetIndex, double xc, double yc )
{
//   x0 = params[0 + offsetIndex];
//   y0 = params[1 + offsetIndex];
//   PA = params[2 + offsetIndex];
//   ell = params[3 + offsetIndex];
//   I_0 = params[4 + offsetIndex ];
//   h = params[5 + offsetIndex ];
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_0 = params[2 + offsetIndex ];
  h = params[3 + offsetIndex ];
// #ifdef DEBUG
//   printf("func_exp: x0 = %g, y0 = %g, PA = %g, ell = %g, I_0 = %g, h = %g\n",
//           x0, y0, PA, ell, I_0, h);
// #endif  
  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
//  printf("func_exp: q = %g, PA_rad = %g, cosPA = %g, sinPA = %g\n",
//         q, PA_rad, cosPA, sinPA);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double Exponential::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  return I_0 * exp(-r/h);
}



/* END OF FILE: func_exp.cpp ------------------------------------------- */
