/* FILE: func_sersic.cpp ----------------------------------------------- */
/* VERSION 0.3
 *
 *   Function object class for a Sersic function, with constant
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
 *     [v0.1]: 19 Nov 2009: Created (as modification of func_exp.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_sersic.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 5;
const char  PARAM_LABELS[][20] = {"PA", "ell", "n", "I_e", "r_e"};
const char  FUNCTION_NAME[] = "Sersic function";
const double  DEG2RAD = 0.017453292519943295;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Sersic::Sersic( )
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

void Sersic::Setup( double params[], int offsetIndex, double xc, double yc )
{
//   x0 = params[0 + offsetIndex];
//   y0 = params[1 + offsetIndex];
//   PA = params[2 + offsetIndex];
//   ell = params[3 + offsetIndex];
//   n = params[4 + offsetIndex ];
//   I_e = params[5 + offsetIndex ];
//   r_e = params[6 + offsetIndex ];
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  n = params[2 + offsetIndex ];
  I_e = params[3 + offsetIndex ];
  r_e = params[4 + offsetIndex ];
  printf("func_sersic: x0 = %g, y0 = %g, PA = %g, ell = %g, n = %g, r_e = %g, I_e = %g\n",
          x0, y0, PA, ell, n, r_e, I_e);
  
  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
//  printf("func_exp: q = %g, PA_rad = %g, cosPA = %g, sinPA = %g\n",
//         q, PA_rad, cosPA, sinPA);
  n2 = n*n;
  /* The following approximation for b_n is good for all n > 0.36 */
  bn = 2*n - 0.333333333333333 + 0.009876543209876543/n
       + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
       - 0.0018028610621203215e-5/(n2*n2);
  invn = 1.0 / n;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double Sersic::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  return ( I_e * exp( -bn * (pow((r/r_e), invn) - 1.0)) );
}



/* END OF FILE: func_sersic.cpp ---------------------------------------- */
