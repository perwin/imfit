/* FILE: func1d_exp.cpp ------------------------------------------------ */
/* VERSION 0.2
 *
 *   Function object class for a 1-D exponential function (output in magnitudes
 * per sq.arcsec).
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.2]: 28 Nov 2009: Updated to new FunctionObject interface.
 *     [v0.1]: 27 Nov 2009: Created (as modification of func_exp.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_exp.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 2;
const char  PARAM_LABELS[][20] = {"mu_0", "h"};
const char FUNCTION_NAME[] = "Exponential-1D function";
#define CLASS_SHORT_NAME  "Exponential-1D"


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Exponential1D::Exponential1D( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = CLASS_SHORT_NAME;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Exponential1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  mu_0 = params[0 + offsetIndex ];
  h = params[1 + offsetIndex ];
//  printf("func_exp: x0 = %g, y0 = %g, PA = %g, ell = %g, I_0 = %g, h = %g\n",
//          x0, y0, PA, ell, I_0, h);
  
  // pre-compute useful things for this round of invoking the function
//  I_0 = pow(10.0, -0.4*mu_0);
  I_0 = pow(10.0, 0.4*(ZP - mu_0));
//  printf("Exponential1D::Setup: mu_0 = %g, h = %g, I_0 = %g\n", mu_0, h, I_0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Exponential1D::GetValue( double x )
{
//  printf("In GetValue: x = %g, I_0 = %g, h = %g\n", x, I_0, h);
//  double  mu = -2.5 * log10(I);
  double  r = fabs(x - x0);
  return (I_0 * exp(-r/h));
}



/* END OF FILE: func1d_exp.cpp ----------------------------------------- */
