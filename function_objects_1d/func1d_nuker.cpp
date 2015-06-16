/* FILE: func1d_nuker.cpp ---------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D Nuker-law function.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 23 Oct 2013: Created (as modification of func1d_core-sersic.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include "func1d_nuker.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 5;
const char  PARAM_LABELS[][20] = {"alpha", "beta", "gamma", "r_b", "mu_b"};
const char FUNCTION_NAME[] = "NukerLaw-1D function";
#define CLASS_SHORT_NAME  "NukerLaw-1D"

// The following is the minimum allowable value of r = |x - x0|, meant to
// avoid blowups due to the fact that the inner, power-law part of the
// Nuker-law function becomes infinite at r = 0.
const double  R_MIN = 0.01;

const char NukerLaw1D::className[] = CLASS_SHORT_NAME;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

NukerLaw1D::NukerLaw1D( )
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

void NukerLaw1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  alpha = params[0 + offsetIndex ];
  beta = params[1 + offsetIndex ];
  gamma = params[2 + offsetIndex ];
  r_b = params[3 + offsetIndex ];
  mu_b = params[4 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
    // code from IDL function
//   I_b = 10.0^(-mu_b*0.4)
//   I1 = I_b*2.0^( (beta - gamma)/alpha ) * (r_b/x)^gamma

  double  I_b = pow(10.0, 0.4*(ZP - mu_b));
  exponent = (beta - gamma)/alpha;
  Iprime = I_b * pow(2.0, exponent);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double NukerLaw1D::GetValue( double x )
{
  double  r = fabs(x - x0);
  // kludge to handle cases when r is very close to zero:
  if (r < R_MIN)
    r = R_MIN;

    // code from IDL function
//   I1 = I_b*2.0^( (beta - gamma)/alpha ) * (r_b/x)^gamma
//   I2 = (1.0 + (x/r_b)^alpha)^( (gamma - beta)/alpha )

  double  I1 = Iprime * pow(r_b/r, gamma);
  double  I2 = pow((1.0 + pow(r/r_b, alpha)), -exponent);
  return (I1*I2);
}



/* END OF FILE: func1d_nuker.cpp --------------------------------------- */
