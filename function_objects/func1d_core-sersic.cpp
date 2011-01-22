/* FILE: func1d_core-sersic.cpp ---------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D Sersic function.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 22 Jan 2011: Created (as modification of func1d_sersic.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include "func1d_core-sersic.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"n", "mu_b", "r_e", "r_b", "alpha", "gamma"};
const char FUNCTION_NAME[] = "Core-Sersic-1D function";
#define CLASS_SHORT_NAME  "Core-Sersic-1D"


/* ---------------- CONSTRUCTOR ---------------------------------------- */

CoreSersic1D::CoreSersic1D( )
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

void CoreSersic1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  n = params[0 + offsetIndex ];
  mu_b = params[1 + offsetIndex ];
  r_e = params[2 + offsetIndex ];
  r_b = params[3 + offsetIndex ];
  alpha = params[4 + offsetIndex ];
  gamma = params[5 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  I_b = pow(10.0, 0.4*(ZP - mu_b));
  n2 = n*n;
  /* The following approximation for b_n is good for all
   * n > 0.36 */
  bn = 2*n - 0.333333333333333 + 0.009876543209876543/n
       + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
       - 7.1510122958919723e-05/(n2*n2);
  invn = 1.0 / n;
  Iprime = I_b * pow(2.0, -gamma/alpha) * exp( bn * pow( pow(2.0, 1.0/alpha) * r_b/r_e, (1.0/n) ));
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double CoreSersic1D::GetValue( double x )
{
  double  r = fabs(x - x0);
  double  powerlaw_part = pow( 1.0 + pow(r_b/r, alpha), gamma/alpha );
  double  exp_part = exp( -bn * pow( ( pow(r, alpha) + pow(r_b, alpha) )/pow(r_e, alpha), 1.0/(alpha*n)) );
  double  I = Iprime * powerlaw_part * exp_part;
  return (I);
}



/* END OF FILE: func1d_core-sersic.cpp --------------------------------- */
