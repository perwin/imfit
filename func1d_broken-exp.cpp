/* FILE: func1d_broken-exp.cpp ----------------------------------------- */
/* VERSION 0.2
 *
 *   Function object class for a 1-D broken-exponential function.
 * Flux-related parameters are in surface-brightness (mag/arcsec^2), but output
 * is flux.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.2]: 28 Nov 2009: Updated to new FunctionObject interface.
 *     [v0.1]: 27 Nov 2009: Created (as modification of func1d_exp.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_broken-exp.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 5;
const char  PARAM_LABELS[][20] = {"mu_0", "h_1", "h_2", "r_b", "alpha"};
const char FUNCTION_NAME[] = "Broken-Exponential-1D function";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

BrokenExponential1D::BrokenExponential1D( )
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

void BrokenExponential1D::Setup( double params[], int offsetIndex )
{
  mu_0 = params[0 + offsetIndex ];
  h_1 = params[1 + offsetIndex ];
  h_2 = params[2 + offsetIndex ];
  r_b = params[3 + offsetIndex ];
  alpha = params[4 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  I_0 = pow(10.0, -0.4*mu_0);
  exponent = (1.0/alpha) * (1.0/h_1 - 1.0/h_2);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  S = pow( (1.0 + exp(-alpha*r_b)), (-exponent) );
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double BrokenExponential1D::GetValue( double x )
{
  double  I;
  
//  printf("In GetValue: x = %g, I_0 = %g, h = %g\n", x, I_0, h);
  if ( alpha*(x - r_b) > 100.0 ) {
    // Outer-exponential approximation:
    I = I_0 * S * exp(r_b/h_2 - r_b/h_1 - x/h_2);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    I = I_0 * S * exp(-x/h_1) * pow( 1.0 + exp(alpha*(x - r_b)), exponent );
  }
  return (I);
}



/* END OF FILE: func1d_broken-exp.cpp ---------------------------------- */
