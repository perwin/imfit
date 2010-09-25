/* FILE: func1d_gaussian.cpp ------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D Gaussian function.
 * Flux-related parameters are in surface-brightness (mag/arcsec^2), but output
 * is flux (calling function -- e.g., ModelObject1D::CreateModelImage -- will
 * convert back to magnitudes for comparison with data
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 18 Jun 2010: Created (as modification of func1d_exp.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_gaussian.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 2;
const char  PARAM_LABELS[][20] = {"mu_0", "sigma"};
const char FUNCTION_NAME[] = "Gaussian-1D function";
#define CLASS_SHORT_NAME  "Gaussian-1D"


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Gaussian1D::Gaussian1D( )
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

void Gaussian1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  mu_0 = params[0 + offsetIndex ];
  sigma = params[1 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  I_0 = pow(10.0, -0.4*mu_0);
//  printf("Gaussian1D::Setup: mu_0 = %g, sigma = %g, I_0 = %g\n", mu_0, sigma, I_0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// NOTE: for compatibility with 2D functions (and base class FunctionObject), we
// include y as an input, but don't use it.
double Gaussian1D::GetValue( double x )
{
//  printf("In GetValue: x = %g, I_0 = %g, h = %g\n", x, I_0, h);
//  double  mu = -2.5 * log10(I);
  double  scaledDeltaR = fabs(x - x0) / sigma;
  double  exponent = -(scaledDeltaR * scaledDeltaR)/2.0;
  return I_0 * exp( exponent );
}



/* END OF FILE: func1d_gaussian.cpp ------------------------------------ */
