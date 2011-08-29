/* FILE: func1d_gaussian2side.cpp -------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for an asymmetric 1-D Gaussian function.
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
 *     [v0.1]: 17 Mar 2011: Created (as modification of func1d_gaussian.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_gaussian2side.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 3;
const char  PARAM_LABELS[][20] = {"mu_0", "sigma_left", "sigma_right"};
const char FUNCTION_NAME[] = "Asymmetric Gaussian-1D function";
#define CLASS_SHORT_NAME  "Gaussian2Side-1D"


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Gaussian2Side1D::Gaussian2Side1D( )
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

void Gaussian2Side1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  mu_0 = params[0 + offsetIndex ];
  sigma_left = params[1 + offsetIndex ];
  sigma_right = params[2 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  I_0 = pow(10.0, 0.4*(ZP - mu_0));
//  printf("Gaussian2Side1D::Setup: mu_0 = %g, sigma = %g, I_0 = %g\n", mu_0, sigma, I_0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Gaussian2Side1D::GetValue( double x )
{
//  printf("In GetValue: x = %g, I_0 = %g, h = %g\n", x, I_0, h);
//  double  mu = -2.5 * log10(I);
  double  deltaR = x - x0;
  double  scaledDeltaR, exponent;
  
  if (deltaR < 0)
    scaledDeltaR = fabs(deltaR) / sigma_left;
  else
    scaledDeltaR = fabs(deltaR) / sigma_right;
  exponent = -(scaledDeltaR * scaledDeltaR)/2.0;
  return I_0 * exp( exponent );
}



/* END OF FILE: func1d_gaussian2side.cpp ------------------------------- */
