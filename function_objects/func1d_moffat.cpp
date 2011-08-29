/* FILE: func1d_moffat.cpp --------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D Moffat function.
 *   This will generate a normalized profile (total flux = 1.0).
 *
 *   NOTES: This function returns the intensity using the Moffat function:
 *      I(r) = I_0 / [1 + (r/alpha)^2]^beta
 *
 *   User inputs are I_0, FWHM and beta; the scale radius alpha is derived from:
 *      alpha = FWHM / (2 [2^(1/beta) - 1]^0.5 );
 *   i.e.,
 *      FWHM = 2 alpha [2^(1/beta) - 1]^0.5
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 12 Aug 2010: Created (as modification of func1d_moffat.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_moffat.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 3;
const char  PARAM_LABELS[][20] = {"I_0", "fwhm", "beta"};
const char FUNCTION_NAME[] = "Moffat-1D function";
#define CLASS_SHORT_NAME  "Moffat-1D"

const char Moffat1D::className[] = CLASS_SHORT_NAME;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Moffat1D::Moffat1D( )
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

void Moffat1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  I_0 = params[0 + offsetIndex ];
  fwhm = params[1 + offsetIndex ];
  beta = params[2 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  // compute alpha:
  double  exponent = pow(2.0, 1.0/beta);
  alpha = 0.5*fwhm/sqrt(exponent - 1.0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Moffat1D::GetValue( double x )
{
  double  scaledR, denominator;
  
  scaledR = fabs(x - x0) / alpha;
  denominator = pow((1.0 + scaledR*scaledR), beta);
  return (I_0 / denominator);
}



/* END OF FILE: func1d_moffat.cpp -------------------------------------- */
