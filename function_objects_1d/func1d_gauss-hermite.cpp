/* FILE: func1d_gauss-hermite.cpp --------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D Gauss-Hermite function.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_gauss-hermite.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 4;
const char  PARAM_LABELS[][20] = {"A", "sigma", "h3", "h4"};
const char FUNCTION_NAME[] = "GaussHermite-1D function";
#define CLASS_SHORT_NAME  "GaussHermite-1D"

const char GaussHermite1D::className[] = CLASS_SHORT_NAME;
const double PI = 3.14159265358979;



/* ---------------- CONSTRUCTOR ---------------------------------------- */

GaussHermite1D::GaussHermite1D( )
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

void GaussHermite1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  A = params[0 + offsetIndex ];
  sigma = params[1 + offsetIndex ];
  h3 = params[2 + offsetIndex ];
  h4 = params[3 + offsetIndex ];
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double GaussHermite1D::GetValue( double x )
{
  double  dx_div_sigma = (x - x0) / fabs(sigma);
  double  dx_div_sigma_squared = dx_div_sigma * dx_div_sigma;
  double  ghScaling = 1.0;
  double  gaussian, tmp;
  
  if (h3 != 0.0) {
  	tmp = (2*sqrt(2)*dx_div_sigma*dx_div_sigma_squared - 3*sqrt(2)*dx_div_sigma) / sqrt(6);
  	ghScaling += tmp;
  }
  if (h4 != 0.0) {
  	tmp = (4*dx_div_sigma_squared*dx_div_sigma_squared - 12*dx_div_sigma_squared + 3) / sqrt(24);
  	ghScaling += tmp;
  }
  gaussian = A * exp(-0.5*dx_div_sigma_squared);
  return ghScaling * gaussian;
}

//   double  scaledDeltaR = fabs(x - x0) / sigma;
//   double  exponent = -(scaledDeltaR * scaledDeltaR)/2.0;
//   return I_0 * exp( exponent );


/* END OF FILE: func1d_gauss-hermite.cpp ------------------------------- */
