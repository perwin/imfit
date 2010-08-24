/* FILE: func_gaussian.cpp --------------------------------------------- */
/* VERSION 0.3
 *
 *   This is the base class for the various function object classes.
 *   It really shouldn't be instantiated by itself.
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
 *   MODIFICATION HISTORY:
 *     [v0.3]: 21 Jan 2010: Modified to treat x0,y0 as separate inputs.
 *     [v0.2]: 28 Nov 2009: Updated to new FunctionObject interface.
 *     [v0.01]: 13--15 Nov 2009: Created (as modification of nonlinfit2's
 *   function_object class).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_gaussian.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 2;
const char  PARAM_LABELS[][20] = {"A", "sigma"};
const char  FUNCTION_NAME[] = "Circular Gaussian function";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Gaussian::Gaussian( )
{
  string  paramName;
  nParams = N_PARAMS;
  
  functionName = FUNCTION_NAME;
  shortFunctionName = CLASS_SHORT_NAME;   // defined in header file

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Gaussian::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  A = params[0 + offsetIndex ];
  sigma = params[1 + offsetIndex ];
// #ifdef DEBUG
//   printf("func_gaussian: x0 = %g, y0 = %g, A = %g, sigma = %g\n",
//           x0, y0, A, sigma);
// #endif
  
  // pre-compute useful things for this round of invoking the function
  twosigma_squared = 2.0 * sigma*sigma;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double Gaussian::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  r_squared = x_diff*x_diff + y_diff*y_diff;
  return A * exp(-r_squared/twosigma_squared);
}



/* END OF FILE: func_gaussian.cpp -------------------------------------- */
