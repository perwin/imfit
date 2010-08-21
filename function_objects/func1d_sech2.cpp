/* FILE: func1d_sech2.cpp ---------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D sech^2 function (output in magnitudes
 * per sq.arcsec).
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 20 Aug 2010: Created (as modification of func1d_exp.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_sech2.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 2;
const char  PARAM_LABELS[][20] = {"mu_0", "h"};
const char FUNCTION_NAME[] = "Sech^2-1D function";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Sech21D::Sech21D( )
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

void Sech21D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  mu_0 = params[0 + offsetIndex ];
  h = params[1 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  I_0 = pow(10.0, -0.4*mu_0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Sech21D::GetValue( double x )
{
  double  r = fabs(x - x0);
  double  sech = (1.0 / cosh(r/h));

  return (I_0 * sech*sech);
}



/* END OF FILE: func1d_sech2.cpp --------------------------------------- */
