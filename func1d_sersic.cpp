/* FILE: func1d_sersic.cpp --------------------------------------------- */
/* VERSION 0.2
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
 *     [v0.2]: 28 Nov 2009: Updated to new FunctionObject interface.
 *     [v0.1]: 27 Nov 2009: Created (as modification of func_sersic.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include "func1d_sersic.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 3;
const char  PARAM_LABELS[][20] = {"n", "mu_e", "r_e"};
const char FUNCTION_NAME[] = "Sersic-1D function";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Sersic1D::Sersic1D( )
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

void Sersic1D::Setup( double params[], int offsetIndex )
{
  n = params[0 + offsetIndex ];
  mu_e = params[1 + offsetIndex ];
  r_e = params[2 + offsetIndex ];
//  printf("func_Sersic1D: x0 = %g, y0 = %g, PA = %g, ell = %g, n = %g, r_e = %g, I_e = %g\n",
//          x0, y0, PA, ell, n, r_e, I_e);
  
  // pre-compute useful things for this round of invoking the function
  I_e = pow(10.0, -0.4*mu_e);
  n2 = n*n;
  /* The following approximation for b_n is good for all
   * n > 0.36 */
  bn = 2*n - 0.333333333333333 + 0.009876543209876543/n
       + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
       - 0.0018028610621203215e-5/(n2*n2);
  invn = 1.0 / n;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double Sersic1D::GetValue( double x )
{
  double  I = I_e * exp( -bn * (pow((x/r_e), 1.0/n) - 1.0) );
//  double  mu = -2.5 * log10(I);
  // printf("\tI_e = %g, n2 = %g, bn = %g, I = %g, mu = %g\n", I_e, n2, bn, I, mu);
  return (I);
}



/* END OF FILE: func1d_sersic.cpp -------------------------------------- */
