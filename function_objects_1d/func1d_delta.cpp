/* FILE: func1d_delta.cpp ---------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D delta function (output in counts).
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 16 Aug 2010: Created (as modification of func1d_exp.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
//#include <math.h>
//#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_delta.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 1;
const char  PARAM_LABELS[][20] = {"I"};
const char FUNCTION_NAME[] = "Delta-1D function";
#define CLASS_SHORT_NAME  "Delta-1D"

const char Delta1D::className[] = CLASS_SHORT_NAME;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Delta1D::Delta1D( )
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

void Delta1D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  I = params[0 + offsetIndex ];
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Delta1D::GetValue( double x )
{
  if (x == x0)
    return I;
  else
    return 0.0;
}



/* END OF FILE: func1d_delta.cpp --------------------------------------- */
