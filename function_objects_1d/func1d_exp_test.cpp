/* FILE: func1d_exp_test.cpp ------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D exponential function (output in counts/pixel)
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
#include <stdlib.h>
#include <map>
#include <string>

#include "utilities_pub.h"
#include "func1d_exp_test.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 2;
const char  PARAM_LABELS[][20] = {"I_0", "h"};
const char FUNCTION_NAME[] = "Exponential-1D function (testing)";
#define CLASS_SHORT_NAME  "Exponential-1D_test"

const char Exponential1D_test::className[] = CLASS_SHORT_NAME;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

Exponential1D_test::Exponential1D_test( )
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
  
  // Default values for extra params:
  floorValue = 0.0;
}


/* ---------------- PUBLIC METHOD: HasExtraParams ---------------------- */

bool Exponential1D_test::HasExtraParams( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: SetExtraParams ---------------------- */

int Exponential1D_test::SetExtraParams( map<string,string>& inputMap )
{
  // check for empty map
  if (inputMap.empty()) {
    printf("   Exponentia1D_test::SetExtraParams: input map is empty!\n");
    return -1;
  }
  // only one possible parameter for this function, so no need to loop
  map<string,string>::iterator iter;
  for( iter = inputMap.begin(); iter != inputMap.end(); iter++) {
    if (iter->first == "floor") {
      if (IsNumeric(iter->second.c_str())) {
        floorValue = atof(iter->second.c_str());
        printf("   Exponential1D_test::SetExtraParams -- setting floor = %f\n", floorValue);
        extraParamsSet = true;
        return 1;
      } else
        return -3;
    }
  }
  return 0;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Exponential1D_test::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
//  mu_0 = params[0 + offsetIndex ];
  I_0 = params[0 + offsetIndex ];
  h = params[1 + offsetIndex ];
//  printf("func_exp: x0 = %g, y0 = %g, PA = %g, ell = %g, I_0 = %g, h = %g\n",
//          x0, y0, PA, ell, I_0, h);
  
  // pre-compute useful things for this round of invoking the function
//  I_0 = pow(10.0, -0.4*mu_0);
//  I_0 = pow(10.0, 0.4*(ZP - mu_0));
//  printf("Exponential1D::Setup: mu_0 = %g, h = %g, I_0 = %g\n", mu_0, h, I_0);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Exponential1D_test::GetValue( double x )
{
//  printf("In GetValue: x = %g, I_0 = %g, h = %g\n", x, I_0, h);
//  double  mu = -2.5 * log10(I);
  double  r = fabs(x - x0);
  return (I_0 * exp(-r/h) + floorValue);
}



/* END OF FILE: func1d_exp.cpp ----------------------------------------- */
