/* FILE: func1d_spline.cpp --------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a 1-D cubic spline interpolation function (output in magnitudes
 * per sq.arcsec).
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 26 Oct 2016: Created (as modification of func1d_exp.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include <gsl/gsl_spline.h>

#include "func1d_spline.h"
#include "utilities_pub.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
// maximum number of data points (always >= 2)
const char  PARAM_LABELS[][20] = {"I_0", "r_1", "I_1", "r_2", "I_2", "r_3", "I_3"};
const char FUNCTION_NAME[] = "Spline-1D function";
#define CLASS_SHORT_NAME  "Spline-1D"

const char Spline1D::className[] = CLASS_SHORT_NAME;



/* ---------------- CONSTRUCTOR ---------------------------------------- */

Spline1D::Spline1D( )
{
  string  paramName;
  
  functionName = FUNCTION_NAME;
  shortFunctionName = CLASS_SHORT_NAME;
  
  maxInterpPoints = MAX_SPLINE_POINTS;   // may be modified by SetExtraParams()
  nParams = 2*maxInterpPoints - 1;  // r,I values for each interpolation point, with first
                                    // point r = 0 by definition (and so not an input parameter)
                                    //  may be modified by SetExtraParams()

  // Set up the (maximum-sized) vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  splineFuncAllocated = false;
  splineCacheAllocated = false;
  for (int i = 0; i < MAX_SPLINE_POINTS; i++) {
    xInterp[i] = 0.0;
    yInterp[i] = 0.0;
  }
  splineCache = gsl_interp_accel_alloc();
  splineCacheAllocated = true;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

Spline1D::~Spline1D( )
{
  // GSL spline cleanup:
  if (splineFuncAllocated)
    gsl_spline_free(splineFunc);
  if (splineCacheAllocated)
    gsl_interp_accel_free(splineCache);
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Spline1D::Setup( double params[], int offsetIndex, double xc )
{
  int  paramIndex;
  double  currentRVal;
  
  x0 = xc;
  I_0 = params[0 + offsetIndex ];
  r_1 = params[1 + offsetIndex ];
  I_1 = params[2 + offsetIndex ];
  
  // Set up "data vector" for interpolation, and determine how many
  // valid interpolation points we actually have
  // r_n < 0 ==> this and any points further out are to be ignored
  xInterp[0] = 0.0;
  yInterp[0] = I_0;
  xInterp[1] = r_1;
  yInterp[1] = I_1;
  nInterpPoints = maxInterpPoints;
  for (int i = 2; i < maxInterpPoints; i++) {
    paramIndex = 2*i - 1 + offsetIndex;
    currentRVal = params[paramIndex];
//     printf("i = %d: paramIndex = %d, currentRVal = %f\n", i, paramIndex, currentRVal);
    if (currentRVal < 0) {
      nInterpPoints = i;
      break;
    } else {
      xInterp[i] = currentRVal;
      yInterp[i] = params[paramIndex + 1];
    }
  }

  gsl_interp_accel_reset(splineCache);
  if (splineFuncAllocated)
    gsl_spline_free(splineFunc);
  splineFunc = gsl_spline_alloc(gsl_interp_cspline, nInterpPoints);
  splineFuncAllocated = true;

//   printf("Calling gsl_spline_init (size = %d) with nInterpPoints = %d\n", (int)splineFunc->size, nInterpPoints);
//   printf("   using following data points: ");
//   for (int i = 0; i < nInterpPoints; i++)
//     printf("(i,r,y) = (%d,%.3f,%.3f)", i, xInterp[i], yInterp[i]);
  gsl_spline_init(splineFunc, xInterp, yInterp, nInterpPoints);
}


/* ---------------- PUBLIC METHOD: HasExtraParams ---------------------- */

bool Spline1D::HasExtraParams( )
{
  return true;
}


/* ---------------- PUBLIC METHOD: SetExtraParams ---------------------- */
// Returns -1 if map is empty, 0 if map is not empty but no valid parameter
// name is found. If map has valid parameter name, returns 1 if parameter
// value is OK, -3 if not.
int Spline1D::SetExtraParams( map<string,string>& inputMap )
{
  // check for empty map
  if (inputMap.empty()) {
    printf("   Spline1D::SetExtraParams: input map is empty!\n");
    return -1;
  }
  // only one possible parameter for this function, so no need to loop
  map<string,string>::iterator iter;
  for( iter = inputMap.begin(); iter != inputMap.end(); iter++) {
    if (iter->first == "N") {
      if (NotANumber(iter->second.c_str(), 0, kPosInt))
        return -3;
      else {
        maxInterpPoints = atoi(iter->second.c_str());
        if (maxInterpPoints > MAX_SPLINE_POINTS)
          return -3;
        printf("   Spline1D::SetExtraParams -- setting maxInterpPoints = %d\n", maxInterpPoints);
        extraParamsSet = true;
        nParams = 2*maxInterpPoints - 1;
        return 1;
      }
    }
  }
  return 0;
}



/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
double Spline1D::GetValue( double x )
{
  double  r = fabs(x - x0);
  return (gsl_spline_eval(splineFunc, r, splineCache));
}



/* END OF FILE: func1d_spline.cpp -------------------------------------- */
