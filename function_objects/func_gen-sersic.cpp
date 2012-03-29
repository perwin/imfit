/* FILE: func_gen-sersic.cpp ------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a generalized Sersic function, with constant
 * (generalized) ellipticity and position angle.
 *
 *   Generalized ellipse is from Athanassoula et al. (1990), with parameterization
 * as in Peng et al. (2002).
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
 *   NOTE: Currently, we assume input PA is in *degrees* [and then we
 * convert it to radians] relative to +x axis.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]  7 June 2010: Created (as modification of func_sersic.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <algorithm>

#include "func_gen-sersic.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "ell", "c0", "n", "I_e", "r_e"};
const char  FUNCTION_NAME[] = "Generalized-ellipse Sersic function";
//const char  SHORT_FUNCTION_NAME[] = "Sersic_GenEllipse";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char GenSersic::className[] = "Sersic_GenEllipse";

const double  A0_M03 = 0.01945;
const double  A1_M03 = -0.8902;
const double  A2_M03 = 10.95;
const double  A3_M03 = -19.67;
const double  A4_M03 = 13.43;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

GenSersic::GenSersic( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void GenSersic::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  c0 = params[2 + offsetIndex ];
  n = params[3 + offsetIndex ];
  I_e = params[4 + offsetIndex ];
  r_e = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  // generalized ellipse exponents
  ellExp = c0 + 2.0;
  invEllExp = 1.0 / ellExp;

  n2 = n*n;
  if (n > 0.36) {
    // Ciotti & Bertin (1999) approximation, good for all n > 0.36
    bn = 2*n - 0.333333333333333 + 0.009876543209876543/n
         + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
         - 7.1510122958919723e-05/(n2*n2);
  } else {
    // MacArthur et al. (2003) approximation for n < 0.36
    bn = A0_M03 + A1_M03*n + A2_M03*n2 + A3_M03*n2*n + A4_M03*n2*n2;
  }
  invn = 1.0 / n;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a Sersic function at radius r,
// with the various parameters and derived values (n, b_n, r_e, etc.)
// pre-calculated by Setup().

double GenSersic::CalculateIntensity( double r )
{
  double  intensity;
  
  intensity = I_e * exp( -bn * (pow((r/r_e), invn) - 1.0));
  return intensity;
}


/* ---------------- PRIVATE METHOD: CalculateRadius -------------------- */
// This function calculates the equivalent radius for a generalized ellipse,
// for a coordinate system where r=0 at deltaX = deltaY = 0.

double GenSersic::CalculateRadius( double deltaX, double deltaY )
{
  double  xp, yp_scaled, powerSum;
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio;
  // convert both to absolute values
  xp = fabs(deltaX*cosPA + deltaY*sinPA);
  yp_scaled = fabs((-deltaX*sinPA + deltaY*cosPA)/q);
  // Compute radius for generalized ellipse
  powerSum = pow(xp, ellExp) + pow(yp_scaled, ellExp);
  return pow(powerSum, invEllExp);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double GenSersic::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  r, totalIntensity;
  int  nSubsamples;
  
  r = CalculateRadius(x_diff, y_diff);
  
  nSubsamples = CalculateSubsamples(r);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-picel
    double deltaSubpix = 1.0 / nSubsamples;
    double x_sub_start = x - 0.5 + 0.5*deltaSubpix;
    double y_sub_start = y - 0.5 + 0.5*deltaSubpix;
    double theSum = 0.0;
    for (int ii = 0; ii < nSubsamples; ii++) {
      double x_ii = x_sub_start + ii*deltaSubpix;
      for (int jj = 0; jj < nSubsamples; jj++) {
        double y_ii = y_sub_start + jj*deltaSubpix;
        x_diff = x_ii - x0;
        y_diff = y_ii - y0;
        r = CalculateRadius(x_diff, y_diff);
        theSum += CalculateIntensity(r);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// Sersic function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int GenSersic::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Subsampling as for basic (pure-ellipse) Sersic (func_sersic.cpp)
  if ((doSubsampling) && (r < 10.0)) {
    if ((r_e <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / r_e));
    else {
      if (r <= 4.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}



/* END OF FILE: func_gen-sersic.cpp ------------------------------------ */
