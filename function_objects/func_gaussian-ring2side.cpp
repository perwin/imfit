/* FILE: gaussian-on-ring2side.cpp ------------------------------------- */
/* VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces an elliptical ring with a 2-sided Gaussian profile, which
 * allows for different Gaussian sigma for r < r_ring vs r > r_ring.
 *
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
 *     [v0.1]  16 Mar 2011: Created as modification of func_edge-on-ring2side.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_gaussian-ring2side.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "ell", "A", "R_ring", "sigma_r_in", "sigma_r_out"};
const char  FUNCTION_NAME[] = "2-sided Gaussian Ring function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char GaussianRing2Side::className[] = "GaussianRing2Side";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

GaussianRing2Side::GaussianRing2Side( )
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

void GaussianRing2Side::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  A = params[2 + offsetIndex ];
  R_ring = params[3 + offsetIndex ];   // major-axis radius of ring
  sigma_r_inner = params[4 + offsetIndex ];
  sigma_r_outer = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  twosigma_rin_squared = 2.0 * sigma_r_inner*sigma_r_inner;
  twosigma_rout_squared = 2.0 * sigma_r_outer*sigma_r_outer;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for an edge-on ring 2D function, at
// cylindrical radius r (along the major-axis of the component) and height z 
// (perpendicular to the major-axis of the component).
// NOTE: This function requires that z be *non-negative*!
// input r is assumed to always be positive
double GaussianRing2Side::CalculateIntensity( double r )
{
  double  I;
  double  r_diff = fabs(r - R_ring);
  double  twosigma_r_squared;
  
  // decide which sigma to use for each Gaussian component
  if (r < R_ring) {  // inside the ring
    twosigma_r_squared = twosigma_rin_squared;
  }
  else {  // outside the ring
    twosigma_r_squared = twosigma_rout_squared;
  }
  
  I = A * exp(-(r_diff*r_diff)/twosigma_r_squared);
  return I;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// R is *non-negative*!
double GaussianRing2Side::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  
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
        xp = x_diff*cosPA + y_diff*sinPA;
        yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
        r = sqrt(xp*xp + yp_scaled*yp_scaled);
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
// given that the current pixel is a (scaled) distance of r away from the center of the
// ring.
// For now, we don't do any subsampling (because it's easier and also because Gaussians
// usually aren't as prone to subsampling effects).
int GaussianRing2Side::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  return nSamples;
}



/* END OF FILE: gaussian-on-ring2side.cpp ------------------------------ */
