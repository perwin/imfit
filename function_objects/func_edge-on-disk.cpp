/* FILE: func_edge-on-disk.cpp ----------------------------------------- */
/* VERSION 0.2
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces a generalized edge-on exponential disk, using Bessel-function
 * solution (van der Kruit & Searle 1981) for radial profile and generalized
 * sech function (van der Kruit 1988) for vertical profile.  My version of this
 * looks like the following (Sigma = surface brightness in counts/pixel):
 *
 *      Sigma(r,z) = Sigma(0,0) * (r/h) * K_1(r/h) * sech^alpha(r/(alpha*z0))
 *
 *    with Sigma(0,0) = 2 * h * L_0
 *
 *   This should corresponds to a face-on disk with:
 *      Sigma(r) = mu_0 * exp(-r/h)
 *   where mu_0 = 2 * z0 * L_0
 *   although this is guaranteed to be true only for alpha = 2 [i.e., that is
 * the most direct match to van der Kruit & Searle 1981]
 *
 *   Note that for the case of a sech^2 vertical profile (alpha = 2), our z0
 * is 1/2 of the usual z0 in e.g. van der Kruit & Searle (1981).
 *
 *   The main notational change from van der Kruit's original version is that
 * our alpha = his n/2
 *
 *   TO-DO:
 *      [] SPECULATIVE: have an option where alpha = 0 invokes calculation of
 * z-profile with exact exponential calculation (note that for continuity,
 * we need to keep using z0 and have it correspond to the correct exponential
 * calculation; see van der Kruit 1988).
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
 *     [v0.2]  25 Sept: Corrected calculation of r to be actual cylindrical
 * radius; changed xp and z_perp to R and z; updated subsample calculations
 * to incorporate z-dependence as well as r-dependence.
 *     [v0.1]  21 Sept 2010: Created as modification of func_broken-exp2d.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_edge-on-disk.h"
#include "gsl/gsl_sf_bessel.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 5;
const char  PARAM_LABELS[][20] = {"PA", "L_0", "h", "alpha", "z_0"};
const char  FUNCTION_NAME[] = "Edge-on Disk function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char EdgeOnDisk::className[] = "EdgeOnDisk";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

EdgeOnDisk::EdgeOnDisk( )
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

void EdgeOnDisk::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  L_0 = params[1 + offsetIndex ];
  h = params[2 + offsetIndex ];
  alpha = params[3 + offsetIndex ];
  z_0 = params[4 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  scaledZ0 = alpha*z_0;
  Sigma_00 = 2.0*h*L_0;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for an edge-on disk 2D function, at
// cylindrical radius r (along the major-axis of the component) and height z 
// (perpendicular to the major-axis of the component).
// NOTE: This function requires that both r and z be *non-negative*!
double EdgeOnDisk::CalculateIntensity( double r, double z )
{
  double  verticalScaling, sech, I_radial;
  
  // special case of r = 0, to short-circuit 0*infty problem
  if (r == 0.0)
    I_radial = Sigma_00;
  else {
    double  scaledR = r/h;
    I_radial = Sigma_00 * scaledR * gsl_sf_bessel_K1(scaledR);
  }

  sech = 1.0 / cosh(z/scaledZ0);
  verticalScaling = pow(sech, alpha);
  return I_radial * verticalScaling;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// NOTE: x and y are orthognal image coordinates; R and z are orthogonal
// coordinates *in the component reference frame* (corresponding to xp and
// yp in other, non-edge-on function objects).
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// R and z are *non-negative*!
double EdgeOnDisk::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  R, z, totalIntensity;
  int  nSubsamples;
  
  // Calculate R,z (= x,y in component reference frame)
  R = fabs(x_diff*cosPA + y_diff*sinPA);    // "R" is x in the component reference frame
  z = fabs(-x_diff*sinPA + y_diff*cosPA);   // "z" is y in the component reference frame

  nSubsamples = CalculateSubsamples(R, z);
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
        R = fabs(x_diff*cosPA + y_diff*sinPA);
        z = fabs(-x_diff*sinPA + y_diff*cosPA);
        theSum += CalculateIntensity(R, z);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(R, z);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// r=0 line and/or z away from the disk plane.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int EdgeOnDisk::CalculateSubsamples( double r, double z )
{
  int  nSamples = 1;
  int  nSr, nSz;
  double  R_abs = fabs(r);
  double  z_abs = fabs(z);
  
  // based on standard exponential-function subsampling
  if ( (doSubsampling) && ((R_abs < 10.0) || (z_abs < 10.0)) ) {
    if ( ((h <= 1.0) && (R_abs <= 1.0)) || ((z_0 <= 1.0) && (z_abs <= 1.0)) ) {
      nSr = min(100, (int)(2 * SUBSAMPLE_R / h));
      nSz = min(100, (int)(2 * SUBSAMPLE_R / z_0));
      nSamples = max(nSr, nSz);
    }
    else {
      if ((R_abs <= 3.0) || (z_abs <= 3.0))
        nSamples = 2 * SUBSAMPLE_R;
      else {
        nSr = min(100, (int)(2 * SUBSAMPLE_R / R_abs));
        nSz = min(100, (int)(2 * SUBSAMPLE_R / z_abs));
        nSamples = max(nSr, nSz);
      }
    }
  }
  return nSamples;
}



/* END OF FILE: func_edge-on-disk.cpp ---------------------------------- */
