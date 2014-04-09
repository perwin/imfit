/* FILE: func_gaussianring3d.cpp -------------------------------------- */
/* VERSION 0.2
 *
 *   Experimental function object class for an elliptical 3D ring (luminosity
 * density = Gaussian centered at radius a_ring along major axis, with width sigma 
 * and vertical exponential with scale heigh h_z); ring has intrinsic (in-plane) 
 * ellipticity ell and position angle PA_ring w.r.t. to line of nodes; system is 
 * seem with line of nodes at angle PA (w.r.t. image +y axis) and inclination inc.
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
 *     [v0.2]: 21 Oct 2012: Modified from circular to elliptical ring shape.
 *     [v0.1]: 24 Aug 2012: Created (as modification of func_exp3d.cpp).
 */





/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <gsl/gsl_errno.h>

#include "func_gaussianring3d.h"
#include "integrator.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 8;
const char  PARAM_LABELS[][20] = {"PA", "inc", "PA_ring", "ell", "J_0", "a_ring", "sigma", "h_z"};
const char  FUNCTION_NAME[] = "GaussianRing3D function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const double  INTEGRATION_MULTIPLIER = 2;

const char GaussianRing3D::className[] = "GaussianRing3D";


/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensityRing( double s, void *params );




/* ---------------- CONSTRUCTOR ---------------------------------------- */

GaussianRing3D::GaussianRing3D( )
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

  // Stuff related to GSL integration
  gsl_set_error_handler_off();
  F.function = LuminosityDensityRing;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void GaussianRing3D::Setup( double params[], int offsetIndex, double xc, double yc )
{
  double  PA_rad, inc_rad, ringPA_rad;
  
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  ringPA = params[2 + offsetIndex];
  ell = params[3 + offsetIndex];
  J_0 = params[4 + offsetIndex ];
  a_ring = params[5 + offsetIndex ];
  sigma = params[6 + offsetIndex ];
  h_z = params[7 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  
  // ring PA rotations are computed relative to +y axis; convert to +x-axis reference
  ringPA_rad = (ringPA + 90.0) * DEG2RAD;
  cosRingPA = cos(ringPA_rad);
  sinRingPA = sin(ringPA_rad);
  
  inc_rad = inclination * DEG2RAD;
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);

  twosigma_squared = 2.0 * sigma*sigma;

  // We could do this here using integLimit as a class data member, but it tends
  // to be marginally slower this way
//  integLimit = INTEGRATION_MULTIPLIER * a_ring;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double GaussianRing3D::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity;
  double  integLimit;
  double  xyParameters[13];
  int  nSubsamples;
  
  // Calculate x,y in component's (projected sky) reference frame: xp,yp
  xp = x_diff*cosPA + y_diff*sinPA;
  yp = -x_diff*sinPA + y_diff*cosPA;

  // Calculate (x,y,z)_start in component's native xyz reference frame, corresponding to
  // intersection of line-of-sight ray with projected sky frame
  x_d0 = xp;
  y_d0 = yp * cosInc;
  z_d0 = yp * sinInc;

  // Set up parameter vector for the integration (everything that stays unchanged
  // for this particular xp,yp location)
  xyParameters[0] = x_d0;
  xyParameters[1] = y_d0;
  xyParameters[2] = z_d0;
  xyParameters[3] = cosInc;
  xyParameters[4] = sinInc;
  xyParameters[5] = cosRingPA;
  xyParameters[6] = sinRingPA;
  xyParameters[7] = q;
  xyParameters[8] = J_0;
  xyParameters[9] = a_ring;
  xyParameters[10] = sigma;
  xyParameters[11] = h_z;
  xyParameters[12] = twosigma_squared;
  F.params = xyParameters;

  // integrate out to +/- integLimit, which is multiple of ring radius
  integLimit = INTEGRATION_MULTIPLIER * a_ring;
  totalIntensity = Integrate(F, -integLimit, integLimit);

  return totalIntensity;
}







/* ----------------------------- OTHER FUNCTIONS -------------------------------- */


/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 */ 
double LuminosityDensityRing( double s, void *params )
{
  double  y_d, z_d, R, deltaR, lumDensity;
  double  x_ring, y_ring, z_ring, y_ring_scaled;
  double  *paramsVect = (double *)params;
  double  x_d0 = paramsVect[0];
  double  y_d0 = paramsVect[1];
  double  z_d0 = paramsVect[2];
  double  cosInc = paramsVect[3];
  double  sinInc = paramsVect[4];
  double  cosRingPA = paramsVect[5];
  double  sinRingPA = paramsVect[6];
  double  q = paramsVect[7];
  double  J_0 = paramsVect[8];
  double  a_ring = paramsVect[9];
  double  sigma = paramsVect[10];
  double  h_z = paramsVect[11];
  double  twosigma_squared = paramsVect[12];
  
  // Given current value of s and the pre-defined parameters, determine our 
  // 3D location (x_d,y_d,z_d) [by construction, x_d = x_d0]
  y_d = y_d0 + s*sinInc;
  z_d = z_d0 - s*cosInc;
  
  // Convert 3D Cartesian coordinate to rotated x_ring,y_ring,z_ring coordinate,
  // where x_ring is along ring major axis
  x_ring = x_d0*cosRingPA + y_d*sinRingPA;
  y_ring = -x_d0*sinRingPA + y_d*cosRingPA;
  z_ring = fabs(z_d);
  
  // NOTE: FOR A CONSTANT-WIDTH RING (i.e., where sigma does *not* scale with ellipticity),
  // we could: scale *a_ring* and compute R without scaling...
  // OR: proceed as normal, but rescale sigma to account for ellipticity...
  
  // Convert x_ring,y_ring to scaled radius R
  y_ring_scaled = y_ring/q;
  R = sqrt(x_ring*x_ring + y_ring_scaled*y_ring_scaled);

  deltaR = R - a_ring;
  lumDensity = J_0 * exp(-deltaR*deltaR/twosigma_squared) * exp(-z_ring/h_z);
  return lumDensity;
}



/* END OF FILE: func_gaussianring3d.cpp -------------------------------- */
