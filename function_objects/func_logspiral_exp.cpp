/* FILE: func_logspiral_exp.cpp ---------------------------------------- */
/* VERSION 0.1
 *
 *   Class (derived from FunctionObject; function_object.h) which produces a simple
 * inner-truncated logarithmic spiral. The basic spiral pattern is that specified by 
 * Eqn. 8 of Junqueira+2013, with the radial amplitude specified by an exponential 
 * for R > R_max and by (R/R_max * a Gaussian * the exponential) for R < R_max.
 *
 *
 *   Eqn. 8 of Junqueira+2013:
 *   Sigma(r,phi) = Sigma_rad(r) * exp( (-R^2/sigma_az^2) * (1 - cos(m*phi - f_m(R))) )
 *
 *      where:
 *         Sigma_rad(r) describes radial variation of spiral amplitude
 *           (in this version, an inner-truncated exponential -- see below)
 *         and the rest of the RHS is an azimuthal Gaussian with m peaks:
 *           (1 - cos(...)) is a term which varies sinusoidally from 0 to 2 pi, with
 *              m peaks (where value = 2)
 *           f_m(R) is "shape function", which describes the logarithmic spiral;
 *              it acts as a phase-angle offset, whose angular position shifts with
 *              radius following. For a logarithmic spiral:
 *              f_m(R) = (m/tan i) * ln(R/R_i) [+ gamma]
 *                 where i = pitch angle, R_i = "point where the spiral crosses the
 *                 coordinate x = 0" (and radius of rings if spiral were infinitely wound)
 *                 and gamma is azimuthal-angle offset for spiral pattern
 *
 *         Sigma_rad(r) = a function which is a simple exponential for r > R_max and
 *            a damped exponential at smaller radius, where the damping is the product
 *              of the scaled radius r/R_max and a Gaussian with dispersion sigma_r:
 *
 *              I_0 * exp(-r/h)   for r >= R_max
 *              (r/R_max) * exp(-(r - R_max)^2 / (2 * sigma_r^2)) * I_0 * exp(-r/h)   for r < R_max
 *
 *              Note that I_0 = I_max / exp(-R_max/h), so that I_max is the intensity
 *              at R = R_max.
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
 *     [v0.1]  14 Sep 2021: Created as modification of func_logspiral2.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
// Use cmath instead of math.h to avoid GCC-5 problems with C++-11 and isnan()
//#include <math.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_logspiral_exp.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 11;
const char  PARAM_LABELS[][20] = {"PA", "ell", "m", "i_pitch", "R_i", "sigma_az", "gamma",
								"I_max", "h", "R_max", "sigma_trunc"};
const char  FUNCTION_NAME[] = "Logarithmic Spiral function (outer exp + inner-Gaussian truncation)";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;
const double  MIN_RADIUS = 0.001;

const char LogSpiralExp::className[] = "LogSpiralExp";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

LogSpiralExp::LogSpiralExp( )
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

void LogSpiralExp::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];        // line of nodes for projected circle
  ell = params[1 + offsetIndex];       // ellipticity of projected circle
  m = params[2 + offsetIndex];         // multiplicity of spiral
  i_pitch = params[3 + offsetIndex];   // pitch angle [degrees]
  R_i = params[4 + offsetIndex];       // radius where spiral crosses x=0 [pix; ring radius for infinite winding]
  sigma_az = params[5 + offsetIndex];  // Gaussian azimuthal width of spiral [pix]
  gamma = params[6 + offsetIndex];     // phase angle (azimuthal offset) for spiral pattern
  I_max = params[7 + offsetIndex];     // intensity at peak of spiral amplitude
  h = params[8 + offsetIndex];         // exponential radial scale length [pix]
  R_max = params[9 + offsetIndex];     // radius where radial intensity profile = I_max
  sigma_trunc = params[10 + offsetIndex];  // inner Gaussian radial sigma (for r < R_max)

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  // this is relative to internal orientation of spiral; corrects for
  // value of m (if user wants to see pattern rotated by 45 deg, we need
  // gamma = m*45 in the actual logSpiralFn equation)
  gamma_true = m * gamma * DEG2RAD;

  m_over_tani = m / tan(i_pitch * DEG2RAD);
  sigma_az_squared = sigma_az*sigma_az;
  twosigma_trunc_squared = 2.0*sigma_trunc*sigma_trunc;
  I_0 = I_max / exp(-R_max / h);   // nominal intensity of exponential at r=0
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for logarithmic spiral, evaluated
// at radius r and azimuthal angle phi; r is assumed to be positive, and phi
// is assumed to be in radians.
double LogSpiralExp::CalculateIntensity( double r, double phi )
{
  double  I_rad, truncationScaling;
  double  logSpiralFn, phi_term, spiralAzimuthalTerm;
  
  if (r <= 0.0) {
    r = MIN_RADIUS;  // to avoid problems with log(0)
  }
  logSpiralFn = m_over_tani * log(r/R_i) + gamma_true;
  phi_term = 1.0 - cos(m*phi - logSpiralFn);
  spiralAzimuthalTerm = exp( (-r*r/sigma_az_squared) * phi_term );

  // radial amplitude function
//   if (r < R_max) {  // inside the "truncation radius"
//     double  r_diff = R_max - r;
//     double  gaussianTerm = exp( -(r_diff*r_diff)/twosigma_trunc_squared );
//     truncationScaling = (r/R_max) * gaussianTerm;
//   }
//   else  // outside the ring
//     truncationScaling = 1.0;
//   I_rad = truncationScaling * I_0 * exp(-r/h);
  
  if (r < R_max) {   // inside R_max
    double  r_diff = R_max - r;
    I_rad = I_max * exp( -(r_diff*r_diff)/twosigma_trunc_squared );
  }
  else   // outside R_max
    I_rad = I_0 * exp(-r/h);

  return I_rad * spiralAzimuthalTerm;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// r is *non-negative*!
double LogSpiralExp::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled;
  double  r, phi, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  
  // NOTE: use atan2, *not* atan(y/x) [latter causes m=odd to be messed up]
  phi = atan2(yp_scaled, xp);

  nSubsamples = CalculateSubsamples(r);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-pixel
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
        r = sqrt(x_diff*x_diff + y_diff*y_diff);
        phi = atan(y_diff/x_diff);
        theSum += CalculateIntensity(r, phi);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r, phi);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a (scaled) distance of r away from the center of the
// ring.
// For now, we don't do any subsampling.
int LogSpiralExp::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  return nSamples;
}



/* END OF FILE: func_logspiral_exp.cpp --------------------------------- */

