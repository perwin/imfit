/* FILE: func_logspiral_gauss.cpp -------------------------------------- */
/* VERSION 0.1
 *
 *   Class derived from FunctionObject; function_object.h) which produces a simple/crude
 * logarithmic spiral.
 *   Right now, we only handle face-on spirals.
 *
 *   Logarithmic spiral equation:
 *   For radius r, theta(peak) = 1/b * ln(r/a)
 *      b = winding/pitch angle
 *      As b -> 0, 1/b -> infinity: spiral -> circle of radius a
 *      As b -> infinity, 1/b -> 0: spiral -> straight line
 *
 *   Eqn. 8 of Junqueira+2013:
 *   I(r,phi) = I_rad(r) * exp( -R^2/sigma^2 * (1 - cos(m*phi - f_m(R))) )
 *
 *      where:
 *         I_rad(r) describes radial decay of spiral amplitude
 *            [currently, I_rad(r) = Gaussian]
 *         (1 - cos(...)) is a term which varies sinusoidally from 0 to 2, with
 *            m peaks (where value = 2)
 *         f_m(R) is "shape function", which describes the logarithmic spiral;
 *            it acts as a phase-angle offset, whose angular position shifts with
 *            radius following.
 *            f_m(R) = (m/tan i) * ln(R/R_i) [+ gamma]
 *               where i = pitch angle, R_i = "point where the spiral crosses the
 *               coordinate x = 0" (and radius of rings if spiral were infinitely wound)
 *               and gamma is azimuthal-angle offset for spiral pattern
 *
 *               for R_i = 100 and m=2, spiral will cross +x axis at x = 100, 51, 26, 13
 *
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
 *     [v0.1]  21 March 2014: Created as modification of func_gaussian-ring.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_logspiral_gauss.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 11;
const char  PARAM_LABELS[][20] = {"PA", "ell", "m", "i_wind", "I_0", "R_i", "sigma", "gamma",
								"R_max", "sigma_max_in", "sigma_max_out"};
const char  FUNCTION_NAME[] = "Logarithmic Spiral function (2-sided Gaussian radial modulation)";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;
const double  MIN_RADIUS = 0.001;

const char LogSpiralGauss::className[] = "LogSpiralGauss";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

LogSpiralGauss::LogSpiralGauss( )
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

void LogSpiralGauss::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];   // line of nodes for projected circle
  ell = params[1 + offsetIndex];  // ellipticity of projected circle
  m = params[2 + offsetIndex];   // multiplicity of spiral (m=2 for standard 2-arm spiral)
  i_pitch = params[3 + offsetIndex];   // pitch angle [degrees]
  I_0 = params[4 + offsetIndex];   // intensity at peak of spiral/ring
  R_i = params[5 + offsetIndex];   // radius where spiral crosses x=0 [ring for infinite winding]
  sigma = params[6 + offsetIndex];  // Gaussian width of spiral
  gamma = params[7 + offsetIndex];  // phase angle (azimuthal offset) for spiral pattern
  R_max = params[8 + offsetIndex];  // radius where intensity is maximum
  sigma_max_in = params[9 + offsetIndex];  // inner sigma for R_max
  sigma_max_out = params[10 + offsetIndex];  // outer sigma for R_max

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  m_over_tani = m / tan(i_pitch * DEG2RAD);
  sigma_squared = sigma*sigma;
  twosigma_max_in_squared = 2.0*sigma_max_in*sigma_max_in;
  twosigma_max_out_squared = 2.0*sigma_max_out*sigma_max_out;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for logarithmic spiral, evaluated
// at radius r and azimuthal angle phi; r is assumed to be positive
double LogSpiralGauss::CalculateIntensity( double r, double phi )
{
  double  I, logSpiralFn, r_scaled;
  double  r_diff = fabs(r - R_max);
  double  twosigma_r_squared;
  
  if (r <= 0.0) {
    r = MIN_RADIUS;
  }
  logSpiralFn = m_over_tani * log(r/R_i) + gamma;
  double  phi_term = 1.0 - cos(m*phi - logSpiralFn);
  double  exp_stuff = exp( (-r*r/sigma_squared) * phi_term );

  // decide which sigma to use for each Gaussian component
  if (r < R_max) {  // inside the ring
    twosigma_r_squared = twosigma_max_in_squared;
  }
  else {  // outside the ring
    twosigma_r_squared = twosigma_max_out_squared;
  }
  r_scaled = r - R_max;
  I = I_0 * exp(-(r_diff*r_diff)/twosigma_r_squared);
  return I * exp_stuff;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// r is *non-negative*!
double LogSpiralGauss::GetValue( double x, double y )
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
  phi = atan(y_diff/x_diff);
  
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
int LogSpiralGauss::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  return nSamples;
}



/* END OF FILE: func_logspiral_gauss.cpp ------------------------------- */
