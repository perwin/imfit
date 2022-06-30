/* FILE: func_logspiral_arc.cpp ---------------------------------------- */
/* VERSION 0.1
 *
 *   Class derived from FunctionObject; function_object.h) which produces a simple/crude
 * 1-arm logarithmic spiral where intensity is at maximum along arm at user-specified
 * angle theta_max and decays at larger and smaller angles following 2-sided Gaussian.
 *
 *   Logarithmic spiral equation:
 *   R_arm = I_max * exp(tan_i_pitch * theta)
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
 *     [v0.1]  18 Nov 2021: Created as modification of func_gaussian-ring.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <tuple>

#include "func_logspiral_arc.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 9;
const char  PARAM_LABELS[][20] = {"PA", "ell", "r_scale", "i_pitch", "theta_max", 
								"I_max", "sigma_r", "sigma_theta_ccw", "sigma_theta_cw"};
const char  FUNCTION_NAME[] = "Logarithmic Spiral function (2-sided Gaussian radial modulation)";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 9;
const double  MIN_RADIUS = 0.001;

const char LogSpiralArc::className[] = "LogSpiralArc";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

LogSpiralArc::LogSpiralArc( )
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

void LogSpiralArc::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];   // line of nodes for projected circle [deg]
  ell = params[1 + offsetIndex];  // ellipticity of projected circle
  r_scale = params[2 + offsetIndex ];   // radial size of spiral arc [pix]
  i_pitch = params[3 + offsetIndex ];   // pitch or winding angle [degrees]
  theta_max = params[4 + offsetIndex];  // angle of maximum intensity [deg]
  I_max = params[5 + offsetIndex];   // maximum intensity of spiral arc
  sigma_r = params[6 + offsetIndex];  // Gaussian width of spiral in radial direction [pix]
  sigma_theta_ccw = params[7 + offsetIndex];  // sigma (in CCW direction) for azimuthal Gaussian [deg]
  sigma_theta_cw = params[8 + offsetIndex];  // sigma (in CW direction) for azimuthal Gaussian [deg]

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  tan_i_pitch = tan(i_pitch * DEG2RAD);

  twosigma_r_squared = 2.0 * sigma_r*sigma_r;
  theta_max_rad = theta_max * DEG2RAD;
  // theta_sp_ref = reference angle for computing log.spiral R_arm
  // (= 180 deg away from discontinuity vector)
  // for theta_max > 0, < 180:
  //	use theta_final = -theta_max
  // for theta_max > 180, < 360   225(135) -> 45, 270(90) -> 90
  // 	use theta_final = 360 - theta_max
  if ((theta_max_rad >= 0.0) && (theta_max_rad <= M_PI))
    theta_sp_ref = -theta_max_rad;
  else
    theta_sp_ref = 2.0*M_PI - theta_max_rad;
  // do this *after* using theta_max_rad to determin theta_sp
  // temporarily comment-out so theta_max stays in std angle coords
  // theta_max_rad = theta_max_rad + M_PI_2;

  sigma_theta_ccw_rad = sigma_theta_ccw * DEG2RAD;
  sigma_theta_cw_rad = sigma_theta_cw * DEG2RAD;
  twosigma_theta_ccw_squared = 2.0 * sigma_theta_ccw_rad*sigma_theta_ccw_rad;
  twosigma_theta_cw_squared = 2.0 * sigma_theta_cw_rad*sigma_theta_cw_rad;

  // Finally, work out correct radial scaling, so that R_arm at theta_std = theta_max
  // = requested r_scale
  double  theta_rel_max = GetThetaRel(theta_max_rad);
  R0 = r_scale / exp(tan_i_pitch * theta_rel_max);

}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a 1-armed logarithmic spiral, evaluated
// at radius r and azimuthal angle theta_rel; r is assumed to be positive
double LogSpiralArc::CalculateIntensity( double r, double theta_rel, double theta_std )
{
  double  R_arm, r_diff, r_squared, theta_diff, theta_diff_squared;
  double  I_theta;
  
  R_arm = R0 * exp(tan_i_pitch * theta_rel);
  r_diff = r - R_arm;
  r_squared = r_diff*r_diff;
  theta_diff = GetThetaDiff(theta_std);
  theta_diff_squared = theta_diff*theta_diff;
  if (theta_diff >= 0)   // offset from theta_max in CCW direction
    I_theta = I_max * exp(-theta_diff_squared/twosigma_theta_ccw_squared);
  else
    I_theta = I_max * exp(-theta_diff_squared/twosigma_theta_cw_squared);
  return I_theta * exp(-r_squared/twosigma_r_squared);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// r is *non-negative*!
double LogSpiralArc::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled;
  double  r, theta_std, theta_rel;
  double  totalIntensity;
  double  yx_ratio;
  int  nSubsamples;
    
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  
  std::tie(theta_std, theta_rel) = GetThetas(x_diff, y_diff);
  
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
        xp = x_diff*cosPA + y_diff*sinPA;
        yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
        r = sqrt(xp*xp + yp_scaled*yp_scaled);
        std::tie(theta_std, theta_rel) = GetThetas(x_diff, y_diff);
        theSum += CalculateIntensity(r, theta_rel, theta_std);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r, theta_rel, theta_std);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a (scaled) distance of r away from the center of the
// ring.
// For now, we don't do any subsampling.
int LogSpiralArc::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  return nSamples;
}


/* ---------------- PROTECTED METHOD: GetStdAngle --------------------------------- */
// Given Cartesian position x_diff,y_diff, returns angle of vector from origin
// to that position, in standard form (measured CCW from +x axis, ranges from
// 0 to 2pi)
double LogSpiralArc::GetStdAngle( double x_diff, double y_diff )
{
  double  theta_atan2 = atan2(y_diff, x_diff);
  // atan2 returns angles in range (0, pi/2, pi, -pi/2, 0), but we want (0, pi, 2pi)
  return fmod(2.0*M_PI + theta_atan2, 2.0*M_PI);
}


/* ---------------- PROTECTED METHOD: GetThetaRel --------------------------------- */
double LogSpiralArc::GetThetaRel( double theta_std )
{
  double  theta_rel = theta_std - (M_PI - theta_sp_ref);
  if (theta_rel < 0.0)
    theta_rel += 2*M_PI;
  return theta_rel;
}


/* ---------------- PROTECTED METHOD: GetThetas ----------------------------------- */
std::tuple<double, double> LogSpiralArc::GetThetas( double x_diff, double y_diff )
{
  double  theta_std = GetStdAngle(x_diff, y_diff);
  double  theta_rel = GetThetaRel(theta_std);
  return std::make_tuple(theta_std, theta_rel);
}


/* ---------------- PROTECTED METHOD: GetThetaDiff -------------------------------- */
// Given angle in image of vector from component center to current pixel coordinate
// in form of standard angle coordinate theta_std, computes relative angular distance
// to theta_max_rad, such that distances in CCW direction from theta_max are positive,
// and distances in CW direction are negative.
double LogSpiralArc::GetThetaDiff( double theta_std )
{
  double theta, theta_diff;
  
  theta = theta_std;
  if (theta_max_rad > theta)
    theta += 2.0*M_PI;
  theta_diff = theta - theta_max_rad;
  if (theta_diff > M_PI)
    theta_diff = -(2.0*M_PI - theta_diff);
  return theta_diff;
}


/* END OF FILE: func_logspiral_arc.cpp --------------------------------- */
