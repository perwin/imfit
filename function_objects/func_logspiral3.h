/*   Class interface definition for func_logspiral2.cpp
 *   VERSION 0.1
 *
 *   Class (derived from FunctionObject; function_object.h) which produces a simple
 * logarithmic spiral with basic surface-brightness specified by Eqn. 8 of Junqueira+2013,
 * with the radial amplitude specified by a broken exponential for R > R_max and the
 * product of R/R_max, a Gaussian, and the inner part of the broken exponential for R < R_max.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];        // line of nodes for projected circle
 * ell = params[1 + offsetIndex];       // ellipticity of projected circle
 * m = params[2 + offsetIndex];         // multiplicity of spiral
 * i_pitch = params[3 + offsetIndex];   // pitch angle [degrees]
 * R_i = params[4 + offsetIndex];       // radius where spiral crosses x=0 [ring for infinite winding]
 * sigma_az = params[5 + offsetIndex];  // Gaussian azimuthal width of spiral
 * gamma = params[6 + offsetIndex];     // phase angle (azimuthal offset) for spiral pattern
 * I_0 = params[7 + offsetIndex];       // intensity at peak of spiral amplitude
 * h1 = params[8 + offsetIndex];        // inner exponential radial scale length (inside r_b)
 * h2 = params[9 + offsetIndex];        // outer exponential radial scale length (outside r_b)
 * r_b = params[10 + offsetIndex];      // break radius for broken-exponential profile
 * alpha = params[11 + offsetIndex];    // sharpness parameter for broken-exponential profile
 * R_max = params[12 + offsetIndex];    // inner truncation radius
 * sigma_trunc = params[13 + offsetIndex];  // inner Gaussian radial sigma (for r < R_max)
 *
 */


// CLASS LogSpiral3:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "LogSpiral3"


class LogSpiral3 : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    LogSpiral3( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r, double phi );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, m, i_pitch, R_i, sigma_az, gamma;   // parameters
    double  I_0, h1, h2, r_b, alpha, R_max, sigma_trunc;   // parameters
    double  q, PA_rad, cosPA, sinPA, gamma_rad;
    double  m_over_tani;   // other useful quantities
    double  sigma_az_squared, twosigma_trunc_squared;
    double  exponent, I_0_times_S, delta_Rb_scaled;
};
