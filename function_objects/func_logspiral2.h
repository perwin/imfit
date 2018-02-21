/*   Class interface definition for func_logspiral2.cpp
 *   VERSION 0.1
 *
 *   Class (derived from FunctionObject; function_object.h) which produces a simple
 * logarithmic spiral with basic surface-brightness specified by Eqn. 8 of Junqueira+2013,
 * with the radial amplitude is specified by an exponential for R > R_max and the
 * product of R/R_max, a Gaussian, and the exponential for R < R_max.
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
 * h = params[8 + offsetIndex];         // exponential radial scale length
 * R_max = params[9 + offsetIndex];     // inner truncation radius
 * sigma_trunc = params[10 + offsetIndex];  // inner Gaussian radial sigma (for r < R_max)
 *
 */


// CLASS LogSpiral2:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "LogSpiral2"


class LogSpiral2 : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    LogSpiral2( );
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
    double  x0, y0, PA, ell, m, i_pitch, R_i, sigma_az, gamma, I_0, h, R_max, sigma_trunc;   // parameters
    double  q, PA_rad, cosPA, sinPA;
    double  m_over_tani, gamma_rad;   // other useful quantities
    double  sigma_az_squared, twosigma_trunc_squared;
};
