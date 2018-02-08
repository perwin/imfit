/*   Class interface definition for func_logspiral_gauss.cpp
 *   VERSION 0.1
 *
 *   Class (derived from FunctionObject; function_object.h) which produces a simple/crude
 * logarithmic spiral with overall amplitude modulated by radial Gaussian function
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   // line of nodes for projected circle
 * ell = params[1 + offsetIndex];  // ellipticity of projected circle
 * m = params[2 + offsetIndex];   // multiplicity of spiral (m=2 for standard 2-arm spiral)
 * i_pitch = params[3 + offsetIndex ];   // pitch or winding angle [degrees]
 * I_0 = params[4 + offsetIndex];   // intensity at peak of spiral/ring
 * R_i = params[5 + offsetIndex];   // radius where spiral crosses x=0 [ring for infinite winding]
 * sigma = params[6 + offsetIndex];  // Gaussian width of spiral
 * gamma = params[7 + offsetIndex];  // phase angle (azimuthal offset) for spiral pattern
 * R_max = params[8 + offsetIndex];  // radius where spiral pattern is maximum
 * sigma_max = params[9 + offsetIndex];  // sigma for R_max
 *
 */


// CLASS LogSpiralGauss:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "LogSpiralGauss"


class LogSpiralGauss : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    LogSpiralGauss( );
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
    double  x0, y0, PA, ell, m, i_pitch, I_0, R_i, sigma, gamma, R_max, sigma_max;   // parameters
    double  q, PA_rad, cosPA, sinPA;
    double  m_over_tani;   // other useful quantities
    double  sigma_squared, twosigma_max_squared;
};
