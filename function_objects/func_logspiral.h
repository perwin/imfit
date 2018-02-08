/*   Class interface definition for func_logspiral.cpp
 *   VERSION 0.1
 *
 *   Class (derived from FunctionObject; function_object.h) which produces a simple/crude
 * logarithmic spiral
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * ell = params[1 + offsetIndex ];   -- ellipticity of component
 * A = params[2 + offsetIndex ];   -- intensity scaling (ADU/pix)
 * R_ring = params[3 + offsetIndex ];
 * sigma = params[4 + offsetIndex ];   -- Gaussian sigma in radial direction
 *
 *
 */


// CLASS LogSpiral:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "LogSpiral"


class LogSpiral : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    LogSpiral( );
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
    double  x0, y0, PA, ell, m, i_pitch, I_0, R_i, sigma, gamma;   // parameters
    double  q, PA_rad, cosPA, sinPA;
    double  m_over_tani;   // other useful quantities
    double  sigma_squared;
};
