/*   Class interface definition for func_gaussianring3d.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h), which produces projected
 * surface intensity for a circular 3D ring (luminosity density = Gaussian centered at 
 * r0 with width sigma and vertical exponential with scale heigh h_z), seen at position 
 * angle PA and inclination inc.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * inclination = params[1 + offsetIndex];  -- inclination to line of sight (i=0 for face-on)
 * I_0 = params[2 + offsetIndex ];  -- central luminosity density (ADU)
 * R_0 = params[3 + offsetIndex ];   -- in-plane radius of circular ring
 * sigma = params[4 + offsetIndex ];   -- width of ring (Gaussian sigma)
 * h_z = params[5 + offsetIndex ];   -- vertical exp. scale height (pixels)
 *
 *
 */


// CLASS GaussianRing3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;

//#define CLASS_SHORT_NAME  "GaussianRing3D"


class GaussianRing3D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    GaussianRing3D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, PA, inclination, I_0, R_0, sigma, h_z;   // parameters
    double  PA_rad, cosPA, sinPA, inc_rad, cosInc, sinInc;   // other useful quantities
    gsl_function  F;
};

