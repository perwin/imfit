/*   Class interface definition for func_expdisk3d.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the integrated intensity of a 3D perfect exponential disk
 * seen at specified inclination.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * inclination = params[1 + offsetIndex];  -- inclination to line of sight (i=0 for face-on)
 * J_0 = params[2 + offsetIndex ];  -- central luminosity density (ADU)
 * h = params[3 + offsetIndex ];   -- radial exp. scale length (pixels)
 * h_z = params[4 + offsetIndex ];   -- vertical exp. scale height (pixels)
 *
 *
 */


// CLASS ExponentialDisk3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;

//#define CLASS_SHORT_NAME  "ExponentialDisk3D"


class ExponentialDisk3D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    ExponentialDisk3D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, PA, inclination, J_0, h, n, z_0;   // parameters
    double  PA_rad, cosPA, sinPA, inc_rad, cosInc, sinInc;   // other useful quantities
    double  scaledZ0, two_to_alpha, alpha;
    gsl_function  F;
};

