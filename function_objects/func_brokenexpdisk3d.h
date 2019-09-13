/*   Class interface definition for func_brokenexpdisk3d.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the integrated intensity of a 3D disk with a broken-exponential radial
 * profile and vertical exponential profile, seen at specified inclination.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];     -- PA of component line of nodes, rel. to image +x axis
 * inclination = params[1 + offsetIndex];  -- inclination to line of sight (i=0 for face-on)
 * J_0 = params[2 + offsetIndex ];   -- central luminosity density (ADU)
 * h1 = params[3 + offsetIndex ];    -- inner exp. scale length (pixels)
 * h2 = params[4 + offsetIndex ];    -- outer exp. scale length (pixels)
 * r_b = params[5 + offsetIndex ];   -- break radius (pixels)
 * alpha = params[6 + offsetIndex ]; -- smoothness/sharpness of break [1/pixels]
 * n = params[7 + offsetIndex ];     -- exponent used in sech vertical function
 * z_0 = params[8 + offsetIndex ];   -- vertical scale height
 *
 */


// CLASS BrokenExponentialDisk3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;



/// \brief Class for image function using LOS integration through 3D model with broken-exponential
///        radial and vertical exponential luminosity-density profiles
class BrokenExponentialDisk3D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    BrokenExponentialDisk3D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, PA, inclination, J_0, h1, h2, r_b, alpha, n, z_0;   // parameters
    double  PA_rad, cosPA, sinPA, inc_rad, cosInc, sinInc;   // other useful quantities
    double  exponent, J_0_times_S, delta_Rb_scaled;
    double  scaledZ0, two_to_alpha, alphaVert;
    gsl_function  F;
};

