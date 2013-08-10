/*   Class interface definition for func_brokenexpdisk3d.cpp
 *   VERSION 0.2
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the integrated intensity of a 3D disk with a broken-exponential radial
 * profile and vertical exponential profile, seen at specified inclination.
 *
 *
 */


// CLASS BrokenExponentialDisk3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;

//#define CLASS_SHORT_NAME  "BrokenExponentialDisk3D"


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

