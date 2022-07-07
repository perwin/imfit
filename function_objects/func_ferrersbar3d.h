/*   Class interface definition for func_ferrersbar3d.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the integrated intensity of a 3D Ferrers ellipsoid ("bar")
 * seen at specified inclination.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component line of nodes, rel. to image +x axis
 * inclination = params[1 + offsetIndex];  -- inclination to line of sight (i=0 for face-on)
 * barPA = params[2 + offsetIndex ];  -- position of bar major-axis relative to line of nodes,
 *                                       in component equatorial plane
 * J_0 = params[3 + offsetIndex ];  -- central luminosity density (ADU)
 * R = params[4 + offsetIndex ];   -- radial length of long axis of bar (semi-major axis a)
 * q = params[5 + offsetIndex ];   -- b/a (ratio of minor planar axis to major planar axis)
 * q_z = params[6 + offsetIndex ];   -- c/a (ratio of z-axis to major planar axis)
 * n = params[7 + offsetIndex ];   -- exponent of Ferrers function
 *
 *
 */


// CLASS FerrersBar3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;


/// \brief Class for image function using LOS integration through 3D Ferrers bar model
class FerrersBar3D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    FerrersBar3D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool ParameterUnitsExist( );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, PA, inclination, barPA, J_0, R_bar, q, q_z, n;   // parameters
    double  PA_rad, cosPA, sinPA, barPA_rad, cosBarPA, sinBarPA;   // other useful quantities
    double  inc_rad, cosInc, sinInc, a2, b2, c2, integrationLimit;
    gsl_function  F;
};

