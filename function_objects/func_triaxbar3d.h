/*   Class interface definition for func_triaxbar3d.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the integrated intensity of a 3D triaxial ellipsoid ("triaxbar")
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
 * h_0 = params[4 + offsetIndex ];  -- radial scale length along major axis of bar
 * q = params[5 + offsetIndex ];   -- b/a (ratio of minor planar axis to major planar axis)
 * q_z = params[6 + offsetIndex ];   -- c/a (ratio of z-axis to major planar axis)
 *
 *
 */


// CLASS TriaxBar3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;


/// \brief Class for image function using LOS integration through 3D Ferrers bar model
class TriaxBar3D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    TriaxBar3D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, PA, inclination, barPA, J_0, sigma, q, q_z;   // parameters
    double  PA_rad, cosPA, sinPA, barPA_rad, cosBarPA, sinBarPA;   // other useful quantities
    double  inc_rad, cosInc, sinInc, twosigma_squared;
    gsl_function  F;
};

