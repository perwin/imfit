/*   Class interface definition for func_gaussianring3d.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h), which produces projected
 * surface intensity for an elliptical 3D ring (luminosity density = Gaussian centered at 
 * r0 with width sigma and vertical exponential with scale heigh h_z), with major-axis
 * position angle PA_ring relative to the line of nodes; the overall orientation is
 * defined by the PA of the line of nodes (PA) and the inclination (inc).
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component line-of-nodes, rel. to +x axis
 * inclination = params[1 + offsetIndex];  -- inclination to line of sight (i=0 for face-on)
 * ringPA = params[2 + offsetIndex];  -- PA of ring major axis, relative to line of nodes
 * ell = params[3 + offsetIndex];  -- ellipticity of ring
 * J_0 = params[4 + offsetIndex ];  -- central luminosity density (ADU)
 * a_ring = params[5 + offsetIndex ];   -- in-plane semi-major axis of circular ring
 * sigma = params[6 + offsetIndex ];   -- width of ring (Gaussian sigma)
 * h_z = params[7 + offsetIndex ];  -- vertical exp. scale height (pixels)
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
    double  x0, y0, PA, inclination, ringPA, ell, J_0, a_ring, sigma, h_z;   // parameters
    double  cosPA, sinPA, cosInc, sinInc;   // other useful quantities
    double  cosRingPA, sinRingPA, q, twosigma_squared;
    gsl_function  F;
};

