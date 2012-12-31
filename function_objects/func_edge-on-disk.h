/*   Class interface definition for func_edge-on-disk.cpp
 *   VERSION 0.3
 *
 *   Class derived from FunctionObject (function_object.h/cpp) which produces
 * surface brightnesses for a generalized edge-on exponential disk, using 
 * Bessel-function solution of van der Kruit & Searle (1981) for radial profile 
 * and generalized sech function (van der Kruit 1988) for vertical profile:
 *
 *      Sigma(r,z) = Sigma(0,) * (r/h) * K_1(r/h) * sech^(2/n)(n*z/(2*z0))
 *
 *    with Sigma(0,0) = 2 * h * L_0
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * L_0 = params[1 + offsetIndex ]; -- intensity scaling (ADU/pix)
 * h = params[2 + offsetIndex ];   -- radial exp. scale length (pixels)
 * n = params[3 + offsetIndex ];     -- exponent used in sech vertical function
 * z0 = params[4 + offsetIndex ] -- vertical scale height
 *
 *
 */


// CLASS EdgeOnDisk:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "EdgeOnDisk"


class EdgeOnDisk : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    EdgeOnDisk( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r, double z );
    int  CalculateSubsamples( double r, double z );


  private:
    double  x0, y0, PA, L_0, h, n, z_0;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  alpha, scaledZ0, Sigma_00, two_to_alpha;   // other useful quantities
};
