/*   Class interface definition for func_edge-on-disk.cpp
 *   VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces a generalized edge-on exponential disk, using Bessel-function
 * solution (van der Kruit & Searle 1981) for radial profile and generalized
 * sech function (van der Kruit 1988) for vertical profile.  My version of this
 * looks like the following (Sigma = surface brightness in counts/pixel):
 *
 *      Sigma(r,z) = I_0 * (r/h) * K_1(r/h) * sech^alpha(r/(alpha*z0))
 *
 *    with Sigma(0,0) = 2 * h * I_0
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * I_0 = params[1 + offsetIndex ]; -- intensity scaling (ADU/pix)
 * h = params[2 + offsetIndex ];   -- radial exp. scale length (pixels)
 * alpha = params[3 + offsetIndex ];     -- exponent for sech vertical function
 * z0 = params[4 + offsetIndex ] -- vertical scale height
 *
 *
 */


// CLASS EdgeOnDisk:

#include "function_object.h"

#define CLASS_SHORT_NAME  "EdgeOnDisk"


class EdgeOnDisk : public FunctionObject
{
  public:
    // Constructors:
    EdgeOnDisk( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = CLASS_SHORT_NAME; };


  protected:
    double CalculateIntensity( double r, double z );
    int  CalculateSubsamples( double r, double z );


  private:
    double  x0, y0, PA, I_0, h, alpha, z_0;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  scaledZ0, Sigma_00;   // other useful quantities
};
