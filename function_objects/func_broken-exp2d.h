/*   Class interface definition for func_broken-exp.cpp
 *   VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces an approximate 2D edge-on broken-exponential (with radial
 * broken-exponential profile and vertical exponential profile).
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * I_0 = params[1 + offsetIndex ]; -- central intensity of inner exponential (ADU/pix)
 * h1 = params[2 + offsetIndex ];   -- inner exp. scale length (pixels)
 * h2 = params[3 + offsetIndex ];   -- outer exp. scale length (pixels)
 * r_break = params[4 + offsetIndex ];   -- break radius (pixels)
 * alpha = params[5 + offsetIndex ];     -- smoothness/sharpness of break
 * h_z = params[6 + offsetIndex ] -- vertical scale height
 *
 *
 */


// CLASS BrokenExponential2D:

#include "function_object.h"

#define CLASS_SHORT_NAME  "BrokenExponential2D"


class BrokenExponential2D : public FunctionObject
{
  public:
    // Constructors:
    BrokenExponential2D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = CLASS_SHORT_NAME; };


  protected:
    double CalculateIntensity( double r, double z );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, I_0, h1, h2, r_b, alpha, h_z;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  exponent, I_0_times_S, delta_Rb_scaled;   // other useful quantities
};
