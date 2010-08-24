/*   Class interface definition for func_gen-exp.cpp
 *   VERSION 0.4
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * exponential; shape is generalized ellipse.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * c0 = params[2 + offsetIndex];    -- ellipse shape parameter (<0 for disky, >0 for boxy)
 * I_0 = params[3 + offsetIndex ]; -- central intensity (ADU/pix)
 * h = params[4 + offsetIndex ];   -- exp. scale length (pixels)
 *
 *
 */


// CLASS GenExponential:

#include "function_object.h"

#define CLASS_SHORT_NAME  "Exponential_GenEllipse"


class GenExponential : public FunctionObject
{
  public:
    // Constructors:
    GenExponential( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = CLASS_SHORT_NAME; };

  protected:
    double CalculateIntensity( double r );
    double CalculateRadius( double deltaX, double deltaY );  // new!
    int  CalculateSubsamples( double r );

  private:
    double  x0, y0, PA, ell, c0, I_0, h;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities (basic geometry)
    double  ellExp, invEllExp;         // more useful quantities
};
