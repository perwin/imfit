/*   Class interface definition for func_gen-sersic.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * Sersic; shape is generalized ellipse.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA =  params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];   -- ellipticity
 * c0 =  params[2 + offsetIndex];   -- ellipse shape parameter (<0 for disky, >0 for boxy)
 * n =   params[3 + offsetIndex ];  -- Sersic index
 * I_e = params[4 + offsetIndex ];  -- half-light-radius intensity (ADU)
 * r_e = params[5 + offsetIndex ];  -- half-light radius (pixels)
 *
 *
 */


// CLASS GenSersic:

#include "function_object.h"


class GenSersic : public FunctionObject
{
  public:
    // Constructors:
    GenSersic( bool subsampling );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now


  protected:
    double CalculateIntensity( double r );
    double CalculateRadius( double deltaX, double deltaY );  // new!
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, c0, n, I_e, r_e;   // parameters
    double  n2, bn, invn;
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities (basic geometry)
    double  ellExp, invEllExp;         // more useful quantities
};
