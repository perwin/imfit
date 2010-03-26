/*   Class interface definition for func_exp.cpp
 *   VERSION 0.4
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * exponential.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * I_0 = params[2 + offsetIndex ]; -- central intensity (ADU/pix)
 * h = params[3 + offsetIndex ];   -- exp. scale length (pixels)
 *
 *
 */


// CLASS Exponential:

#include "function_object.h"


class Exponential : public FunctionObject
{
  public:
    // Constructors:
    Exponential( bool subsampling );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now


  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, I_0, h;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities
};
