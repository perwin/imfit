/*   Class interface definition for func_gen-exp2.cpp
 *
 *   Experimental version of GeneralizedExponential with linear
 * interpolation of c0 parameter.
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


// CLASS GenExponential2:

#include "function_object.h"



/// \brief Class for image function with generalized-elliptical isophotes and exponential profile
class GenExponential2 : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    GenExponential2( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };

  protected:
    double CalculateIntensity( double r );
    double CalculateRadius( double deltaX, double deltaY );  // new!
    int  CalculateSubsamples( double r );

  private:
    double  x0, y0, PA, ell, I_0, h;   // parameters
    double  c0_1, c0_2, r1, r2;            // extra parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities (basic geometry)
};
