/*   Class interface definition for func_flat-exp.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * "flat/exponential" profile (constant intensity for r < r_break, exponential
 * decline for r > r_break).
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * I_0 = params[2 + offsetIndex ]; -- central intensity (ADU/pix)
 * h = params[3 + offsetIndex ];   -- exp. scale length (pixels)
 * r_break = params[4 + offsetIndex ];   -- break radius (pixels)
 *
 *
 */


// CLASS FlatExponential:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "FlatExponential"


class FlatExponential : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    FlatExponential( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, I_0, h, r_b, alpha;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities
    double  exponent, I_0_times_S;   // other useful quantities
};
