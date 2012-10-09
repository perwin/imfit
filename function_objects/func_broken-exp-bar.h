/*   Class interface definition for func_broken-exp-bar.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * component with a broken-exponential profile.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * I_0 = params[2 + offsetIndex ]; -- central intensity of inner exponential (ADU/pix)
 * h1 = params[3 + offsetIndex ];   -- inner exp. scale length (pixels)
 * h2 = params[4 + offsetIndex ];   -- outer exp. scale length (pixels)
 * r_break = params[5 + offsetIndex ];   -- break radius (pixels)
 * alpha = params[6 + offsetIndex ];     -- smoothness/sharpness of break
 *
 *
 */


// CLASS BrokenExponentialBar:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "BrokenExponentialBar"


class BrokenExponentialBar : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    BrokenExponentialBar( );
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
    double CalculateIntensity( double x_in, double y_in );

    double  x0, y0, PA, I_0, h1, h2, r_b, alpha, sigma;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  exponent, I_0_times_S, delta_Rb_scaled, twosigma_squared;   // other useful quantities
};
