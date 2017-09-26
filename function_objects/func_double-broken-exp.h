/*   Class interface definition for func_double-broken-exp.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * component with a double broken-exponential profile.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * I_0 = params[2 + offsetIndex ]; -- central intensity of inner exponential (ADU/pix)
 * h1 = params[3 + offsetIndex ];   -- inner exp. scale length (pixels)
 * h2 = params[4 + offsetIndex ];   -- middle exp. scale length (pixels)
 * h3 = params[5 + offsetIndex ];   -- outer exp. scale length (pixels)
 * r_break1 = params[6 + offsetIndex ];   -- inner break radius (pixels)
 * r_break2 = params[7 + offsetIndex ];   -- outer break radius (pixels)
 * alpha1 = params[8 + offsetIndex ];     -- smoothness/sharpness of break
 * alpha2 = params[9 + offsetIndex ];     -- smoothness/sharpness of break
 *
 *
 */


// CLASS DoubleBrokenExponential:

#include "function_object.h"


/// Class for image function with elliptical isophotes and broken-exponential profile
class DoubleBrokenExponential : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    DoubleBrokenExponential( );
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
    double  x0, y0, PA, ell, I_0, h1, h2, h3, r_b1, r_b2, alpha1, alpha2;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities
    double  exponent2, exponent3, I_0_times_S, delta_Rb1_scaled, delta_Rb2_scaled;   // other useful quantities
};
