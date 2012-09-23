/*   Class interface definition for func_core-sersic.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * Core-Sersic function.
 *
 * PARAMETERS:
 * x0 =  params[0 + offsetIndex];   -- center of component (pixels, x)
 * y0 =  params[1 + offsetIndex];   -- center of component (pixels, y)
 * PA =  params[2 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[3 + offsetIndex];   -- ellipticity
 * n =   params[4 + offsetIndex ];  -- Sersic index
 * I_b = params[5 + offsetIndex ];  -- break-radius surf. brightness (mag/arcsec^2)
 * r_e =  params[6 + offsetIndex ];  -- half-light radius
 * r_b =  params[7 + offsetIndex ];  -- break radius
 * alpha =  params[8 + offsetIndex ];  -- sharpness of break
 * gamma =  params[9 + offsetIndex ];  -- inner power-law slope
 *
 *
 */


// CLASS CoreSersic:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "CoreSersic"


class CoreSersic : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];

  public:
    // Constructors:
    CoreSersic( );
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
  double  x0, y0, PA, ell, n, I_b, r_e, r_b, alpha, gamma;   // parameters
  double  n2, bn, invn, Iprime;
  double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};
