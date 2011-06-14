/*   Class interface definition for func_sersic.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * Sersic.
 *
 * PARAMETERS:
 * x0 =  params[0 + offsetIndex];   -- center of component (pixels, x)
 * y0 =  params[1 + offsetIndex];   -- center of component (pixels, y)
 * PA =  params[2 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[3 + offsetIndex];   -- ellipticity
 * n =   params[4 + offsetIndex ];  -- Sersic index
 * I_e = params[5 + offsetIndex ];  -- half-light-radius intensity (ADU)
 * r_e = params[6 + offsetIndex ];  -- half-light radius (pixels)
 *
 *
 */


// CLASS Sersic:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "Sersic"


class Sersic : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];

  public:
    // Constructors:
    Sersic( );
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
  double  x0, y0, PA, ell, n, I_e, r_e;   // parameters
  double  n2, bn, invn;
  double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};
