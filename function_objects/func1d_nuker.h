/*   Class interface definition for func_core-sersic.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D
 * Nuker-Law profile (Lauer+1995, Byun+1996).
 *
 * PARAMETERS:
 * alpha =    params[0 + offsetIndex ];  -- sharpness of break
 * beta = params[1 + offsetIndex ];  -- inner power-law slope
 * gamma =  params[2 + offsetIndex ];  -- outer power-law slope
 * r_b =  params[3 + offsetIndex ];  -- break radius
 * mu_b =  params[4 + offsetIndex ];  -- break-radius surf. brightness (mag/arcsec^2)
 *
 *
 */


// CLASS NukerLaw1D:

#include "function_object.h"


class NukerLaw1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    NukerLaw1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, alpha, beta, gamma, r_b, mu_b;   // parameters
    double  Iprime, exponent;
};
