/*   Class interface definition for func_sersic.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D
 * Sersic.
 *
 * PARAMETERS:
 * n =    params[0 + offsetIndex ];  -- Sersic index
 * mu_e = params[1 + offsetIndex ];  -- half-light-radius surf. brightness (mag/arcsec^2)
 * r_e =  params[2 + offsetIndex ];  -- half-light radius
 *
 *
 */


// CLASS Sersic1D:

#include "function_object.h"


class Sersic1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Sersic1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, n, mu_e, r_e;   // parameters
    double  I_e, n2, bn, invn;
};
