/*   Class interface definition for func_core-sersic.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D
 * Core-Sersic profile (Graham+2003, Trujillo+2004).
 *
 * PARAMETERS:
 * n =    params[0 + offsetIndex ];  -- Sersic index
 * mu_b = params[1 + offsetIndex ];  -- break-radius surf. brightness (mag/arcsec^2)
 * r_e =  params[2 + offsetIndex ];  -- half-light radius
 * r_b =  params[3 + offsetIndex ];  -- break radius
 * alpha =  params[4 + offsetIndex ];  -- sharpness of break
 * gamma =  params[5 + offsetIndex ];  -- inner power-law slope
 *
 *
 */


// CLASS CoreSersic1D:

#include "function_object.h"


class CoreSersic1D : public FunctionObject
{
  public:
    // Constructors:
    CoreSersic1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, n, mu_b, r_e, r_b, alpha, gamma;   // parameters
    double  I_b, n2, bn, invn, Iprime;
};
