/*   Class interface definition for func_sersic.cpp
 *   VERSION 0.2
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
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
  public:
    // Constructors:
    Sersic1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  n, mu_e, r_e;   // parameters
    double  I_e, n2, bn, invn;
};
