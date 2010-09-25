/*   Class interface definition for func1d_vdksech2.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius (or height) for a 1-D
 * profile following van der Kruit's (1988, A&A) modified sech family of
 * functions.
 *   alpha = exponent for sech function
 *   as alpha --> 0, we approach an exponential profile
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ];  -- central surf. brightness (mag/arcsec^2)
 * z0 = params[1 + offsetIndex ];    -- scale height
 * alpha = params[2 + offsetIndex];  -- exponent (= 2/n in van der Kruit 1988)
 *
 *
 */


// CLASS vdKSech1D:

#include "function_object.h"


class vdKSech1D : public FunctionObject
{
  public:
    // Constructors:
    vdKSech1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, mu_0, z0, alpha, scaledZ0;   // parameters
    double  I_0;   // other useful quantities
};
