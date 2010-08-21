/*   Class interface definition for func1d_sech2.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius (or height) for a 1-D
 * sech^2 profile.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ]; -- central surf. brightness (mag/arcsec^2)
 * h = params[1 + offsetIndex ];    -- scale length
 *
 *
 */


// CLASS Sech21D:

#include "function_object.h"


class Sech21D : public FunctionObject
{
  public:
    // Constructors:
    Sech21D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, mu_0, h;   // parameters
    double  I_0;   // other useful quantities
};
