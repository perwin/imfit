/*   Class interface definition for func1d_exp.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D
 * exponential.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ]; -- central surf. brightness (mag/arcsec^2)
 * h = params[1 + offsetIndex ];    -- exp. scale length
 *
 *
 */


// CLASS Exponential1D:

#include "function_object.h"


class Exponential1D : public FunctionObject
{
  public:
    // Constructors:
    Exponential1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, mu_0, h;   // parameters
    double  I_0;   // other useful quantities
};
