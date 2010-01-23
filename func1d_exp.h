/*   Class interface definition for func1d_exp.cpp
 *   VERSION 0.2
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
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
    void  Setup( double params[], int offsetIndex );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  mu_0, h;   // parameters
    double  I_0;   // other useful quantities
};
