/*   Class interface definition for func1d_gaussian2side.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D Gaussian.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ];   -- central surf. brightness (mag/arcsec^2)
 * sigma_left = params[1 + offsetIndex ];  -- sigma of the Gaussian's left side
 * sigma_right = params[2 + offsetIndex ];  -- sigma of the Gaussian's right side
 *
 *
 */


// CLASS Gaussian2Side1D:

#include "function_object.h"


class Gaussian2Side1D : public FunctionObject
{
  public:
    // Constructors:
    Gaussian2Side1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, mu_0, sigma_left, sigma_right;   // parameters
    double  I_0;
};
