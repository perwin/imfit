/*   Class interface definition for func_gaussian.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a circular
 * Gaussian.  This will generate a normalized profile (total flux = 1.0).
 *
 */


// CLASS MoffatPSF:

#include "function_object.h"


class Gaussian : public FunctionObject
{
  public:
    // Constructors:
    Gaussian( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now


  private:
    double  x0, y0, A, sigma, twosigma_squared;
};
