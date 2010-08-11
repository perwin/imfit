/*   Class interface definition for func1d_moffat.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D Moffat profile.
 *
 * PARAMETERS:
 * I_0 = params[0 + offsetIndex ];   -- central intensity (counts)
 * fwhm = params[1 + offsetIndex ];
 * beta = params[2 + offsetIndex ];
 *
 *
 */


// CLASS Moffat1D:

#include "function_object.h"


class Moffat1D : public FunctionObject
{
  public:
    // Constructors:
    Moffat1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, I_0, fwhm, beta;   // parameters
    double  alpha;
};
