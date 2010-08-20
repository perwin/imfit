/*   Class interface definition for func1d_delta.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces a 1-D delta function.  (The location of the
 * delta function maximum is described by x0.)
 *
 * PARAMETERS:
 * I = params[0 + offsetIndex ];   -- intensity (counts)
 *
 *
 */


// CLASS Delta1D:

#include "function_object.h"


class Delta1D : public FunctionObject
{
  public:
    // Constructors:
    Delta1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now


  private:
    double  x0, I;   // parameters
};
