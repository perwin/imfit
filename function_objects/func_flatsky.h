/*   Class interface definition for func_flatsky.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces constant intensity per pixel.
 *
 */


// CLASS FlatSky:

#include "function_object.h"


class FlatSky : public FunctionObject
{
  public:
    // Constructors:
    FlatSky( bool subsampling );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now


  private:
    double  x0, y0, I_sky;
};
