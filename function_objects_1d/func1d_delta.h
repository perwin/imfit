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
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Delta1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, I;   // parameters
};
