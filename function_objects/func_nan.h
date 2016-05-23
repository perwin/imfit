/*   Class interface definition for func_nan.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which always returns NaN values (for testing purposes).
 *
 */


// CLASS NaNFunc:

#include "function_object.h"



/// Class for image function generating NaN for all pixels (for testing)
class NaNFunc : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    NaNFunc( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, I_sky;
};
