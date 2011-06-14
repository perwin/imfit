/*   Class interface definition for func_flatsky.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces constant intensity per pixel.
 *
 */


// CLASS FlatSky:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "FlatSky"


class FlatSky : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    FlatSky( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, I_sky;
};
