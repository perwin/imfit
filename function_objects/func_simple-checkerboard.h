/*   Class interface definition for func_simple-checkerboard.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces a checkerboard grid.
 *
 */


// CLASS SimpleCheckerboard:

#include "function_object.h"



/// Class for image function with constant background level
class SimpleCheckerboard : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    SimpleCheckerboard( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool  IsBackground( );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, I_pos;
    int  step;
};
