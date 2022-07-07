/*   Class interface definition for func_tilted-sky-plane.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces intensity per pixel as a tilted plane.
 *
 */


// CLASS TiltedSkyPlane:

#include "function_object.h"



/// Class for image function with constant background level
class TiltedSkyPlane : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    TiltedSkyPlane( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool  IsBackground( );
    bool ParameterUnitsExist( );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, I_0, m_x, m_y;
};
