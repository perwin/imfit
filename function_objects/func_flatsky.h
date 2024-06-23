/*   Class interface definition for func_flatsky.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces constant intensity per pixel.
 *
 */


#ifndef _FUNC_FLATSKY_H_
#define _FUNC_FLATSKY_H_

// CLASS FlatSky:

#include "function_object.h"



/// Class for image function with constant background level
class FlatSky : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    FlatSky( );
    // redefined method/member function:
    void AdjustParametersForImage( const double inputFunctionParams[], 
									double adjustedFunctionParams[], int offsetIndex );
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, I_sky;
};


#endif /* _FUNC_FLATSKY_H_ */
