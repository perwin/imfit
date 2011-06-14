/*   Class interface definition for func_gaussian.cpp
 *   VERSION 0.4
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a circular
 * Gaussian.  This will generate a normalized profile (total flux = 1.0).
 *
 */


// CLASS Gaussian:

#include "function_object.h"
#include <string>

//#define CLASS_SHORT_NAME  "Gaussian"


class Gaussian : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Gaussian( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now
    
    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, I_0, sigma;   // parameters
    double  twosigma_squared;
    double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};
