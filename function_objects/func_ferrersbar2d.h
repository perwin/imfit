/*   Class interface definition for func_ferrersbar2d.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the surface brightness of a 2D Ferrers ellipsoid (as in
 * GALFIT -- really a 2D analog of the 3D Ferrers ellipsoid). Isophote
 * shapes are generalized ellipses.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA =  params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];   -- ellipticity
 * c0 =  params[2 + offsetIndex];   -- ellipse shape parameter (<0 for disky, >0 for boxy)
 * n =   params[3 + offsetIndex ];  -- profile shape parameter
 * I_0 = params[4 + offsetIndex ];  -- central intensity (ADU)
 * a_bar = params[5 + offsetIndex ];  -- semi-major axis (pixels)
 *
 */


// CLASS FerrersBar2D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;


/// \brief Class for image function using LOS integration through 3D Ferrers bar model
class FerrersBar2D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    FerrersBar2D( );
    // redefined method/member function:
    void AdjustParametersForImage( const double inputFunctionParams[], 
									double adjustedFunctionParams[], int offsetIndex );
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r );
    double CalculateRadius( double deltaX, double deltaY );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, c0, n, I_0, a_bar;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities (basic geometry)
    double  ellExp, invEllExp, a_bar2;         // more useful quantities
};

