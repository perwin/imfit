/*   Class interface definition for func_flatbar.cpp
 * 
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity for the outer (vertically thin) part of
 * an early-type "flat" bar, seen face-on or at low inclination
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * deltaPA_max = params[2 + offsetIndex];
 * I_0 = params[3 + offsetIndex ]; -- central intensity (ADU/pix)
 * h1 = params[4 + offsetIndex ];   -- inner exp. scale length (pixels)
 * h2 = params[5 + offsetIndex ];   -- outer exp. scale length (pixels)
 * r_break = params[7 + offsetIndex ];   -- break radius (pixels)
 * alpha = params[7 + offsetIndex ];     -- smoothness/sharpness of break
 *
 *
 */


// CLASS FlatBar:

#include "function_object.h"
#include <string>
#include <tuple>

using namespace std;



/// Class for image function with elliptical isophotes and FlatBar profile
class FlatBar : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    FlatBar( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool ParameterUnitsExist( );
   // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };

  protected:
    double CalculateIntensity( double r, double h2_adj, double r_b_adj );
    int  CalculateSubsamples( double r );
    // extra functions (not redefinitions):
    double  ParabolicInterpolation( double x, double x1, double x2, double y1, double y2 );
    std::tuple<double, double>  GetAdjustedRbh2( double xp, double yp, double r, double r_circ  );


  private:
    double  x0, y0, PA, ell, I_0, deltaPA_max, h1, h2, r_b, alpha;   // parameters
    double  q, PA_rad, cosPA, sinPA, deltaPA_max_rad, b_interp;   // other useful quantities
};

