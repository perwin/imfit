/*   Class interface definition for func1d_spline.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of a cubic spline interpolation.
 *
 * PARAMETERS:
 * I_0 = params[0 + offsetIndex ]; -- central surf. brightness (mag/arcsec^2)
 * r_1 = params[1 + offsetIndex ];    -- radius of first interpolation data point
 * I_1 = params[1 + offsetIndex ];    -- surf. brightness of first interpolation data point
 * r_2 = params[1 + offsetIndex ];    -- radius of first interpolation data point
 * I_2 = params[1 + offsetIndex ];    -- surf. brightness of first interpolation data point
 *
 *
 */


// CLASS Spline1D:
#include <gsl/gsl_spline.h>

#include "function_object.h"

// maximum number of data points (always >= 2)
const int  MAX_POINTS = 4;


class Spline1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Spline1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // Destructor (handles deallocation of GSL spline structures)
    ~Spline1D( );

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, I_0, r_1, I_1;   // parameters
    double  xInterp[MAX_POINTS];
    double  yInterp[MAX_POINTS];
    int  nInterpPoints;
    bool  splineFuncAllocated;
    bool  splineCacheAllocated;
    gsl_interp_accel *splineCache;
    gsl_spline *splineFunc;
};
