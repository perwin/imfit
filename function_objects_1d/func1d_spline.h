/*   Class interface definition for func1d_spline.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of a cubic spline interpolation.
 *
 * PARAMETERS:
 * PARAM_LABELS[][20] = {"I_0", "r_1", "I_1", "r_2", "I_2", "r_3", "I_3"};
 * I_0 = params[0 + offsetIndex ]; -- central surf. brightness (1st interpolation data point)
 * r_1 = params[1 + offsetIndex ];    -- radius of 2nd interpolation data point
 * I_1 = params[2 + offsetIndex ];    -- surf. brightness of 2nd interpolation data point
 * r_2 = params[3 + offsetIndex ];    -- radius of 3rd interpolation data point
 * I_2 = params[4 + offsetIndex ];    -- surf. brightness of 3rd interpolation data point
 * r_3 = params[5 + offsetIndex ];    -- radius of 4th interpolation data point
 * I_3 = params[6 + offsetIndex ];    -- surf. brightness of 4th interpolation data point
 *
 * Note that this function includes one optional/extra parameter:
      OPTIONAL_PARAMS_START
      N    100.0
      OPTIONAL_PARAMS_END
 * where N = maximum number of interpolation points for the spline interpolation.
 * This lets us specify that the interpolation should only use e.g. 2 or 3 points rather
 * than the full 4 possible points (the config file must still include lines for all
 * possible interpolation points)
 *
 */


// CLASS Spline1D:
#include <gsl/gsl_spline.h>

#include "function_object.h"

// maximum number of data points for interpolation function (always >= 2)
const int  MAX_SPLINE_POINTS = 4;


class Spline1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Spline1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    bool HasExtraParams( );
    int SetExtraParams( map<string, string>& inputMap );
    double  GetValue( double x );
    // Destructor (handles deallocation of GSL spline structures)
    ~Spline1D( );

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, I_0, r_1, I_1;   // parameters
    double  xInterp[MAX_SPLINE_POINTS];
    double  yInterp[MAX_SPLINE_POINTS];
    int  nInterpPoints, maxInterpPoints;
    bool  splineFuncAllocated;
    bool  splineCacheAllocated;
    gsl_interp_accel *splineCache;
    gsl_spline *splineFunc;
};
