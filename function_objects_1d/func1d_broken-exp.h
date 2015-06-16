/*   Class interface definition for func1d_broken-exp.cpp
 *   VERSION 0.2
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * exponential.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ];    -- central surf. brightness (mag/arcsec^2)
 * h_1 = params[1 + offsetIndex ];    -- inner exp. scale length
 * h_2 = params[2 + offsetIndex ];    -- outer exp. scale length
 * r_b = params[3 + offsetIndex ];    -- break radius
 * alpha = params[4 + offsetIndex ];  -- sharpness of break
 *
 *
 */


// CLASS BrokenExponential1D:

#include "function_object.h"


class BrokenExponential1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    BrokenExponential1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, mu_0, h_1, h_2, r_b, alpha;   // parameters
    double  I_0, exponent, S;   // other useful quantities
};
