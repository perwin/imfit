/*   Class interface definition for func1d_gauss-hermite.cpp
 *   VERSION 0.3
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces y-values as a function of x for the sum of two
 * Gauss-Hermite functions.
 *
 * PARAMETERS:
 * A = params[0 + offsetIndex ]; -- relative weight of G-H component
 * sigma = params[1 + offsetIndex ];    -- sigma of G-H component
 * h3 = params[2 + offsetIndex ];    -- h_3 of G-H component
 * h4 = params[3 + offsetIndex ];    -- h_4 of G-H component
 *
 */


// CLASS GaussHermite1D:

#include "function_object.h"


class GaussHermite1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    GaussHermite1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, A, sigma, h3, h4;   // parameters
};
