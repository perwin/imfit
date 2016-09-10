/*   Class interface definition for func1d_gaussian_linear.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D Gaussian.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ];   -- central surf. brightness (mag/arcsec^2)
 * sigma = params[1 + offsetIndex ];  -- sigma of the Gaussian
 *
 *
 */


// CLASS LinearGaussian1D:

#include "function_object.h"


class LinearGaussian1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    LinearGaussian1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, I_0, sigma;   // parameters
};
