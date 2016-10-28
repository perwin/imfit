/*   Class interface definition for func1d_exp_test.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a 1-D
 * exponential.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ]; -- central surf. brightness (mag/arcsec^2)
 * h = params[1 + offsetIndex ];    -- exp. scale length
 *
 *
 */


// CLASS Exponential1D_test:
#include <map>
#include <string>
using namespace std;

#include "function_object.h"


class Exponential1D_test : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Exponential1D_test( );
    // redefined method/member functions:
    bool HasExtraParams( );
    int  SetExtraParams( map<string, string> inputMap );
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, mu_0, h;   // parameters
    double  I_0, floorValue;   // other useful quantities
};
