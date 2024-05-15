/*   Class interface definition for func_peanut_dattathri.cpp
 * 
 *   A class derived from FunctionObject (function_object.h).
 *
 *
 */


// CLASS DattathriPeanut3D:

#include <string>
#include "gsl/gsl_integration.h"
#include "function_object.h"

using namespace std;



/// \brief Class for image function using LOS integration through 3D model
class DattathriPeanut3D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructor
    DattathriPeanut3D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, y0, PA, inc, psi_bar, J0_bar, R_bar_x, q, q_z, R_peanut, A_bar, 
    		sigma_peanut, c_bar_par, c_bar_perp;   // parameters
    double  PA_rad, cosPA, sinPA, inc_rad, cosInc, sinInc, cosPsi, sinPsi;   // other useful quantities
    double  R_bar_y, R_bar_z;   // other useful quantities
    gsl_function  F;
};

