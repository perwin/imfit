/*   Class interface definition for func_nuker.cpp
 * 
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * Nuker-Law profile (Lauer+1995, Byun+1996) function.
 *
 * PARAMETERS:
 * PA = params[0 + offsetIndex];      -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];     -- ellipticity
 * I_b = params[2 + offsetIndex ];    -- break-radius surf. brightness (counts/pixel)
 * r_b = params[3 + offsetIndex ];    -- break radius (pixels)
 * alpha = params[4 + offsetIndex ];  -- sharpness of break
 * beta = params[5 + offsetIndex ];   -- outer power-law slope
 * gamma = params[6 + offsetIndex ];  -- inner power-law slope
 *
 *
 */


// CLASS NukerLaw:

#include "function_object.h"



/// \brief Class for image function with elliptical isophotes and Core-Sersic profile
class NukerLaw : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];

  public:
    // Constructors:
    NukerLaw( );
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
  double  x0, y0, PA, ell, n, I_b, r_b, alpha, beta, gamma;   // parameters
  double  Iprime, exponent;
  double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};
