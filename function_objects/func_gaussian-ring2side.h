/*   Class interface definition for func_gaussian-ring2side.cpp
 *   VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces an elliptical ring with a 2-sided Gaussian profile, which
 * allows for different Gaussian sigma for r < r_ring vs r > r_ring.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * ell = params[1 + offsetIndex ];   -- ellipticity of component
 * A = params[2 + offsetIndex ];   -- intensity scaling (ADU/pix)
 * R_ring = params[3 + offsetIndex ];
 * sig_r_in = params[4 + offsetIndex ];   -- Gaussian sigma in radial direction (inner side)
 * sig_r_out = params[5 + offsetIndex ];   -- Gaussian sigma in radial direction (outer side)
 *
 *
 */


// CLASS GaussianRing2Side:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "GaussianRing2Side"


class GaussianRing2Side : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    GaussianRing2Side( );
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
    double  x0, y0, PA, A, ell, R_ring, sigma_r_inner, sigma_r_outer;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful quantities
    double  twosigma_rin_squared, twosigma_rout_squared;
};
