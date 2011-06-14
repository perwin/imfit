/*   Class interface definition for func_edge-on-ring.cpp
 *   VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces a crude edge-on ring model; this version allows for different
 * Gaussian sigma for r < r_ring vs r > r_ring.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * A = params[1 + offsetIndex ];   -- intensity scaling (ADU/pix)
 * r1 = params[2 + offsetIndex ];
 * r2 = -r1;
 * sig_r_in = params[3 + offsetIndex ];   -- Gaussian sigma in radial direction (inner side)
 * sig_r_out = params[4 + offsetIndex ];   -- Gaussian sigma in radial direction (outer side)
 * sig_z = params[5 + offsetIndex ];   -- Gaussian sigma in vertical direction
 *
 *
 */


// CLASS EdgeOnRing2Side:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "EdgeOnRing2Side"


class EdgeOnRing2Side : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    EdgeOnRing2Side( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r, double z );
    int  CalculateSubsamples( double r, double z );


  private:
    double  x0, y0, PA, A, r1, r2, sigma_r_inner, sigma_r_outer, sigma_z;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  twosigma_rin_squared, twosigma_rout_squared, twosigma_z_squared;
};
