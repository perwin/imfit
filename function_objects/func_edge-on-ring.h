/*   Class interface definition for func_edge-on-ring.cpp
 *   VERSION 0.1
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces a crude edge-on ring model.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * A = params[1 + offsetIndex ];   -- intensity scaling (ADU/pix)
 * r1 = params[2 + offsetIndex ];
 * r2 = -r1;
 * sig_r = params[3 + offsetIndex ];   -- Gaussian sigma in radial direction
 * sig_z = params[4 + offsetIndex ];   -- Gaussian sigma in vertical direction
 *
 *
 */


// CLASS EdgeOnRing:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "EdgeOnRing"


class EdgeOnRing : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    EdgeOnRing( );
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
    double  x0, y0, PA, A, r1, r2, sigma_r, sigma_z;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  twosigma_r_squared, twosigma_z_squared;
};
