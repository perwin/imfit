/*   Class interface definition for func_n4608disk.cpp
 *
 *   A class derived from FunctionObject (function_object.h), which produces
 * the intensity from a broken-exponential profile (BrokenExp component) added
 * to a Gaussian ring with azimuthally varying intensity (GaussianRingAz
 * component) -- intended for use in modeling WFC3-UVIS images of NGC 4608.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of BrokenExp component, rel. to +y axis
 * ell = params[1 + offsetIndex];  -- BrokenExp ellipticity
 * I_0 = params[2 + offsetIndex ]; -- BrokenExp: central intensity of inner exponential (ADU/pix)
 * h1 = params[3 + offsetIndex ];   -- BrokenExp: inner exp. scale length (pixels)
 * h2 = params[4 + offsetIndex ];   -- BrokenExp: outer exp. scale length (pixels)
 * r_break = params[5 + offsetIndex ];   -- BrokenExp: break radius (pixels)
 * alpha = params[6 + offsetIndex ];     -- BrokenExp: smoothness/sharpness of break [1/pixels]
 * PA_ring = params[7 + offsetIndex];   -- PA of ring, rel. to +y axis
 * ell_ring = params[8 + offsetIndex ];   -- ellipticity of ring
 * A_maj = params[9 + offsetIndex ];   -- intensity scaling RELATIVE TO I_0 [dimensionless]
 * A_min = params[10 + offsetIndex ];   -- intensity scaling RELATIVE TO I_0 [dimensionless]
 * R_ring = params[11 + offsetIndex ];  -- radius of ring (= center of Gaussian) along PA_ring direction
 * sigma = params[12 + offsetIndex ];   -- Gaussian sigma of ring in radial direction
 *
 */


// CLASS N4608Disk:

#include "function_object.h"
#include "func_broken-exp.h"
#include "func_gaussian-ring-az.h"


/// Class for image function with elliptical isophotes and broken-exponential profile
class N4608Disk : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    N4608Disk( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now
    void SetSubsampling( bool subsampleFlag );

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
  	double  x0, y0;   // general parameters: center of component
  	// BrokenExp object and parameters
  	BrokenExponential  funcBrokenExp;
    double  PA, ell, I_0, h1, h2, r_b, alpha;
    // GaussianRingAz object and parameters
  	GaussianRingAz  funcGaussianRingAz;
  	double  ringParams[6];
    double  PA_ring, A_maj_rel, A_min_rel, ell_ring, R_ring, sigma_r;
};
