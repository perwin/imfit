/*   Class interface definition for func_gaussian-ring-az.cpp
 * 
 *   Class (derived from FunctionObject; function_object.h) which produces an elliptical 
 * ring with a Gaussian radial profile, with an azimuthal variation of the Gaussian
 * amplitude.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * ell = params[1 + offsetIndex ];   -- ellipticity of component
 * A_maj = params[2 + offsetIndex ];   -- intensity scaling (ADU/pix)
 * A_min = params[3 + offsetIndex ];   -- intensity scaling along ellipse major axis (ADU/pix)
 * R_ring = params[4 + offsetIndex ];  -- intensity scaling along ellipse minor axis (ADU/pix)
 * sigma = params[5 + offsetIndex ];   -- Gaussian sigma in radial direction
 *
 *
 */


#ifndef _FUNC_GAUSSIANRINGAZ_H_
#define _FUNC_GAUSSIANRINGAZ_H_


// CLASS GaussianRingAz:

#include "function_object.h"



/// \brief Class for image function using 2D elliptical ring with Gaussian radial profile
/// and azimuthal variation of Gaussian amplitude
class GaussianRingAz : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    GaussianRingAz( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r, double theta );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, A_maj, A_min, ell, R_ring, sigma_r;   // parameters
    double  q, PA_rad, cosPA, sinPA, A_mid, deltaA;   // other useful quantities
    double  twosigma_squared;
};

#endif   // _FUNC_GAUSSIANRINGAZ_H_
