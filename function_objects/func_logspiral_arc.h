/*   Class interface definition for func_logspiral_arc.cpp
 *   VERSION 0.1
 *
 *   Class (derived from FunctionObject; function_object.h) which produces a simple/crude
 * 1-arm logarithmic spiral with overall amplitude modulated by azimuthal, 2-sided 
 * Gaussian function
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   // line of nodes for projected circle [deg]
 * ell = params[1 + offsetIndex];  // ellipticity of projected circle
 * r_scale = params[2 + offsetIndex ];   // radial size of spiral arc [pix]
 * i_pitch = params[3 + offsetIndex ];   // pitch or winding angle [degrees]
 * theta_max = params[4 + offsetIndex];  // angle of maximum intensity [deg]
 * I_max = params[5 + offsetIndex];   // maximum intensity of spiral arc
 * sigma_r = params[6 + offsetIndex];  // Gaussian width of spiral in radial direction [pix]
 * sigma_theta_ccw = params[7 + offsetIndex];  // sigma (in CCW direction) for azimuthal Gaussian [deg]
 * sigma_theta_cw = params[8 + offsetIndex];  // sigma (in CW direction) for azimuthal Gaussian [deg]
 *
 */


// CLASS LogSpiralArc:

#include "function_object.h"
#include <tuple>

//#define CLASS_SHORT_NAME  "LogSpiralArc"


class LogSpiralArc : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    LogSpiralArc( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r, double theta_rel, double theta_std );
    int  CalculateSubsamples( double r );
    double GetStdAngle( double x_diff, double y_diff );
    double  GetThetaRel( double theta_std );
    std::tuple<double, double>  GetThetas( double x_diff, double y_diff );
    double  GetThetaDiff( double theta_std );


  private:
    double  x0, y0, PA, ell, r_scale, i_pitch, theta_max, I_max, sigma_r;
    double  sigma_theta_ccw, sigma_theta_cw;
    double  q, PA_rad, cosPA, sinPA, tan_i_pitch, theta_max_rad, twosigma_r_squared;
    double  sigma_theta_ccw_rad, sigma_theta_cw_rad;
    double  twosigma_theta_ccw_squared, twosigma_theta_cw_squared, theta_sp_ref, R0;
};
