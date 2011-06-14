/*   Class interface definition for func_edge-on-disk_n4762.cpp
 *   VERSION 0.2
 *
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces an alternate form of edge-on disk, where the vertical profile
 * is exponential and the horizontal profile is a flat-exponential, with a break
 * radius which is determined by a linear function of z (r_b = a + b*z).
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +y axis
 * I_0 = params[1 + offsetIndex ]; -- central intensity scaling (ADU/pix)
 * h2 = params[2 + offsetIndex];  -- horizontal exp. scale length outside break radius
 * a_rb = params[3 + offsetIndex];  -- intercept for horizontal break radius
 * b_rb = params[4 + offsetIndex];  -- slope for horizontal break radius
 * alpha = params[5 + offsetIndex ];   -- sharpness of break
 * h_z = params[6 + offsetIndex ] -- vertical exp. scale height
 *
 *
 */


// CLASS EdgeOnDiskN4762v2:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "EdgeOnDisk_n4762v2"


class EdgeOnDiskN4762v2 : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    EdgeOnDiskN4762v2( );
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
    double  x0, y0, PA, I_0, h2, alpha, a_rb, b_rb, h_z;   // parameters
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    double  I_0_times_S, exponent, r_b, delta_Rb_scaled;
};
