/*   Class interface definition for func1d_broken-exp.cpp
 *   VERSION 0.2
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * exponential.
 *
 * PARAMETERS:
 * mu_0 = params[0 + offsetIndex ];    -- central surf. brightness (mag/arcsec^2)
 * h_1 = params[1 + offsetIndex ];    -- inner exp. scale length
 * h_2 = params[2 + offsetIndex ];    -- outer exp. scale length
 * r_b = params[3 + offsetIndex ];    -- break radius
 * alpha = params[4 + offsetIndex ];  -- sharpness of break
 *
 *
 */


// CLASS n1543MajMinCircBulge1D:

#include "function_object.h"


class n1543MajMinCircBulge1D : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    n1543MajMinCircBulge1D( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc );
    double  GetValue( double x );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  private:
    double  x0, n, mu_e, r_e, mu_0_bar, h_bar1, h_bar2, r_b, alpha;   // parameters
    double  mu_0_disk, h_disk_maj, h_disk_min, sigma_bar;   // more parameters
    double  I_e, n2, bn, invn, I_0_disk, I_0_bar, exponent, S;   // other useful quantities

	double  Gaussian( double x, double I_0, double sigma );
	double  Exponential( double x, double I_0, double h );
	double  Sersic( double x, double n, double Ie, double re );

};

