/*   Class interface definition for func_gaussian.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a circular
 * Gaussian.  This will generate a normalized profile (total flux = 1.0).
 *
 */


#ifndef _FUNC_GAUSSIAN_H_
#define _FUNC_GAUSSIAN_H_

// CLASS Gaussian:

#include "function_object.h"
#include <string>



/// Class for image function with elliptical isophotes and Gaussian profile
class Gaussian : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Gaussian( );
    // redefined method/member function:
    void AdjustParametersForImage( const double inputFunctionParams[], 
									double adjustedFunctionParams[], int offsetIndex );
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool CanCalculateTotalFlux(  );
    double TotalFlux( );
    // No destructor for now
    
    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, I_0, sigma;   // parameters
    double  twosigma_squared;
    double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};


#endif /* _FUNC_GAUSSIAN_H_ */
