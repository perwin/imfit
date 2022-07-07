/*   Class interface definition for func_pointsource.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces a shifted and interpolated copy of a PSF image.
 *
 */


// CLASS PointSource:

#include "function_object.h"
#include "psf_interpolators.h"
#include <string>



/// Class for image function with shifted and interpolated PSF image
class PointSource : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    PointSource( );
    // Need a destructor to dispose of PsfInterpolator object
    ~PointSource( );

    // redefined method/member function:
    bool IsPointSource( );
    void AddPsfInterpolator( PsfInterpolator *theInterpolator );
    void AddPsfData( double *psfPixels, int nColumns_psf, int nRows_psf );
    string GetInterpolationType( );
    void SetOversamplingScale( int oversampleScale );
    bool HasExtraParams( );
    int SetExtraParams( map<string, string>& inputMap );

    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool CanCalculateTotalFlux(  );
    double TotalFlux( );
    bool ParameterUnitsExist( );
    
    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };



  private:
    double  x0, y0, I_tot;   // parameters
    double  I_scaled;
    string  interpolationType = "bicubic";
    int  oversamplingScale;
    PsfInterpolator *psfInterpolator;
    bool interpolatorAllocated = false;
};
