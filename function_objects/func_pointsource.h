/*   Class interface definition for func_pointsource.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces XXX
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
    void AddPsfInterpolator( PsfInterpolator *theInterpolator );
    bool IsPointSource( );
    void AddPsfData( double *psfPixels, int nColumns_psf, int nRows_psf );
    bool HasExtraParams( );
    int SetExtraParams( map<string, string>& inputMap );

    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    bool CanCalculateTotalFlux(  );
    double TotalFlux( );
    
    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };



  private:
    double  x0, y0, I_tot;   // parameters
    double  I_scaled;
    string  interpolatorType = "bicubic";
    PsfInterpolator *psfInterpolator;
    bool interpolatorAllocated = false;
};
