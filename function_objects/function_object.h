/*   Class interface definition for function_object.cpp [imfit]
 *
 * This is intended to be an abstract base class for the various
 * function objects (e.g., Sersic function, broken-exponential
 * function, etc.).
 * 
 */


// CLASS FunctionObject [base class]:

#ifndef _FUNCTION_OBJ_H_
#define _FUNCTION_OBJ_H_

#include <map>
#include <string>
#include <vector>

using namespace std;


/// Virtual base class for function objects (i.e., 2D image functions)
class FunctionObject
{
  public:
    // Constructors:
    FunctionObject( );

    // override in derived classes only if said class *does* have user-settable
    // extra parameters
    virtual bool HasExtraParams( ) { return(false); }
    virtual int SetExtraParams( map<string, string>& )  { return -1; }
    virtual bool ExtraParamsSet( ) { return(extraParamsSet); }

    // override in derived classes only if said class is PointSource or similar
    virtual bool IsPointSource( ) { return(false); };
    virtual void AddPsfData( double *psfPixels, int nColumns_psf, int nRows_psf ) { ; };

    // probably no need to modify this:
    virtual void SetSubsampling( bool subsampleFlag );

    // probably no need to modify this:
    virtual void SetZeroPoint( double zeroPoint );

    // derived classes will almost certainly modify this, which
    // is used for pre-calculations and convolutions, if any:
    virtual void Setup( double params[], int offsetIndex, double xc, double yc );

    // all derived classes working with 1D data must override this:
    virtual void Setup( double params[], int offsetIndex, double xc );

    // all derived classes working with 2D images must override this:
    virtual double GetValue( double x, double y );

    // all derived classes working with 1D data must override this:
    virtual double GetValue( double x );

    // override in derived classes only if said class *can* calcluate total flux
    virtual bool CanCalculateTotalFlux(  ) { return(false); }
    // override in derived classes only if said class *can* calcluate total flux
    virtual double TotalFlux( ) { return -1.0; }

    // no need to modify this:
    virtual string GetDescription( );

    // no need to modify this:
    virtual string& GetShortName( );

    // no need to modify this:
    virtual void GetParameterNames( vector<string> &paramNameList );

    // no need to modify this:
    virtual int GetNParams( );

    // Destructor (doesn't have to be modified, but MUST be declared
    // virtual in order for this to be a sensible base object
    // [see e.g. Scott Meyers, Effective C++]; otherwise behavior is 
    // undefined when a derived class is deleted)
    virtual ~FunctionObject();


  private:
  
  protected:  // same as private, except accessible to derived classes
    int  nParams;
    bool  doSubsampling;
    bool  extraParamsSet;
    vector<string>  parameterLabels;
    string  functionName, shortFunctionName;
    double  ZP;

    // class member (constant char-vector string) which will hold name of
    // individual class in derived classes
    static const char  shortFuncName[];  ///< Class data member holding name of individual class
  
};

#endif   // _FUNCTION_OBJ_H_
