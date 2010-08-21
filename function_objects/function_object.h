/*   Class interface definition for function_object.cpp [imfit]
 *   VERSION 0.2
 *
 * This is intended to be an abstract base class for the various
 * function objects (e.g., Sersic function, broken-exponential
 * function, etc.).
 * 
 */


// CLASS FunctionObject [base class]:

#ifndef _FUNCTION_OBJ_H_
#define _FUNCTION_OBJ_H_

#include <string>
#include <vector>

#include "definitions.h"

using namespace std;


class FunctionObject
{
  public:
    // Constructors:
    FunctionObject( );

    // derived classes will almost certainly modify this, which
    // is used for pre-calculations and convolutions, if any:
    virtual void  Setup( double params[], int offsetIndex, double xc, double yc );

    // all derived classes working with 1D data must override this:
    virtual void  Setup( double params[], int offsetIndex, double xc );

    // all derived classes working with 2D images must override this:
    virtual double GetValue( double x, double y );

    // all derived classes working with 1D data must override this:
    virtual double GetValue( double x );

    // no need to modify this:
    virtual string GetDescription( );

    // no need to modify this:
    virtual string& GetShortName( );

    // no need to modify this:
    virtual void GetParameterNames( vector<string> &paramNameList );

    // no need to modify this:
    virtual int GetNParams( );

    // Destructor (doesn't have to be modified):
    virtual ~FunctionObject();


  private:
  
  protected:  // same as private, except accessible to derived classes
    int  nParams;
    bool  doSubsampling;
    vector<string>  parameterLabels;
    string  functionName, shortFunctionName;
  
};

#endif   // _FUNCTION_OBJ_H_
