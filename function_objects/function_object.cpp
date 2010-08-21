/* FILE: function_object.cpp ------------------------------------------- */
/* VERSION 0.25
 *
 *   This is the base class for the various function object classes.
 *   It really shouldn't be instantiated by itself.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *      So for a vector of x[], we expect the user to call Setup() once at
 *      the start, then loop through x[], calling GetValue() to compute
 *      the function results for each x.
 *
 *   MODIFICATION HISTORY:
 *     [v0.25]: 5 Dec: Minor cleanup.
 *     [v0.2]: 26---27 Nov: Modifications to accomodate derived classes
 * for working with 1-D fitting; also changed uses of char * to uses of
 * C++ string class.
 *     [v0.01]: 13 Nov 2009: Created (as modification of nonlinfit2's
 *   function_object class.
 */


/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <string.h>
#include <string>

#include "function_object.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */


/* ---------------- CONSTRUCTOR ---------------------------------------- */

FunctionObject::FunctionObject( )
{
  functionName = "Base (undefined) function";
  shortFunctionName = "BaseFunction";
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */
// Base method for 2D functions
void FunctionObject::Setup( double params[], int offsetIndex, double xc, double yc )
{
  // This method is for passing the current parameters into the function
  // object, storing them for when GetValue() is called, and for pre-computing
  // various useful quantities.
  // The parameter array params contains *all* parameters for *all* components
  // in the overall model; offsetIndex is used to select the correct starting point
  // for *this* component's parameters.
  ;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */
// Base method for 1D functions
void FunctionObject::Setup( double params[], int offsetIndex, double xc )
{
  // This method is for passing the current parameters into the function
  // object, storing them for when GetValue() is called, and for pre-computing
  // various useful quantities.
  // The parameter array params contains *all* parameters for *all* components
  // in the overall model; offsetIndex is used to select the correct starting point
  // for *this* component's parameters.
  ;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Base method for 2D functions
double FunctionObject::GetValue( double x, double y )
{
  // This method is for computing the actual function value at the 
  // pixel coordinates (x,y).  It will be called once per pixel.
  return -1.0;   // dummy (bad) value indicating the user has accidentally
                 // accessed the base class, which should never happen...
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// Base method for 1D functions
double FunctionObject::GetValue( double x )
{
  // This method is for computing the actual function value at for the
  // specified independent variable value x.
  return -1.0;   // dummy (bad) value indicating the user has accidentally
                 // accessed the base class, which should never happen...
}


/* ---------------- PUBLIC METHOD: GetDescription ---------------------- */

string FunctionObject::GetDescription( )
{
  return functionName;
}


/* ---------------- PUBLIC METHOD: GetShortName ----------------------- */

string& FunctionObject::GetShortName( )
{
  return shortFunctionName;
}


/* ---------------- PUBLIC METHOD: GetParameterNames ------------------- */
// Add this function's parameter names to a vector of strings
void FunctionObject::GetParameterNames( vector<string> &paramNameList )
{
  for (int i = 0; i < nParams; i++)
    paramNameList.push_back(parameterLabels[i]);
}


/* ---------------- PUBLIC METHOD: GetNParams -------------------------- */

int FunctionObject::GetNParams( )
{
  return nParams;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

FunctionObject::~FunctionObject()
{
  ;   // nothing to do in base class
}



/* END OF FILE: function_object.cpp ------------------------------------ */
