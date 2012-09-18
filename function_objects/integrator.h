/*   Public interfaces for function(s) which takes a list of user-specified
 * function objects and adds them to the ModelObject.
 */

#ifndef _GSL_INTEGRATOR_H_
#define _GSL_INTEGRATOR_H_

#include "gsl/gsl_integration.h"

using namespace std;


double  Integrate( gsl_function F, double s1, double s2 );


#endif  // _GSL_INTEGRATOR_H_
