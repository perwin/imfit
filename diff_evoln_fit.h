/*   Public interfaces for function(s) which deal with fitting via
 * differential evolution.
 */

#ifndef _DIFF_EVOLN_FIT_H_
#define _DIFF_EVOLN_FIT_H_

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"


int DiffEvolnFit( int nParamsTot, double *initialParams, mp_par *parameterLimits, 
									ModelObject *theModel, double ftol, int verbose );


#endif  // _DIFF_EVOLN_FIT_H_
