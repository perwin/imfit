/*   Public interfaces for function(s) which deal with fitting via
 * differential evolution.
 */

#ifndef _DIFF_EVOLN_FIT_H_
#define _DIFF_EVOLN_FIT_H_

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"


// Note on possible return values for DiffEvolnFit: these are meant to be similar to
// the return values of LevMarFit and NMSimplexFit (see levmar_fit.h and nmsimplex_fit.h). 
//    value < 0   --> FAILURE
//    value = 0   --> FAILURE: input parameter error
//    value = 1   --> generic success
//    value = 5   --> max iterations reached

int DiffEvolnFit( int nParamsTot, double *initialParams, mp_par *parameterLimits, 
									ModelObject *theModel, double ftol, int verbose );


#endif  // _DIFF_EVOLN_FIT_H_
