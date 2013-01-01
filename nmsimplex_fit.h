/*   Public interfaces for function(s) which deal with fitting via Nelder-Mead Simplex
 */

#ifndef _NM_SIMPLEX_FIT_H_
#define _NM_SIMPLEX_FIT_H_

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"


int NMSimplexFix(int nParamsTot, double *initialParams, mp_par *parameterLimits, 
									ModelObject *theModel, double ftol );


#endif  // _NM_SIMPLEX_FIT_H_
