/*   Public interfaces for function(s) which deal with fitting via Nelder-Mead Simplex
 */

#ifndef _LEVMAR_FIT_H_
#define _LEVMAR_FIT_H_

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"


int LevMarFit( int nParamsTot, int nFreeParams, int nDataVals, double *paramVector, 
				mp_par *parameterLimits, ModelObject *theModel, double ftol, 
				bool paramLimitsExist, int verbose );


#endif  // _LEVMAR_FIT_H_
