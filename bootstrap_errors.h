/*   Public interfaces for function(s) which deal with estimating parameter
 * errors via bootstrap resampling
 */

#ifndef _BOOTSTRAP_ERRORS_H_
#define _BOOTSTRAP_ERRORS_H_

#include <string>

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"


void BootstrapErrors( double *bestfitParams, mp_par *parameterLimits, bool paramLimitsExist, 
					ModelObject *theModel, double ftol, int nIterations, int nFreeParams,
					int whichStatistic, double **outputParamArray );


#endif  // _BOOTSTRAP_ERRORS_H_
