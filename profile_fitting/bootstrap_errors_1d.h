/*   Public interfaces for function(s) which deal with estimating parameter
 * errors via bootstrap resampling
 */

#ifndef _BOOTSTRAP_ERRORS_1D_H_
#define _BOOTSTRAP_ERRORS_1D_H_

#include <vector>
#include "param_struct.h"   // for mp_par structure
#include "model_object.h"


void BootstrapErrors( double *bestfitParams, std::vector<mp_par> parameterLimits, 
						bool paramLimitsExist, ModelObject *theModel, double ftol,
						int nIterations, int nFreeParams );


#endif  // _BOOTSTRAP_ERRORS_1D_H_
