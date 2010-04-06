/*   Public interfaces for function(s) which deal with fitting via
 * differential evolution.
 */

#ifndef _DIFF_EVOLN_FIT_H_
#define _DIFF_EVOLN_FIT_H_

#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "model_object.h"


int DiffEvolnFit(int nParamsTot, double *initialParams, mp_par *parameterLimits, 
									ModelObject *theModel, int maxGenerations);


#endif  // _DIFF_EVOLN_FIT_H_
