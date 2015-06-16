/*   Public interfaces for function(s) which deal with fitting via Nelder-Mead Simplex
 */

#ifndef _LEVMAR_FIT_H_
#define _LEVMAR_FIT_H_

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"
#include "solver_results.h"


// Note on possible return values for LevMarFit: these are the same as the return values
// of mpfit (see mpfit_cpp.h). In summary:
//    values <= 0: error of some kind
//    values = 1--4: general convergence success of different types
//    value = 5: max number of iterations
//    value = 6--8: ftol,xtol,gtol too small, no further improvement possible
// #define MP_OK_CHI (1)            /* Convergence in chi-square value */
// #define MP_OK_PAR (2)            /* Convergence in parameter value */
// #define MP_OK_BOTH (3)           /* Both MP_OK_PAR and MP_OK_CHI hold */
// #define MP_OK_DIR (4)            /* Convergence in orthogonality */
// #define MP_MAXITER (5)           /* Maximum number of iterations reached */
// #define MP_FTOL (6)              /* ftol is too small; no further improvement*/
// #define MP_XTOL (7)              /* xtol is too small; no further improvement*/
// #define MP_GTOL (8)              /* gtol is too small; no further improvement*/

int LevMarFit( int nParamsTot, int nFreeParams, int nDataVals, double *paramVector, 
				mp_par *parameterLimits, ModelObject *theModel, const double ftol, 
				const bool paramLimitsExist, const int verbose, SolverResults *solverResults=0 );


#endif  // _LEVMAR_FIT_H_
