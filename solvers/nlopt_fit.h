/** @file
 * \brief Public functions for setting up and running various minimizers from NLOpt lirary
 *
 *   Public interfaces for function(s) which deal with fitting via the various
 * algorithms in the NLopt library.
 */

#ifndef _NM_NLOPT_FIT_H_
#define _NM_NLOPT_FIT_H_

#include <string>
#include "param_struct.h"   // for mp_par structure
#include "model_object.h"
#include "solver_results.h"

using namespace std;

// Note on possible return values for NMSimplexFit: these are the same as the return
// values of nlopt_optimize (see nlopt.h):
//      NLOPT_FAILURE = -1, /* generic failure code */
//      NLOPT_INVALID_ARGS = -2,
//      NLOPT_OUT_OF_MEMORY = -3,
//      NLOPT_ROUNDOFF_LIMITED = -4,
//      NLOPT_FORCED_STOP = -5,
//      NLOPT_SUCCESS = 1, /* generic success code */
//      NLOPT_STOPVAL_REACHED = 2,
//      NLOPT_FTOL_REACHED = 3,
//      NLOPT_XTOL_REACHED = 4,
//      NLOPT_MAXEVAL_REACHED = 5,
//      NLOPT_MAXTIME_REACHED = 6

bool ValidNLOptSolverName( string solverName );

int NLOptFit( int nParamsTot, double *initialParams, vector<mp_par> parameterLimits, 
					ModelObject *theModel, double ftol, int verbose, string solverName,
					SolverResults *solverResults=0 );

void GetInterpretation_NLOpt( int resultValue, string& outputString );


#endif  // _NM_NLOPT_FIT_H_
