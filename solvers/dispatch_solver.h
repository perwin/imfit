/** @file
 * \brief Code for calling specific solvers (minimization algorithms)
 *
 * 
 *
 */

#ifndef _DISPATCH_SOLVER_H_
#define _DISPATCH_SOLVER_H_

#include <string>
#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "solver_results.h"


/// Function which handles selecting and calling appropriate solver
int DispatchToSolver( int solverID, int nParametersTot, int nFreeParameters, int nPixelsTot,
					double *parameters, vector<mp_par> parameterInfo, ModelObject *modelObj, double fracTolerance,
					bool paramLimitsExist, int verboseLevel, SolverResults *solverResults,
					string& solverName, unsigned long rngSeed=0 );


#endif /* _DISPATCH_SOLVER_H_ */
