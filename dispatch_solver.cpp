#include <stdio.h>

#include "definitions.h"
#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "solver_results.h"
#include "dispatch_solver.h"

// Solvers (optimization algorithms)
#include "levmar_fit.h"
#include "diff_evoln_fit.h"
#ifndef NO_NLOPT
#include "nmsimplex_fit.h"
#include "nlopt_fit.h"
#endif

int DispatchToSolver( int solverID, int nParametersTot, int nFreeParameters, int nPixelsTot,
					double *parameters, mp_par *parameterInfo, ModelObject *modelObj, double fracTolerance,
					bool paramLimitsExist, int verboseLevel, SolverResults *solverResults,
					string& solverName )
{
  int  fitStatus = -100;
  
  switch (solverID) {
    case MPFIT_SOLVER:
      printf("Calling Levenberg-Marquardt solver ...\n");
      fitStatus = LevMarFit(nParametersTot, nFreeParameters, nPixelsTot, parameters, parameterInfo, 
      						modelObj, fracTolerance, paramLimitsExist, verboseLevel, solverResults);
      break;
    case DIFF_EVOLN_SOLVER:
      printf("Calling Differential Evolution solver ..\n");
      fitStatus = DiffEvolnFit(nParametersTot, parameters, parameterInfo, modelObj, fracTolerance, 
      							verboseLevel, solverResults);

      break;
#ifndef NO_NLOPT
    case NMSIMPLEX_SOLVER:
      printf("Calling Nelder-Mead Simplex solver ..\n");
      fitStatus = NMSimplexFit(nParametersTot, parameters, parameterInfo, modelObj, fracTolerance, 
      							verboseLevel, solverResults);
      break;
    case GENERIC_NLOPT_SOLVER:
      printf("\nCalling NLOpt solver %s ..\n", solverName.c_str());
      fitStatus = NLOptFit(nParametersTot, parameters, parameterInfo, modelObj, fracTolerance, 
      						verboseLevel, solverName, solverResults);
      break;
#endif
  }

  return fitStatus;
}



// in main():
// 
//   fitStatus = DispatchSolver(options.solver, nParamsTot, nFreeParams, nPixels_tot, paramsVect,
//   							 parameterInfo, theModel, options.ftol, paramLimitsExist,
// 							 options.verbose, &resultsFromSolver, options.nloptSolverName);
// 
//   PrintResults(paramsVect, 0, theModel, nFreeParams, parameterInfo, fitStatus);
// 
