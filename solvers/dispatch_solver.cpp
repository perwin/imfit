/** @file
 * \brief DispatchToSolver function.
 */
/* FILE: dispatch_solver.cpp --------------------------------------------- */

// Copyright 2010--2018 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.


// * Return values:
// 1 = Generic success return value (e.g., fit-statistic convergence)
// 5 = Optimization stopped because maximum number of generations was reached
// Error codes (negative return values)
// -1  = Generic failure code [NOT USED YET]
// -2 = Missing fitting bounds for at least one parameter.
// -5 = Halted because user objective function returned NaN.

#include <stdio.h>

#include "definitions.h"
#include "model_object.h"
#include "param_struct.h"   // for vector<mp_par>
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
					double *parameters, vector<mp_par> parameterInfo, ModelObject *modelObj, 
					double fracTolerance, bool paramLimitsExist, int verboseLevel, 
					SolverResults *solverResults, string& solverName, 
					unsigned long rngSeed, bool useLHS )
{
  int  fitStatus = -100;
  
  switch (solverID) {
    case MPFIT_SOLVER:
      if (verboseLevel >= 0)
        printf("Calling Levenberg-Marquardt solver ...\n");
      fitStatus = LevMarFit(nParametersTot, nFreeParameters, nPixelsTot, parameters, parameterInfo, 
      						modelObj, fracTolerance, paramLimitsExist, verboseLevel, solverResults);
      break;
    case DIFF_EVOLN_SOLVER:
      if (verboseLevel >= 0)
        printf("Calling Differential Evolution solver ..\n");
      fitStatus = DiffEvolnFit(nParametersTot, parameters, parameterInfo, modelObj, fracTolerance, 
      							verboseLevel, solverResults, rngSeed, useLHS);

      break;
#ifndef NO_NLOPT
    case NMSIMPLEX_SOLVER:
      if (verboseLevel >= 0)
        printf("Calling Nelder-Mead Simplex solver ..\n");
      fitStatus = NMSimplexFit(nParametersTot, parameters, parameterInfo, modelObj, fracTolerance, 
      							verboseLevel, solverResults);
      break;
    case GENERIC_NLOPT_SOLVER:
      if (verboseLevel >= 0)
        printf("\nCalling NLOpt solver %s ..\n", solverName.c_str());
      fitStatus = NLOptFit(nParametersTot, parameters, parameterInfo, modelObj, fracTolerance, 
      						verboseLevel, solverName, solverResults);
      break;
#endif
  }

  return fitStatus;
}


/* END OF FILE: dispatch_solver.cpp -------------------------------------- */
