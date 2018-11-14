/* FILE: levmar_fit.cpp -------------------------------------------------- */

// Copyright 2012-2018 by Peter Erwin.
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



// Note : the following are the default tolerance values we are currently using
// in mpfitfun.cpp:
//  conf.ftol = 1e-10;   [relative changes in chi^2]
//  conf.xtol = 1e-10;   [relative changes in parameter values]

#include <string.h>   // for memset
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "mpfit.h"
#include "print_results.h"
#include "solver_results.h"

const int  MAX_ITERATIONS = 1000;
const double  FTOL = 1.0e-8;
const double  XTOL = 1.0e-8;


/* ------------------- Function Prototypes ----------------------------- */

int myfunc_mpfit( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );




/// This is the function used by mpfit() to compute the vector of deviates.
/// In our case, it's a wrapper which tells the ModelObject to compute 
/// and return the deviates.
int myfunc_mpfit( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *theModel )
{

  theModel->ComputeDeviates(deviates, params);
  return 0;
}





int LevMarFit( int nParamsTot, int nFreeParams, int nDataVals, double *paramVector, 
				vector<mp_par> parameterLimits, ModelObject *theModel, const double ftol, 
				const bool paramLimitsExist, const int verbose, SolverResults *solverResults )
{
  double  *paramErrs;
  mp_par  *mpfitParameterConstraints;
  bool  parameterConstraintsAllocated = false;
  mp_result  mpfitResult;
  mp_config  mpConfig;
  int  status;


  // Since we now use vector<mp_par> in main and elsewhere, we need to allocate
  // and construct a corresponding mp_par * array, if parameter limits actually exist
  if (! paramLimitsExist) {
    // If parameters are unconstrained, then mpfit() expects a NULL mp_par array
    mpfitParameterConstraints = NULL;
  } else {
    mpfitParameterConstraints = (mp_par *) calloc((size_t)nParamsTot, sizeof(mp_par));
    parameterConstraintsAllocated = true;
    for (int i = 0; i < nParamsTot; i++) {
      mpfitParameterConstraints[i].fixed = parameterLimits[i].fixed;
      mpfitParameterConstraints[i].limited[0] = parameterLimits[i].limited[0];
      mpfitParameterConstraints[i].limited[1] = parameterLimits[i].limited[1];
      mpfitParameterConstraints[i].limits[0] = parameterLimits[i].limits[0];
      mpfitParameterConstraints[i].limits[1] = parameterLimits[i].limits[1];
    }
  }
  
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  memset(&mpfitResult, 0, sizeof(mpfitResult));       /* Zero results structure */
  mpfitResult.xerror = paramErrs;
  memset(&mpConfig, 0, sizeof(mpConfig));
  mpConfig.maxiter = MAX_ITERATIONS;
  mpConfig.ftol = ftol;
  mpConfig.verbose = verbose;

  status = mpfit(myfunc_mpfit, nDataVals, nParamsTot, paramVector, mpfitParameterConstraints,
					&mpConfig, theModel, &mpfitResult);

  // Store information about the optimization, if SolverResults object was supplied
  if (solverResults != NULL) {
    solverResults->SetSolverType(MPFIT_SOLVER);
    solverResults->AddMPResults(mpfitResult);
  }

  if (parameterConstraintsAllocated)
    free(mpfitParameterConstraints);
  free(paramErrs);
  return status;
}



/* END OF FILE: levmar_fit.cpp ------------------------------------------- */
