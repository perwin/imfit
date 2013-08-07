/* FILE: levmar_fit.cpp -------------------------------------------------- */

// Copyright 2012, 2013 by Peter Erwin.
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

#include <strings.h>   // for bzero
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "print_results.h"

const int  MAX_ITERATIONS = 1000;
const double  FTOL = 1.0e-8;
const double  XTOL = 1.0e-8;


/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
int myfunc_mpfit( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );




/* This is the function used by mpfit() to compute the vector of deviates.
 * In our case, it's a wrapper which tells the ModelObject to compute the deviates.
 */
int myfunc_mpfit( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *theModel )
{

  theModel->ComputeDeviates(deviates, params);
  return 0;
}





int LevMarFit( int nParamsTot, int nFreeParams, int nDataVals, double *paramVector, 
				mp_par *parameterLimits, ModelObject *theModel, double ftol, 
				bool paramLimitsExist, int verbose )
{
  double  *paramErrs;
  mp_par  *mpfitParameterConstraints;
  mp_result  mpfitResult;
  mp_config  mpConfig;
  int  status;

  if (! paramLimitsExist) {
    // If parameters are unconstrained, then mpfit() expects a NULL mp_par array
    mpfitParameterConstraints = NULL;
  } else {
    mpfitParameterConstraints = parameterLimits;
  }

  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero results structure */
  mpfitResult.xerror = paramErrs;
  bzero(&mpConfig, sizeof(mpConfig));
  mpConfig.maxiter = MAX_ITERATIONS;
  mpConfig.ftol = ftol;
  mpConfig.verbose = verbose;
//   if (verbose)
//     mpConfig.verbose = 1;
//   else
//     mpConfig.verbose = 0;

  status = mpfit(myfunc_mpfit, nDataVals, nParamsTot, paramVector, mpfitParameterConstraints,
					&mpConfig, theModel, &mpfitResult);

  if (verbose >= 0) {
    printf("\n");
    PrintResults(paramVector, 0, &mpfitResult, theModel, nFreeParams, parameterLimits, status);
    printf("\n");
  }

  return status;
}




/* END OF FILE: levmar_fit.cpp ------------------------------------------- */
