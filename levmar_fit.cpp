/* FILE: levmar_fit.cpp -------------------------------------------------- */


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





int LevMarFit( int nParamsTot, int nFreeParams, int nPixelsTot, double *paramVector, 
				mp_par *parameterLimits, ModelObject *theModel, double ftol, 
				bool paramLimitsExist, bool verbose )
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
  if (verbose)
    mpConfig.verbose = 1;
  else
    mpConfig.verbose = 0;

  status = mpfit(myfunc_mpfit, nPixelsTot, nParamsTot, paramVector, mpfitParameterConstraints,
					&mpConfig, theModel, &mpfitResult);

  printf("\n");
  PrintResults(paramVector, 0, &mpfitResult, theModel, nFreeParams, parameterLimits, status);
  printf("\n");

  return status;
}




/* END OF FILE: levmar_fit.cpp ------------------------------------------- */
