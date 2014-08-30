/* FILE: new_levmar_fit.cpp ---------------------------------------------- */

// Copyright 2014 by Peter Erwin.
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



// * Return values from nlopt library:
// NLOPT_SUCCESS = 1
// Generic success return value.
// NLOPT_STOPVAL_REACHED = 2
// Optimization stopped because stopval (above) was reached.
// NLOPT_FTOL_REACHED = 3
// Optimization stopped because ftol_rel or ftol_abs (above) was reached.
// NLOPT_XTOL_REACHED = 4
// Optimization stopped because xtol_rel or xtol_abs (above) was reached.
// NLOPT_MAXEVAL_REACHED = 5
// Optimization stopped because maxeval (above) was reached.
// NLOPT_MAXTIME_REACHED = 6
// Optimization stopped because maxtime (above) was reached.
// [edit]
// Error codes (negative return values)
// NLOPT_FAILURE = -1
// Generic failure code.
// NLOPT_INVALID_ARGS = -2
// Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).
// NLOPT_OUT_OF_MEMORY = -3
// Ran out of memory.
// NLOPT_ROUNDOFF_LIMITED = -4
// Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)
// NLOPT_FORCED_STOP = -5
// Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints.


// Note : the following are the default tolerance values we are currently using
// in mpfitfun.cpp:
//  conf.ftol = 1e-10;   [relative changes in chi^2]
//  conf.xtol = 1e-10;   [relative changes in parameter values]

#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// following is for DBL_MAX
#include <float.h>

#include "model_object.h"
#include "model_object_1d.h"
#include "param_struct.h"   // for mp_par structure
#include "levmar-2.6/levmar.h"
#include "parameter_utils.h"
//#include "newlevmar/newlevmar.h"

const double  FTOL = 1.0e-8;
const double  XTOL = 1.0e-8;
const int  FUNCS_PER_REPORTING_STEP = 20;
const int  REPORT_STEPS_PER_VERBOSE_OUTPUT = 5;


// Module variables -- used in myfunc
static int  verboseOutput;
static int  dimensionality;   // 1D or 2D fitting?
static double *origInputParams;
//static double *fullParams;
static bool  *fixedParamFlags;
static bool  fixedParamsExist = false;
static int  nTotParams;
static int  nFreeParams;
// the following is used *only* inside myfunc_levmar, but we declare it as a
// module variable so we don't have to continually allocate and free it each
// time myfunc_levmar is called
double  *currentFullParams;



// NEW: Object function with interface for levmar code
// p = parameter vector
// e = (unsquared) deviate values = x - hx
// x = data values
// hx = model values
// n = number of data points (pixels)
// m = number of parameters

// The purpose of this function is to populate hx[] with the model values
// corresponding to the parameter vector p
// Note that in the case of one or more fixed parameters, len(p) < len(full_parameter_set)
// (since the main LM function can't handle missing parameters, so we feed it a
// condensed parameter vector with the fixed values removed), and then fill in the
// missing values using the original parameter vector (which was input to NewLevMarFit 
// and stored in the module variable origInputParams)
void myfunc_levmar( double *p, double *hx, int m, int n, void *my_func_data )
{
  double  *localPointerToModelVals;
  double  *paramsToUse;   // never allocated; used to point to other (allocated) vectors
  int  i;
  ModelObject *theModel = (ModelObject *)my_func_data;

  if (fixedParamsExist) {
    ExpandParamVector(origInputParams, p, currentFullParams, nTotParams, m, fixedParamFlags);
    paramsToUse = currentFullParams;
  }
  else
    paramsToUse = p;
  
  if (dimensionality == 1) {
    theModel->CreateModelImage(paramsToUse);
    // the following copies the model values into hx:
    int dummy = theModel->GetModelVector(hx);
  }
  else {
    // for 2D case, we would do something like this:
    // NOTE: hx is allocated & freed within the levmar function, so we can't just
    // pass a pointer to the internal ModelObject array
    theModel->CreateModelImage(paramsToUse);
    localPointerToModelVals = theModel->GetModelImageVector();
    for (i = 0; i < n; i++)
      hx[i] = localPointerToModelVals[i];
   } 
}



void InterpretResult_newlevmar( double *info )
{
  double  initialChiSquare, finalChiSquare;
  int  nIterations, terminationState, nFuncEvals, nJacobianEvals, nLinearSolves;
  string  description;
  
  initialChiSquare = info[0];
  finalChiSquare = info[1];
  nIterations = (int)info[5];
  terminationState = (int)info[6];
  nFuncEvals = (int)info[7];
  nJacobianEvals = (int)info[8];
  nLinearSolves = (int)info[9];
  
  switch (terminationState) {
    case 1:
      description = "stopped by small gradient J^T e";
      break;
    case 2:
      description = "stopped by small delta-p (change in parameter value)";
      break;
    case 3:
      description = "Maximum # iterations reached.";
      break;
    case 4:
      description = "Singular matrix. Restart from current p with increased mu.";
      break;
    case 5:
      description = "No further error reduction is possible. Restart with increased mu.";
      break;
    case 6:
      description = "Stopped by small chi^2";
      break;
    case 7:
      description = "Stopped due to invalid model values (NaN, Inf, etc.)";
      break;
  }

  printf("\nLevenberg-Marquardt solver dlevmar_bc_dif finished.\n");
  printf("%d iterations, %d function evaluations, %d Jacobian evaluations.\n", 
  		nIterations, nFuncEvals, nJacobianEvals);
  printf("   Initial fit-statistic = %g\n", initialChiSquare);
  printf("   Final fit-statistic = %g\n", finalChiSquare);
  printf("Terminated with code %d = %s\n", terminationState, description.c_str());
  
}




int NewLevMarFit( int nParamsTot, double *paramVector, mp_par *parameterLimits, 
                  ModelObject *theModel, double ftol, int verbose )
{
  double  *minParamValues;
  double  *maxParamValues;
  double  *freeParams;
  double  *lowerLimits;
  double  *upperLimits;
  double  *covarianceMatrix;
  bool  paramLimitsExist = false;
  int  result;
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];  // levmar stuff
  int  nFixed = 0;
  int  i;

  // figure out if we're doing 1D or 2D fitting
  dimensionality = theModel->Dimensionality();
  
  nTotParams = nParamsTot;
  // store the original parameter values so we can recall the fixed-param values
  origInputParams = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  for (i = 0; i < nParamsTot; i++)
    origInputParams[i] = paramVector[i];
  // Figure out if any parameters are being held fixed
  nFixed = CountFixedParams(parameterLimits, nParamsTot);
  if (nFixed > 0) {
    fixedParamsExist = true;
  }    
  nFreeParams = nParamsTot - nFixed;
  
  // Stuff mainly useful inside myfunc_levmar
  fixedParamFlags = (bool *)malloc((size_t)(nParamsTot*sizeof(bool)));
  currentFullParams = (double *)calloc( (size_t)nTotParams, sizeof(double) );

  minParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  maxParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  
  freeParams = (double *)calloc( (size_t)nFreeParams, sizeof(double) );
  int  mm = nFreeParams*nFreeParams;
  covarianceMatrix = (double *)calloc( (size_t)mm, sizeof(double) );

  // Populate full-size parameter-limits vector, if parameter limits exist
  if (parameterLimits != NULL) {
    for (i = 0; i < nParamsTot; i++) {
      // default state is to have no limits on a parameter
      minParamValues[i] = -DBL_MAX;
      maxParamValues[i] = DBL_MAX;
      fixedParamFlags[i] = false;
      // check to see if user specified a fixed value for this parameter
      if (parameterLimits[i].fixed == 1)
        fixedParamFlags[i] = true;
      else if ((parameterLimits[i].limited[0] == 1) && (parameterLimits[i].limited[1] == 1)) {
        // user specified parameter limits for this parameter
        minParamValues[i] = parameterLimits[i].limits[0];
        maxParamValues[i] = parameterLimits[i].limits[1];
        paramLimitsExist = true;
      }
    }
  }

  CondenseParamVector(paramVector, freeParams, nParamsTot, nFreeParams, fixedParamFlags);
  if (paramLimitsExist) {
    lowerLimits = (double *)calloc( (size_t)nFreeParams, sizeof(double) );
    upperLimits = (double *)calloc( (size_t)nFreeParams, sizeof(double) );
    CondenseParamLimits(minParamValues, maxParamValues, lowerLimits, 
    					upperLimits, nParamsTot, nFreeParams, fixedParamFlags);
  }
  else {
    lowerLimits = NULL;
    upperLimits = NULL;
  }


  // Specify level of verbosity and start the optimization
  verboseOutput = verbose;
  
  // setup for dlevmar_bc_dif
  opts[0] = LM_INIT_MU; 
  opts[1] = 1E-15; 
  opts[2] = 1E-15; 
  opts[3] = 1E-20;
  opts[4] = LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite 
  
  int maxIters = 1000;
  int nDataVals = theModel->GetNDataValues();
  // test setup for 1D fitting of data without errors
//  int fitType = LM_CHISQ_EQUAL_WT;
  double *dataVals = theModel->GetDataVector();
  
// modified Laurence & Chromy code -- currently doesn't work
//   result = dlevmar_mle_dif(myfunc_levmar, paramVector, dataVals, nParamsTot, nDataVals, 1000, opts, 
//   						info, NULL, NULL, (void *)theModel, fitType);  // no Jacobian

//   result = dlevmar_dif(myfunc_levmar, paramVector, dataVals, nParamsTot, nDataVals, 1000, opts, 
//   						info, NULL, NULL, (void *)theModel);  // no Jacobian
  result = dlevmar_bc_dif(myfunc_levmar, freeParams, dataVals, nFreeParams, nDataVals,
  						lowerLimits, upperLimits, NULL,
  						maxIters, opts, info, NULL, covarianceMatrix, (void *)theModel); 

// int LEVMAR_BC_DIF(
//   void (*func)(LM_REAL *p, LM_REAL *hx, int m, int n, void *adata), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
//   LM_REAL *p,         /* I/O: initial parameter estimates. On output has the estimated solution */
//   LM_REAL *x,         /* I: measurement vector. NULL implies a zero vector */
//   int m,              /* I: parameter vector dimension (i.e. #unknowns) */
//   int n,              /* I: measurement vector dimension */
//   LM_REAL *lb,        /* I: vector of lower bounds. If NULL, no lower bounds apply */
//   LM_REAL *ub,        /* I: vector of upper bounds. If NULL, no upper bounds apply */
//   LM_REAL *dscl,      /* I: diagonal scaling constants. NULL implies no scaling */
//   int itmax,          /* I: maximum number of iterations */
//   LM_REAL opts[5],    /* I: opts[0-4] = minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \delta]. Respectively the
//                        * scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2 and
//                        * the step used in difference approximation to the Jacobian. Set to NULL for defaults to be used.
//                        * If \delta<0, the Jacobian is approximated with central differences which are more accurate
//                        * (but slower!) compared to the forward differences employed by default. 
//                          epsilon2 = relative tolerance for change in parameter value
//                             (quit if (delta-p)^2 < (eps2)^2 * p^2)
//                          epsilon3 = absolute tolerance for chi^2 (quit if chi^2 < eps3)
//                        */
//   LM_REAL info[LM_INFO_SZ],
// 					           /* O: information regarding the minimization. Set to NULL if don't care
//                       * info[0]= ||e||_2 at initial p.
//                       * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
//                       * info[5]= # iterations,
//                       * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
//                       *                                 2 - stopped by small Dp
//                       *                                 3 - stopped by itmax
//                       *                                 4 - singular matrix. Restart from current p with increased mu 
//                       *                                 5 - no further error reduction is possible. Restart with increased mu
//                       *                                 6 - stopped by small ||e||_2
//                       *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
//                       * info[7]= # function evaluations
//                       * info[8]= # Jacobian evaluations
//                       * info[9]= # linear systems solved, i.e. # attempts for reducing error
//   LM_REAL *work,      * working memory at least LM_BC_DIF_WORKSZ() reals large, allocated if NULL */
//   LM_REAL *covar,     * O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */
//   void *adata)        * pointer to possibly additional data, passed uninterpreted to func.

  printf("\nMinimization info:\n");
  for (i = 0; i < LM_INFO_SZ; ++i)
    printf("%g ", info[i]);
  printf("\n");
  
  // Copy best-fitting parameters into input/output vector paramVector, filling in
  // fixed-parameter values from original input if necessary
  ExpandParamVector(origInputParams, freeParams, paramVector, nTotParams, nFreeParams, 
  					fixedParamFlags);

  if (verbose >= 0)
    InterpretResult_newlevmar(info);

  // Don't know how to interpret/convert this yet, so we'll skip showing it
//   printf("Covariance matrix:\n");
//   for (i = 0; i < nFreeParams; i++) {
//     for (int j = 0; j < nFreeParams; j++) {
//       printf("%12g", covarianceMatrix[i*nFreeParams + j]);
//     }
//     printf("\n");
//   }

  // Free arrays:
  free(origInputParams);
  free(fixedParamFlags);
  free(currentFullParams);
  free(minParamValues);
  free(maxParamValues);
  free(freeParams);
  free(covarianceMatrix);
  if (paramLimitsExist) {
    free(lowerLimits);
    free(upperLimits);
  }

  return result;
}




/* END OF FILE: new_levmar_fit.cpp --------------------------------------- */
