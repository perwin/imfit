
// Copyright 2015 by Peter Erwin.
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


#include <stdlib.h>

#include "definitions.h"
#include "mpfit_cpp.h"
#include "solver_results.h"

using namespace std;


SolverResults::SolverResults( )
{
  whichSolver = MPFIT_SOLVER;
  whichFitStatistic = FITSTAT_CHISQUARE;
  bestFitValue = 0.0;
  paramSigmasAllocated = false;
  paramSigmas = NULL;
  mpResultsPresent = false;
  
  nFuncEvals = 0;
  
  mpResult.bestnorm = 0.0;
  mpResult.orignorm = 0.0;
  mpResult.niter = 0;
  mpResult.nfev = 0;
  mpResult.status = 0;
  mpResult.npar = 0;
  mpResult.nfree = 0;
  mpResult.npegged = 0;
  mpResult.nfunc = 0;
  mpResult.resid = NULL;
  mpResult.xerror = NULL;
  mpResult.covar = NULL;
  
}


SolverResults::~SolverResults( )
{
  if (paramSigmasAllocated)
    free(paramSigmas);
}


// The only place where an mp_result structure is used is as a local variable
// in LevMarFit(). Thus, we need to *copy* that structure into this object
// if we want to access its values elsewhere, since the original will go out
// of scope and vanish when LevMarFit() returns.
// 
// struct mp_result_struct {
//   double bestnorm;     /* Final chi^2 */
//   double orignorm;     /* Starting value of chi^2 */
//   int niter;           /* Number of iterations */
//   int nfev;            /* Number of function evaluations */
//   int status;          /* Fitting status code */
//   
//   int npar;            /* Total number of parameters */
//   int nfree;           /* Number of free parameters */
//   int npegged;         /* Number of pegged parameters */
//   int nfunc;           /* Number of residuals (= num. of data points) */
// 
//   double *resid;       /* Final residuals
// 			  nfunc-vector, or 0 if not desired */
//   double *xerror;      /* Final parameter uncertainties (1-sigma)
// 			  npar-vector, or 0 if not desired */
//   double *covar;       /* Final parameter covariance matrix
// 			  npar x npar array, or 0 if not desired */

void SolverResults::AddMPResults( mp_result& mpResult_input )
{
  mpResult.bestnorm = mpResult_input.bestnorm;
  bestFitValue = mpResult_input.bestnorm;
  mpResult.orignorm = mpResult_input.orignorm;
  mpResult.niter = mpResult_input.niter;
  mpResult.nfev = mpResult_input.nfev;
  nFuncEvals = mpResult_input.nfev;
  mpResult.status = mpResult_input.status;
  mpResult.npar = mpResult_input.npar;
  nParameters = mpResult_input.npar;
  mpResult.nfree = mpResult_input.nfree;
  mpResult.npegged = mpResult_input.npegged;
  mpResult.nfunc = mpResult_input.nfunc;
  // if parameter uncertainties are present, copy those into the paramSigmas array;
  // don't both keeping an extra copy within the mp_result struct
  if (mpResult_input.xerror != NULL) {
    paramSigmas = (double *)calloc(nParameters, sizeof(double));
    paramSigmasAllocated = true;
    for (int i = 0; i < nParameters; i++)
      paramSigmas[i] = mpResult_input.xerror[i];
  }
  // ignore mpResult_input.resid and mpResult_input.covar
  
  mpResultsPresent = true;
}


void SolverResults::SetSolverType( int solverType )
{
  whichSolver = solverType;
}

int SolverResults::GetSolverType( )
{
  return whichSolver;
}


void SolverResults::SetFitStatisticType( int fitStatType )
{
  whichFitStatistic = fitStatType;
}

int SolverResults::GetFitStatisticType( )
{
  return whichFitStatistic;
}


void SolverResults::SetBestfitStatisticValue( double fitStatValue )
{
  bestFitValue = fitStatValue;
}

double SolverResults::GetBestfitStatisticValue( )
{
  return bestFitValue;
}



void SolverResults::SetNFunctionEvals( int nFunctionEvals )
{
  nFuncEvals = nFunctionEvals;
}

int SolverResults::GetNFunctionEvals( )
{
  return nFuncEvals;
}


void SolverResults::SetErrors( double *errors, int nParams )
{
  nParameters = nParams;
  paramSigmas = (double *)calloc(nParameters, sizeof(double));
  paramSigmasAllocated = true;
  for (int i = 0; i < nParameters; i++)
    paramSigmas[i] = errors[i];
}
