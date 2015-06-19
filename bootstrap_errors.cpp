/* FILE: bootstrap_errors.cpp ------------------------------------------ */
/* 
 * Code for estimating errors on fitted parameters (for a 1D profile fit via
 * profilefit) via bootstrap resampling.
 *
 *     [v0.1]: 11 Jan 2013: Created; initial development.
 *
 * Note that some of this code was taken from bootstrap2.cpp, part of
 * nonlinfit (imfit's conceptual predecessor), so yay for reuse!
 */

// Copyright 2013-2015 by Peter Erwin.
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


/* ------------------------ Include Files (Header Files )--------------- */

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "definitions.h"
#include "model_object.h"
#include "mpfit_cpp.h"
#include "levmar_fit.h"
#ifndef NO_NLOPT
#include "nmsimplex_fit.h"
#endif
#include "diff_evoln_fit.h"
#include "mersenne_twister.h"
#include "bootstrap_errors.h"
#include "statistics.h"
#include "print_results.h"

using namespace std;



/* ------------------- Function Prototypes ----------------------------- */

int BootstrapErrorsBase( double *bestfitParams, mp_par *parameterLimits, bool paramLimitsExist, 
					ModelObject *theModel, double ftol, int nIterations, int nFreeParams,
					int whichStatistic, double **outputParamArray, FILE *outputFile_ptr );




/* ---------------- FUNCTION: BootstrapErrors -------------------------- */
// Primary wrapper function (meant to be called from main() in imfit_main.cpp, etc.).
// If saving of all best-fit parameters to file is requested, then outputFile_ptr
// should be non-NULL (i.e., should point to a file object opened for writing, possibly
// with header information already written).
int BootstrapErrors( double *bestfitParams, mp_par *parameterLimits, bool paramLimitsExist, 
					ModelObject *theModel, double ftol, int nIterations, int nFreeParams,
					int whichStatistic, FILE *outputFile_ptr )
{
  double  *paramSigmas;
  double **outputParamArray;
  double  lower, upper, plus, minus, halfwidth;
  int  i, nSuccessfulIterations;
  int  nParams = theModel->GetNParams();

  // Allocate 2D array to hold bootstrap results for each parameter
  outputParamArray = (double **)calloc( (size_t)nParams, sizeof(double *) );
  for (i = 0; i < nParams; i++)
    outputParamArray[i] = (double *)calloc( (size_t)nIterations, sizeof(double) );
  
  
  // write column header info to file, if user requested saving to file
  if (outputFile_ptr != NULL) {
    string  headerLine = theModel->GetParamHeader();
    fprintf(outputFile_ptr, "#\n# Bootstrap resampling output (%d iterations requested):\n%s\n", 
   			nIterations, headerLine.c_str());
  }    
  // do the bootstrap iterations (saving to file if user requested it)
  nSuccessfulIterations = BootstrapErrorsBase(bestfitParams, parameterLimits, paramLimitsExist, 
					theModel, ftol, nIterations, nFreeParams, whichStatistic, 
					outputParamArray, outputFile_ptr);
  
  
  // Calculate sigmas and 68% confidence intervals for the parameters
  // vector to hold estimated sigmas for each parameter
  paramSigmas = (double *)calloc( (size_t)nParams, sizeof(double) );
  /* Determine dispersions for parameter values */
  for (i = 0; i < nParams; i++) {
    paramSigmas[i] = StandardDeviation(outputParamArray[i], nSuccessfulIterations);
  }
  /* Print parameter values + standard deviations: */
  /* (note that calling ConfidenceInterval() sorts the vectors in place!) */
  printf("\nStatistics for parameter values from bootstrap resampling");
  printf(" (%d successful iterations):\n", nSuccessfulIterations);
  printf("Best-fit\t\t Bootstrap      [68%% conf.int., half-width]; (mean +/- standard deviation)\n");
  for (i = 0; i < nParams; i++) {
    if (parameterLimits[i].fixed == 0) {
      // OK, this parameter was not fixed
      ConfidenceInterval(outputParamArray[i], nSuccessfulIterations, &lower, &upper);
      plus = upper - bestfitParams[i];
      minus = bestfitParams[i] - lower;
      halfwidth = (upper - lower)/2.0;
      printf("%s = %g  +%g, -%g    [%g -- %g, %g];  (%g +/- %g)\n", 
             theModel->GetParameterName(i).c_str(), 
             bestfitParams[i], plus, minus, lower, upper, halfwidth,
             Mean(outputParamArray[i], nSuccessfulIterations), paramSigmas[i]);
    }
    else {
      printf("%s = %g     [fixed parameter]\n", theModel->GetParameterName(i).c_str(),
                  bestfitParams[i]);
    }
  }

  for (i = 0; i < nParams; i++)
    free(outputParamArray[i]);
  free(outputParamArray);
  free(paramSigmas);

  return nSuccessfulIterations;
}



/* ---------------- FUNCTION: BootstrapErrorsArrayOnly ----------------- */
// Alternate wrapper function: returns array of best-fit parameters in outputParamArray;
// doesn't print any summary statistics (e.g., sigmas, confidence intervals). 
// Note that outputParamArray will be allocated here; it should be de-allocated by 
// whatever function is calling this function.
int BootstrapErrorsArrayOnly( double *bestfitParams, mp_par *parameterLimits, bool paramLimitsExist, 
					ModelObject *theModel, double ftol, int nIterations, int nFreeParams,
					int whichStatistic, double **outputParamArray )
{
  int  i, nSuccessfulIterations;
  int  nParams = theModel->GetNParams();

  // Allocate 2D array to hold bootstrap results for each parameter
  outputParamArray = (double **)calloc( (size_t)nParams, sizeof(double *) );
  for (i = 0; i < nParams; i++)
    outputParamArray[i] = (double *)calloc( (size_t)nIterations, sizeof(double) );
  
  // do the bootstrap iterations
  nSuccessfulIterations = BootstrapErrorsBase(bestfitParams, parameterLimits, paramLimitsExist, 
					theModel, ftol, nIterations, nFreeParams, whichStatistic, 
					outputParamArray, NULL);
  
  return nSuccessfulIterations;
}



/* ---------------- FUNCTION: BootstrapErrorsBase ---------------------- */
// Base function called by the wrapper functions (above), which does the main work
// of overseeing the bootstrap resampling.
// Saving individual best-fit vales to file is done *if* outputFile_ptr != NULL.
int BootstrapErrorsBase( double *bestfitParams, mp_par *parameterLimits, bool paramLimitsExist, 
					ModelObject *theModel, double ftol, int nIterations, int nFreeParams,
					int whichStatistic, double **outputParamArray, FILE *outputFile_ptr )
{
  double  *paramsVect;
  int  i, status, nIter, nSuccessfulIters;
  int  nParams = theModel->GetNParams();
  int  nValidPixels = theModel->GetNValidPixels();
  int  verboseLevel = -1;   // ensure minimizer stays silent
  bool  saveToFile = false;
  
  if (outputFile_ptr != NULL)
    saveToFile = true;
  
  /* seed random number generators with current time */
  init_genrand((unsigned long)time((time_t *)NULL));

  paramsVect = (double *) malloc(nParams * sizeof(double));

  theModel->UseBootstrap();

  if ((whichStatistic == FITSTAT_CHISQUARE) || (whichStatistic == FITSTAT_POISSON_MLR))
    printf("\nStarting bootstrap iterations (L-M solver): ");
  else
#ifndef NO_NLOPT
    printf("\nStarting bootstrap iterations (N-M simplex solver): ");
#else
    printf("\nStarting bootstrap iterations (DE solver): ");
#endif

  // Bootstrap iterations:
  nSuccessfulIters = 0;
  for (nIter = 0; nIter < nIterations; nIter++) {
    printf("%d...  ", nIter + 1);
    fflush(stdout);
    theModel->MakeBootstrapSample();
    for (i = 0; i < nParams; i++)
      paramsVect[i] = bestfitParams[i];
    if ((whichStatistic == FITSTAT_CHISQUARE) || (whichStatistic == FITSTAT_POISSON_MLR)) {
      status = LevMarFit(nParams, nFreeParams, nValidPixels, paramsVect, parameterLimits, 
      					theModel, ftol, paramLimitsExist, verboseLevel);
    } else {
#ifndef NO_NLOPT
      status = NMSimplexFit(nParams, paramsVect, parameterLimits, theModel, ftol,
      						verboseLevel);
#else
      status = DiffEvolnFit(nParams, paramsVect, parameterLimits, theModel, ftol,
      						verboseLevel);
#endif
    }
    // Store parameters in array (and optionally write them to file) if fit was successful
    if (status > 0) {
      for (i = 0; i < nParams; i++)
        outputParamArray[i][nSuccessfulIters] = paramsVect[i];
      if (saveToFile) {
        for (int i = 0; i < (nParams - 1); i++)
          fprintf(outputFile_ptr, "%f\t\t", paramsVect[i]);
        fprintf(outputFile_ptr, "%f\n", paramsVect[nParams - 1]);
      }
      nSuccessfulIters += 1;
    }
  }

 
  free(paramsVect);

  return nSuccessfulIters;
}



/* END OF FILE: bootstrap_errors.cpp ----------------------------------- */
