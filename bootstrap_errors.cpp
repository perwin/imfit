/* FILE: bootstrap_errors.cpp ------------------------------------------ */
/* VERSION 0.2
 *
 * Code for estimating errors on fitted parameters (for a 1D profile fit via
 * profilefit) via bootstrap resampling.
 *
 *     [v0.1]: 11 Jan 2013: Created; initial development.
 *
 */

// Copyright 2013 by Peter Erwin.
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "definitions.h"
#include "model_object.h"
#include "mpfit_cpp.h"
#include "levmar_fit.h"
#include "nmsimplex_fit.h"
#include "mersenne_twister.h"
#include "bootstrap_errors.h"
#include "statistics.h"
#include "print_results.h"


/* ------------------- Function Prototypes ----------------------------- */



void BootstrapErrors( double *bestfitParams, mp_par *parameterLimits, bool paramLimitsExist, 
					ModelObject *theModel, double ftol, int nIterations, int nFreeParams,
					bool usingCashStatistic )
{
  double  *paramsVect, *paramSigmas;
  double  **paramArray;
  double  lower, upper, plus, minus, halfwidth;
  int  i, status, nIter;
  int  nParams = theModel->GetNParams();
  int  nStoredDataVals = theModel->GetNDataValues();
  int  nValidPixels = theModel->GetNValidPixels();
  int  verboseLevel = -1;
  
  /* seed random number generators with current time */
  init_genrand( (unsigned long)time(NULL) );

  paramsVect = (double *) malloc(nParams * sizeof(double));
  // Allocate 2D array to hold bootstrap results for each parameter
  paramArray = (double **)calloc( (size_t)nParams, sizeof(double *) );
  for (i = 0; i < nParams; i++)
    paramArray[i] = (double *)calloc( (size_t)nIterations, sizeof(double) );
  // vector to hold estimated sigmas for each parameter
  paramSigmas = (double *)calloc( (size_t)nParams, sizeof(double) );


  theModel->UseBootstrap();

  if (! usingCashStatistic)
    printf("\nStarting bootstrap iterations (L-M solver): ");
  else
    printf("\nStarting bootstrap iterations (N-M simplex solver): ");

  for (nIter = 0; nIter < nIterations; nIter++) {
    printf("%d...  ", nIter + 1);
    fflush(stdout);
    theModel->MakeBootstrapSample();
    for (i = 0; i < nParams; i++)
      paramsVect[i] = bestfitParams[i];
    if (! usingCashStatistic) {
      status = LevMarFit(nParams, nFreeParams, nValidPixels, paramsVect, parameterLimits, 
      					theModel, ftol, paramLimitsExist, verboseLevel);
    } else {
      status = NMSimplexFit(nParams, paramsVect, parameterLimits, theModel, ftol,
      						verboseLevel);
    }
    for (i = 0; i < nParams; i++) {
      paramArray[i][nIter] = paramsVect[i];
    }
  }


  /* Determine dispersions for parameter values */
  for (i = 0; i < nParams; i++) {
    paramSigmas[i] = StandardDeviation(paramArray[i], nIterations);
  }
  
  /* Print parameter values + standard deviations: */
  /* (note that calling ConfidenceInterval() sorts the vectors in place!) */
  printf("\nStatistics for parameter values from bootstrap resampling");
  printf(" (%d rounds):\n", nIterations);
  printf("Best-fit\t\t Bootstrap      [68%% conf.int., half-width]; (mean +/- standard deviation)\n");
  for (i = 0; i < nParams; i++) {
    if ((paramLimitsExist) && (parameterLimits[i].fixed == 0)) {
      // OK, this parameter was not fixed
      ConfidenceInterval(paramArray[i], nIterations, &lower, &upper);
      plus = upper - bestfitParams[i];
      minus = bestfitParams[i] - lower;
      halfwidth = (upper - lower)/2.0;
      printf("%s = %g  +%g, -%g    [%g -- %g, %g];  (%g +/- %g)\n", 
             theModel->GetParameterName(i).c_str(), 
             bestfitParams[i], plus, minus, lower, upper, halfwidth,
             Mean(paramArray[i], nIterations), paramSigmas[i]);
    }
    else {
      printf("%s = %g     [fixed parameter]\n", theModel->GetParameterName(i).c_str(),
                  bestfitParams[i]);
    
    }
  }


  free(paramsVect);
  free(paramSigmas);
  for (i = 0; i < nParams; i++)
    free(paramArray[i]);
  free(paramArray);

}



/* END OF FILE: bootstrap_errors.cpp ----------------------------------- */
