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

// Copyright 2013-2019,2024 by Peter Erwin.
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
#include <tuple>

#include "definitions.h"
#include "definitions_multimage.h"
#include "model_object.h"
#include "levmar_fit.h"
#ifndef NO_NLOPT
#include "nmsimplex_fit.h"
#endif
#include "diff_evoln_fit.h"
#include "mersenne_twister.h"
#include "bootstrap_errors.h"
#include "statistics.h"
#include "print_results.h"
#include "utilities_pub.h"

using namespace std;


const int MIN_ITERATIONS_FOR_STATISTICS = 3;
const int PROGRESS_BAR_WIDTH = 80;

const vector<string> imageParamLabels = {"PIXEL_SCALE", "ROTATION", "FLUX_SCALE",
										"X0", "Y0"};


/* ------------------- Function Prototypes ----------------------------- */

int BootstrapErrorsBase( const double *bestfitParams, vector<mp_par> parameterLimits, 
					const bool paramLimitsExist, ModelObject *theModel, const double ftol, 
					const int nIterations, const int nFreeParams, const int whichStatistic, 
					double **outputParamArray, FILE *outputFile_ptr, unsigned long rngSeed=0 );
void PrintBootstrapSummary( const double *bestfitParams, double **outputParamArray, int nParams, 
				int nIterations, vector<mp_par> parameterLimits, ModelObject *theModel );
void PrintBootstrapSummaryMultimfit( const double *bestfitParams, double **outputParamArray, 
				int nParams, int nIterations, vector<mp_par> parameterLimits, 
				ModelObject *theModel );




/* ---------------- FUNCTION: BootstrapErrors -------------------------- */
/// Primary wrapper function (meant to be called from main() in imfit_main.cpp, etc.).
/// If saving of all best-fit parameters to file is requested, then outputFile_ptr
/// should be non-NULL (i.e., should point to a file object opened for writing, possibly
/// with header information already written).
/// Returns the number of (successful) bootstrap iterations, or returns -1 if error
/// encountered.
int BootstrapErrors( const double *bestfitParams, vector<mp_par> parameterLimits, 
					const bool paramLimitsExist, ModelObject *theModel, const double ftol, 
					const int nIterations, const int nFreeParams, const int whichStatistic, 
					FILE *outputFile_ptr, unsigned long rngSeed, bool multimfitMode )
{
  double  *paramSigmas;
  double **outputParamArray;
//   double  lower, upper, plus, minus, halfwidth;
  int  i, nSuccessfulIterations;
  int  nParams = theModel->GetNParams();

  // Allocate 2D array to hold bootstrap results for each parameter.
  // Because it's convenient to be able to pass a single *column* of values
  // (i.e., values of a single parameter from all the iterations) to various
  // statistical functions, we construct this array as [nParams][nIterations];
  // thus, we can pass all the values for a particular parameter as
  // outputParamArray[i_param]
  outputParamArray = (double **)calloc( (size_t)nParams, sizeof(double *) );
  for (i = 0; i < nParams; i++)
    outputParamArray[i] = (double *)calloc( (size_t)nIterations, sizeof(double) );

//   paramOffsets = (double *) calloc(nParams, sizeof(double));
//   bestfitParams_offsetCorrected = (double *) calloc(nParams, sizeof(double));
  
  // write column header info to file, if user requested saving to file
  if (outputFile_ptr != nullptr) {
    string  headerLine = theModel->GetParamHeader();
    fprintf(outputFile_ptr, "#\n# Bootstrap resampling output (%d iterations requested):\n%s\n", 
   			nIterations, headerLine.c_str());
  }    
  // do the bootstrap iterations (saving to file if user requested it)
  nSuccessfulIterations = BootstrapErrorsBase(bestfitParams, parameterLimits, paramLimitsExist, 
					theModel, ftol, nIterations, nFreeParams, whichStatistic, 
					outputParamArray, outputFile_ptr, rngSeed);
  
  if (nSuccessfulIterations < MIN_ITERATIONS_FOR_STATISTICS) {
    printf("\nNot enough successful bootstrap iterations (%d) for meaningful statistics!\n",
    		nSuccessfulIterations);
  }
  else {
    // Apply image-offset corrections to best-fit parameters, so their values
    // get printed corrected
//     theModel->GetImageOffsets(paramOffsets);
//     for (i = 0; i < nParams; i++)
//       bestfitParams_offsetCorrected[i] = bestfitParams[i] + paramOffsets[i];
    
    // Calculate sigmas and 68% confidence intervals for the parameters
    // vector to hold estimated sigmas for each parameter
    if (multimfitMode)
      PrintBootstrapSummaryMultimfit(bestfitParams, outputParamArray, nParams, nIterations, 
    								parameterLimits, theModel);
    else
      PrintBootstrapSummary(bestfitParams, outputParamArray, nParams, nIterations, 
    						parameterLimits, theModel);
  }

  for (i = 0; i < nParams; i++)
    free(outputParamArray[i]);
  free(outputParamArray);
//   free(paramOffsets);
//   free(bestfitParams_offsetCorrected);

  return nSuccessfulIterations;
}



/* ---------------- FUNCTION: BootstrapErrorsBase ---------------------- */
/// Base function called by the wrapper functions (above), which does the main work
/// of overseeing the bootstrap resampling.
/// Saving individual best-fit vales to file is done *if* outputFile_ptr != NULL.
/// Returns the number of successful iterations performed (-1 if an error was
/// encountered)
int BootstrapErrorsBase( const double *bestfitParams, vector<mp_par> parameterLimits, 
					const bool paramLimitsExist, ModelObject *theModel, const double ftol, 
					const int nIterations, const int nFreeParams, const int whichStatistic, 
					double **outputParamArray, FILE *outputFile_ptr, unsigned long rngSeed )
{
  double  *paramsVect, *paramOffsets;
  int  i, status, nIter, nDone, nSuccessfulIters;
  int  nParams = theModel->GetNParams();
  int  nValidPixels = theModel->GetNValidPixels();
  int  verboseLevel = -1;   // ensure minimizer stays silent
  bool  saveToFile = false;
  string  outputLine, iterTemplate;

  if (outputFile_ptr != nullptr)
    saveToFile = true;
  
  if (rngSeed > 0)
    init_genrand(rngSeed);
  else
    init_genrand((unsigned long)time((time_t *)NULL));

  paramsVect = (double *) calloc(nParams, sizeof(double));
  paramOffsets = (double *) calloc(nParams, sizeof(double));

  status = theModel->UseBootstrap();
  if (status < 0) {
    fprintf(stderr, "Error encountered during bootstrap setup!\n");
    free(paramsVect);
    return -1;
  }

  if ((whichStatistic == FITSTAT_CHISQUARE) || (whichStatistic == FITSTAT_POISSON_MLR))
    printf("Starting bootstrap iterations (L-M solver):\n");
  else
#ifndef NO_NLOPT
    printf("Starting bootstrap iterations (N-M simplex solver):\n");
#else
    printf("Starting bootstrap iterations (DE solver):\n");
#endif


  int  nDigits = floor(log10(nIterations)) + 1;
  iterTemplate = PrintToString("] %%%dd", nDigits) + " (%3.1f%%)\r";

  // Bootstrap iterations:
  nSuccessfulIters = 0;
  for (nIter = 0; nIter < nIterations; nIter++) {
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
    // Note that paramsVect has subsection-relative values of X0,Y0, so we need to
    // correct them with paramOffsets
    theModel->GetImageOffsets(paramOffsets);
    if (status > 0) {
      for (i = 0; i < nParams; i++) {
        outputParamArray[i][nSuccessfulIters] = paramsVect[i] + paramOffsets[i];
      }
      if (saveToFile) {
        // use paramsVect because PrintModelParamsHorizontalString will automatically
        // apply image-offset corrections
        outputLine = theModel->PrintModelParamsHorizontalString(paramsVect);
        fprintf(outputFile_ptr, "%s\n", outputLine.c_str());
      }
      nSuccessfulIters += 1;
    }
    
    // print/update progress bar
    nDone = nIter + 1;
    PrintProgressBar(nDone, nIterations, iterTemplate, PROGRESS_BAR_WIDTH);
  }
  printf("\n");

 
  free(paramsVect);
  free(paramOffsets);

  return nSuccessfulIters;
}



/* ---------------- FUNCTION: BootstrapErrorsArrayOnly ----------------- */
/// This is the same as BootstrapErrorsBase *except* that outputParamArray is
/// a 1-D array, for ease of transfer to and from Numpy arrays in the Cython
/// wrapper code in PyImfit, *and* that there is no saving to file; printing
/// a progress bar is one only if verboseFlag = true.
int BootstrapErrorsArrayOnly( const double *bestfitParams, vector<mp_par> parameterLimits, 
					const bool paramLimitsExist, ModelObject *theModel, const double ftol, 
					const int nIterations, const int nFreeParams, const int whichStatistic, 
					double *outputParamArray, unsigned long rngSeed, bool verboseFlag )
{
  double  *paramsVect, *paramOffsets;
  int  i, status, nIter, nDone, nSuccessfulIters;
  int  nParams = theModel->GetNParams();
  int  nValidPixels = theModel->GetNValidPixels();
  int  verboseLevel = -1;   // ensure minimizer stays silent
  string  iterTemplate;

  if (rngSeed > 0)
    init_genrand(rngSeed);
  else
    init_genrand((unsigned long)time((time_t *)NULL));

  paramsVect = (double *) calloc(nParams, sizeof(double));
  paramOffsets = (double *) calloc(nParams, sizeof(double));

  status = theModel->UseBootstrap();
  if (status < 0) {
    fprintf(stderr, "Error encountered during bootstrap setup!\n");
    free(paramsVect);
    return -1;
  }

  int  nDigits = floor(log10(nIterations)) + 1;
  iterTemplate = PrintToString("] %%%dd", nDigits) + " (%3.1f%%)\r";

  // Bootstrap iterations:
  if (verboseFlag)
    printf("Starting %d rounds of bootstrap resampling:\n", nIterations);
  nSuccessfulIters = 0;
  for (nIter = 0; nIter < nIterations; nIter++) {
  	if (verboseFlag)
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
    // Store parameters in array if fit was successful.
    // Note that paramsVect has image-subsection-relative values of X0,Y0, 
    // so we need to correct them with paramOffsets
    theModel->GetImageOffsets(paramOffsets);
    if (status > 0) {
      for (int j = 0; j < nParams; j++) {   // j = column number
        outputParamArray[nSuccessfulIters*nParams + j] = paramsVect[j] + paramOffsets[j];
      }
      nSuccessfulIters += 1;
    }

	if (verboseFlag) {
      // print/update progress bar
      nDone = nIter + 1;
      PrintProgressBar(nDone, nIterations, iterTemplate, PROGRESS_BAR_WIDTH);
    }
  }
 
  if (verboseFlag)
    printf("\n");

  free(paramsVect);
  free(paramOffsets);

  return nSuccessfulIters;
}


void PrintBootstrapSummary( const double *bestfitParams, double **outputParamArray, int nParams, 
				int nIterations, vector<mp_par> parameterLimits, ModelObject *theModel )
{
  double  lower, upper, plus, minus, halfwidth;
  double  *paramSigmas;
  int  i;
  
  paramSigmas = (double *)calloc( (size_t)nParams, sizeof(double) );
  for (i = 0; i < nParams; i++)
    paramSigmas[i] = StandardDeviation(outputParamArray[i], nIterations);
  // Print parameter values + standard deviations, for non-fixed parameters
  // (note that calling ConfidenceInterval() sorts the vectors in place!)
  printf("\nStatistics for parameter values from bootstrap resampling");
  printf(" (%d successful iterations):\n", nIterations);
  printf("Best-fit\t\t Bootstrap      [68%% conf.int., half-width]; (mean +/- standard deviation)\n");
  for (i = 0; i < nParams; i++) {
    if (parameterLimits[i].fixed == 0) {
      std::tie(lower, upper) = ConfidenceInterval(outputParamArray[i], nIterations);
      plus = upper - bestfitParams[i];
      minus = bestfitParams[i] - lower;
      halfwidth = (upper - lower)/2.0;
      printf("%s = %g  +%g, -%g    [%g -- %g, %g];  (%g +/- %g)\n", 
             theModel->GetParameterName(i).c_str(), 
             bestfitParams[i], plus, minus, lower, upper, halfwidth,
             Mean(outputParamArray[i], nIterations), paramSigmas[i]);
    }
    else {
      printf("%s = %g     [fixed parameter]\n", theModel->GetParameterName(i).c_str(),
                  bestfitParams[i]);
    }
  }
  free(paramSigmas);
}



void PrintBootstrapSummaryMultimfit( const double *bestfitParams, double **outputParamArray, 
				int nParams, int nIterations, vector<mp_par> parameterLimits, 
				ModelObject *theModel )
{
  double  lower, upper, plus, minus, halfwidth;
  double  *paramSigmas;
  int  i, i_last, nImages;
  
  nImages = theModel->GetNImages();
  
  paramSigmas = (double *)calloc( (size_t)nParams, sizeof(double) );
  for (i = 0; i < nParams; i++)
    paramSigmas[i] = StandardDeviation(outputParamArray[i], nIterations);
  // Print parameter values + standard deviations, for non-fixed parameters
  // (note that calling ConfidenceInterval() sorts the vectors in place!)
  printf("\nStatistics for parameter values from bootstrap resampling");
  printf(" (%d successful iterations):\n", nIterations);
  printf("Best-fit\t\t Bootstrap      [68%% conf.int., half-width]; (mean +/- standard deviation)\n");
  
  i = -1;
  // print results for image-description parameters (for 2nd and subsequent images)
  for (int n = 1; n < nImages; n++) {
    printf("# Image-description parameters for image %d:\n", n + 1);
    for (int j = 0; j < N_IMAGE_PARAMS; j++) {
      i += 1;
      if (parameterLimits[i].fixed == 0) {
        std::tie(lower, upper) = ConfidenceInterval(outputParamArray[i], nIterations);
        plus = upper - bestfitParams[i];
        minus = bestfitParams[i] - lower;
        halfwidth = (upper - lower)/2.0;
        printf("%s = %g  +%g, -%g    [%g -- %g, %g];  (%g +/- %g)\n", imageParamLabels[i].c_str(),
      	       	bestfitParams[i], plus, minus, lower, upper, halfwidth,
        		Mean(outputParamArray[i], nIterations), paramSigmas[i]);
      }
      else
        printf("%s = %g     [fixed parameter]\n", imageParamLabels[i].c_str(),
        		bestfitParams[i]);
    }
  }
  
  // now print results for model parameters
  i_last = i;
  printf("# Model parameters:\n");
  for (i = i_last; i < nParams; i++) {
    if (parameterLimits[i].fixed == 0) {
      std::tie(lower, upper) = ConfidenceInterval(outputParamArray[i], nIterations);
      plus = upper - bestfitParams[i];
      minus = bestfitParams[i] - lower;
      halfwidth = (upper - lower)/2.0;
      printf("%s = %g  +%g, -%g    [%g -- %g, %g];  (%g +/- %g)\n", 
      	 		theModel->GetParameterName(i).c_str(), bestfitParams[i], plus, minus, 
        		lower, upper, halfwidth, Mean(outputParamArray[i], nIterations), paramSigmas[i]);
    }
    else
      printf("%s = %g     [fixed parameter]\n", theModel->GetParameterName(i).c_str(),
               bestfitParams[i]);
  }

  free(paramSigmas);
}


/* END OF FILE: bootstrap_errors.cpp ----------------------------------- */
