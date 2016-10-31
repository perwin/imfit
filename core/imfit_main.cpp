/* FILE: imfit_main.cpp -------------------------------------------------- */
/*
 * This is the main program file for imfit.
 *
 * Useful reminder about FITS image sizes -- the proper translations are:
 * NAXIS1 = naxes[0] = nColumns = sizeX;
 * NAXIS2 = naxes[1] = nRows = sizeY.
 *
 *
 * HISTORY
 *    10 Nov--2 Dec 2009: Early stages of development
*/

// Copyright 2009--2016 by Peter Erwin.
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
#include <strings.h>
#include <string>
#include <vector>
#include <sys/time.h>
#include "fftw3.h"

#include "definitions.h"
#include "utilities_pub.h"
#include "image_io.h"
#include "model_object.h"
#include "add_functions.h"
#include "param_struct.h"   // for mp_par structure
#include "solver_results.h"
#include "bootstrap_errors.h"
#include "option_struct_imfit.h"

// Solvers (optimization algorithms)
#include "dispatch_solver.h"
#include "levmar_fit.h"
#include "diff_evoln_fit.h"
#ifndef NO_NLOPT
#include "nmsimplex_fit.h"
#include "nlopt_fit.h"
#endif

#include "commandline_parser.h"
#include "config_file_parser.h"
#include "print_results.h"
#include "estimate_memory.h"
#include "sample_configs.h"


/* ---------------- Definitions & Constants ----------------------------- */

// Option names for use in config files
static string  kGainString = "GAIN";
static string  kReadNoiseString = "READNOISE";
static string  kExpTimeString = "EXPTIME";
static string  kNCombinedString = "NCOMBINED";
static string  kOriginalSkyString = "ORIGINAL_SKY";


#ifdef USE_OPENMP
#define VERSION_STRING      "1.4.0 (OpenMP-enabled)"
#else
#define VERSION_STRING      "1.4.0"
#endif



/* ------------------- Function Prototypes ----------------------------- */
void DetermineImageOffset( const std::string &fullImageName, double *x_offset,
					double *y_offset);
void ProcessInput( int argc, char *argv[], imfitCommandOptions *theOptions );
bool RequestedFilesPresent( imfitCommandOptions &theOptions );
void HandleConfigFileOptions( configOptions *configFileOptions, 
								imfitCommandOptions *mainOptions );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */





/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  int  nColumns, nRows;
  long  nPixels_tot;
  int  nRows_psf = 0;
  int  nColumns_psf = 0;
  long  nPixels_psf = 0;
  int  nErrColumns, nErrRows, nMaskColumns, nMaskRows;
  long  nDegFreedom;
  int  nParamsTot, nFreeParams;
  double  *allPixels;
  double  *psfPixels;
  double  *allErrorPixels;
  bool  errorPixels_allocated = false;
  double  *allMaskPixels;
  bool  maskAllocated = false;
  double  *psfOversampledPixels;
  int  nColumns_psf_oversampled, nRows_psf_oversampled;
  long  nPixels_psf_oversampled;
  double  *paramsVect;
  double  X0_offset = 0.0;
  double  Y0_offset = 0.0;
  int  x1_oversample = 0;
  int  x2_oversample = 0;
  int  y1_oversample = 0;
  int  y2_oversample = 0;
  std::string  noiseImage;
  ModelObject  *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  FunctionBlockIndices;
  bool  paramLimitsExist = false;
  bool  parameterInfo_allocated = false;
  mp_par  *parameterInfo;
  int  status, fitStatus, nSucessfulIterations;
  SolverResults  resultsFromSolver;
  vector<string>  imageCommentsList;
  imfitCommandOptions  options;
  configOptions  userConfigOptions;
  const std::string  X0_string("X0");
  const std::string  Y0_string("Y0");
  string  progNameVersion = "imfit ";
  vector<string> programHeader;
  FILE  *bootstrapSaveFile_ptr = NULL;
  bool  didBootstrap = false;
  // timing-related
  struct timeval  timer_start_all, timer_end_all;
  struct timeval  timer_start_fit, timer_end_fit;
  struct timeval  timer_start_bootstrap, timer_end_bootstrap;


  gettimeofday(&timer_start_all, NULL);

  progNameVersion += VERSION_STRING;
  MakeOutputHeader(&programHeader, progNameVersion, argc, argv);

 
  /* Define default options, then process the command line */
  SetDefaultImfitOptions(&options);
  ProcessInput(argc, argv, &options);

  // Check for presence of user-requested files; if any are missing, quit.
  // (Appropriate error messages regarding which files are missing will be printed
  // to stderr by RequestedFilesPresent)
  if (! RequestedFilesPresent(options)) {
    fprintf(stderr, "\n");
    exit(-1);
  }


  /* Read configuration file, parse & process user-supplied (non-function-related) values */
  status = ReadConfigFile(options.configFileName, true, functionList, parameterList, 
  							paramLimits, FunctionBlockIndices, paramLimitsExist, userConfigOptions);
  if (status != 0) {
    fprintf(stderr, "\n*** ERROR: Failure reading configuration file!\n\n");
    return -1;
  }
  HandleConfigFileOptions(&userConfigOptions, &options);

  
  if (options.noImage) {
    fprintf(stderr, "*** ERROR: No image to fit!\n\n");
    return -1;
  }

  /* Get image data and sizes */
  printf("Reading data image (\"%s\") ...\n", options.imageFileName.c_str());
  allPixels = ReadImageAsVector(options.imageFileName, &nColumns, &nRows);
  if (allPixels == NULL) {
    fprintf(stderr,  "\n*** ERROR: Unable to read image file \"%s\"!\n\n", 
    			options.imageFileName.c_str());
    exit(-1);
  }
  // Reminder: nColumns = n_pixels_per_row = x-size; nRows = n_pixels_per_column = y-size
  nPixels_tot = (long)nColumns * (long)nRows;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns, nRows, nPixels_tot);
  // Determine X0,Y0 pixel offset values if user specified an image section
  DetermineImageOffset(options.imageFileName, &X0_offset, &Y0_offset);

  /* Get and check mask image */
  if (options.maskImagePresent) {
    printf("Reading mask image (\"%s\") ...\n", options.maskFileName.c_str());
    allMaskPixels = ReadImageAsVector(options.maskFileName, &nMaskColumns, &nMaskRows);
    if (allMaskPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read mask file \"%s\"!\n\n", 
    			options.maskFileName.c_str());
      exit(-1);
    }
    if ((nMaskColumns != nColumns) || (nMaskRows != nRows)) {
      fprintf(stderr, "\n*** ERROR: Dimensions of mask image (%s: %d columns, %d rows)\n",
             options.maskFileName.c_str(), nMaskColumns, nMaskRows);
      fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      exit(-1);
    }
    maskAllocated = true;
  }
           
  /* Get and check error image, if supplied */
  if (options.noiseImagePresent) {
    printf("Reading noise image (\"%s\") ...\n", options.noiseFileName.c_str());
    allErrorPixels = ReadImageAsVector(options.noiseFileName, &nErrColumns, &nErrRows);
    if (allErrorPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read noise-image file \"%s\"!\n\n", 
    			options.noiseFileName.c_str());
      exit(-1);
    }
    errorPixels_allocated = true;
    if ((nErrColumns != nColumns) || (nErrRows != nRows)) {
      fprintf(stderr, "\n*** ERROR: Dimensions of error image (%s: %d columns, %d rows)\n",
             noiseImage.c_str(), nErrColumns, nErrRows);
      fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      exit(-1);
    }
  }
  
  /* Read in PSF image, if supplied */
  if (options.psfImagePresent) {
    printf("Reading PSF image (\"%s\") ...\n", options.psfFileName.c_str());
    psfPixels = ReadImageAsVector(options.psfFileName, &nColumns_psf, &nRows_psf);
    if (psfPixels == NULL) {
      fprintf(stderr,  "\n*** ERROR: Unable to read PSF image file \"%s\"!\n\n", 
    			options.psfFileName.c_str());
      exit(-1);
    }
    nPixels_psf = (long)nColumns_psf * (long)nRows_psf;
    printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns_psf, nRows_psf, nPixels_psf);
  }
  else
    printf("* No PSF image supplied -- no image convolution will be done!\n");

  /* Read in oversampled PSF image, if supplied */
  if (options.psfOversampledImagePresent) {
    if (options.psfOversamplingScale < 1) {
      fprintf(stderr, "\n*** ERROR: the oversampling scale for the oversampled PSF was not supplied!\n");
      fprintf(stderr, "           (use --overpsf_scale to specify the scale)\n\n");
      exit(-1);
    }
    if (! options.oversampleRegionSet) {
      fprintf(stderr, "\n*** ERROR: the oversampling region was not defined!\n");
      fprintf(stderr, "           (use --overpsf_region to specify the region)\n\n");
      exit(-1);
    }
    printf("Reading oversampled PSF image (\"%s\") ...\n", options.psfOversampledFileName.c_str());
    psfOversampledPixels = ReadImageAsVector(options.psfOversampledFileName, 
    							&nColumns_psf_oversampled, &nRows_psf_oversampled);
    if (psfOversampledPixels == NULL) {
      fprintf(stderr, "\n*** ERROR: Unable to read oversampled PSF image file \"%s\"!\n\n", 
    			options.psfOversampledFileName.c_str());
      exit(-1);
    }
    nPixels_psf_oversampled = (long)nColumns_psf_oversampled * (long)nRows_psf_oversampled;
    printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns_psf_oversampled, nRows_psf_oversampled, nPixels_psf_oversampled);
    // Determine oversampling region; correct for X0,Y0 pixel offset values if user 
    // specified an image section (if not, X0_offset and Y0_offset will be = 0)
    GetAllCoordsFromBracket(options.psfOversampleRegion, &x1_oversample, &x2_oversample, 
    						&y1_oversample, &y2_oversample);
    x1_oversample -= (int)X0_offset;
    x2_oversample -= (int)X0_offset;
    y1_oversample -= (int)Y0_offset;
    y2_oversample -= (int)Y0_offset;
    string  msg, msg0;
    if ((X0_offset != 0.0) || (Y0_offset != 0.0)) {
      msg0 = "\tRegion to be oversampled within full image:";
      msg = "\tRegion to be oversampled within fitted subset image:";
      printf("%s [%d:%d,%d:%d]\n", msg0.c_str(),
      		x1_oversample + (int)X0_offset,x2_oversample + (int)X0_offset,
      		y1_oversample + (int)Y0_offset,y2_oversample + (int)Y0_offset);
    }
    else
      msg = "\tRegion to be oversampled within image:";
    printf("%s [%d:%d,%d:%d]\n", msg.c_str(), x1_oversample,x2_oversample,
    		y1_oversample,y2_oversample);
  }

  if (! options.subsamplingFlag)
    printf("* Pixel subsampling has been turned OFF.\n");


  /* Create the model object */
  theModel = new ModelObject();
  
  // Put limits on number of FFTW and OpenMP threads, if requested by user
  if (options.maxThreadsSet)
    theModel->SetMaxThreads(options.maxThreads);

  /* Add functions to the model object */
  status = AddFunctions(theModel, functionList, FunctionBlockIndices, 
  						options.subsamplingFlag, options.verbose);
  if (status < 0) {
  	fprintf(stderr, "*** ERROR: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know total # parameters
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  if (nParamsTot != (int)parameterList.size()) {
  	fprintf(stderr, "*** ERROR: number of input parameters (%d) does not equal", 
  	       (int)parameterList.size());
  	fprintf(stderr, " required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	exit(-1);
  }
  
  // Add PSF image vector, if present (needs to be added prior to image data, so that
  // ModelObject can figure out proper internal model-image size
  if (options.psfImagePresent) {
    status = theModel->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddPSFVector!\n\n");
  	  exit(-1);
    }
  }

  // Add image data and useful information about image (gain, read noise, t_exp, etc.)
  status = theModel->AddImageDataVector(allPixels, nColumns, nRows);
  if (status < 0) {
    // Possible failure if attempt to allocate memory for model image fails
    fprintf(stderr, "*** ERROR: Failure in ModelObject::AddImageDataVector!\n\n");
    exit(-1);
  }
  theModel->AddImageCharacteristics(options.gain, options.readNoise, options.expTime, options.nCombined,
  							options.originalSky);
  theModel->PrintDescription();
  if (options.printImages)
    theModel->PrintInputImage();

  // Add oversampled PSF image vector and corresponding info, if present
  // (this operates on a sub-region of the main image, so ModelObject does not need
  // to know about this prior to the image data)
  if (options.psfOversampledImagePresent) {
    status = theModel->AddOversampledPSFVector(nPixels_psf_oversampled, nColumns_psf_oversampled, 
    			nRows_psf_oversampled, psfOversampledPixels, options.psfOversamplingScale,
    			x1_oversample, x2_oversample, y1_oversample, y2_oversample);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddOversampledPSFVector!\n\n");
  	  exit(-1);
    }
  }

  // If user supplied a mask image, add it and apply it to the internal weight image
  if (maskAllocated) {
    status = theModel->AddMaskVector(nPixels_tot, nColumns, nRows, allMaskPixels,
                             options.maskFormat);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddMaskVector!\n\n");
  	  exit(-1);
    }
  }
  
  // Specify which fit statistic we'll use, and add user-supplied error image if
  // it exists and we're using chi^2; also catch special case of standard Cash
  // statistic + L-M minimizer
  if (options.useCashStatistic) {
    if ((options.solver == MPFIT_SOLVER) && (! options.printFitStatisticOnly)) {
      fprintf(stderr, "*** ERROR -- Cash statistic cannot be used with L-M solver!\n\n");
      return -1;
    }
    status = theModel->UseCashStatistic();
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::UseCashStatistic!\n\n");
      exit(-1);
    }
  } 
  else if (options.usePoissonMLR) {
    theModel->UsePoissonMLR();
  }
  else {
    // normal chi^2 statistics, so we either add error/noise image, or calculate it
    if (options.noiseImagePresent)
      theModel->AddErrorVector(nPixels_tot, nColumns, nRows, allErrorPixels,
                               options.errorType);
    else {
      if (options.useModelForErrors) {
        printf("* No noise image supplied ... will generate noise image from model image.\n");
        status = theModel->UseModelErrors();
        if (status < 0) {
          fprintf(stderr, "*** ERROR: Failure in ModelObject::UseModelErrors!\n\n");
          exit(-1);
        }
      }
      else {
        // default mode
        printf("* No noise image supplied ... will generate noise image from input data image.\n");
      }
    }
  }
  
  // Final fitting-oriented setup for ModelObject instance (generates data-based error
  // vector if needed, created final weight vector from mask and optionally from
  // error vector)
  status = theModel->FinalSetupForFitting();
  if (status < 0) {
    fprintf(stderr, "*** ERROR: Failure in ModelObject::FinalSetupForFitting!\n\n");
    exit(-1);
  }
    


  
  /* START OF MINIMIZATION-ROUTINE-RELATED CODE */
  // Parameter limits and other info:
  // First we create a C-style array of mp_par structures, containing parameter constraints
  // (if any) *and* any other useful info (like X0,Y0 offset values).  This will be used
  // by DiffEvolnFit (if called) and by PrintResults.  We also decrement nFreeParams for
  // each *fixed* parameter.
  printf("Setting up parameter information array ...\n");
  parameterInfo = (mp_par *) calloc((size_t)nParamsTot, sizeof(mp_par));
  parameterInfo_allocated = true;
  for (int i = 0; i < nParamsTot; i++) {
    parameterInfo[i].fixed = paramLimits[i].fixed;
    if (parameterInfo[i].fixed == 1) {
      nFreeParams--;
    }
    parameterInfo[i].limited[0] = paramLimits[i].limited[0];
    parameterInfo[i].limited[1] = paramLimits[i].limited[1];
    parameterInfo[i].limits[0] = paramLimits[i].limits[0];
    parameterInfo[i].limits[1] = paramLimits[i].limits[1];
    // specify different offsets if using image subsection, and apply them to
    // user-specified X0,Y0 limits
    if (theModel->GetParameterName(i) == X0_string) {
      parameterInfo[i].offset = X0_offset;
      parameterInfo[i].limits[0] -= X0_offset;
      parameterInfo[i].limits[1] -= X0_offset;
    } else if (theModel->GetParameterName(i) == Y0_string) {
      parameterInfo[i].offset = Y0_offset;
      parameterInfo[i].limits[0] -= Y0_offset;
      parameterInfo[i].limits[1] -= Y0_offset;
    }
  }
  nDegFreedom = theModel->GetNValidPixels() - nFreeParams;
  printf("%d free parameters (%ld degrees of freedom)\n", nFreeParams, nDegFreedom);


  // Now that we know all about the model (including nFreeParams), estimate the
  // memory usage and warn if it will be large
  long  estimatedMemory;
  double  nGBytes;
  bool  usingLevMar, usingCashTerms;
  if (options.solver == MPFIT_SOLVER)
    usingLevMar = true;
  else
    usingLevMar = false;
  if ((options.useCashStatistic) || (options.usePoissonMLR))
    usingCashTerms = true;
  else
    usingCashTerms = false;
  if (options.psfOversampledImagePresent) {
    int  deltaX_oversampled = x2_oversample - x1_oversample + 1;
    int  deltaY_oversampled = y2_oversample - y1_oversample + 1;
    estimatedMemory = EstimateMemoryUse(nColumns, nRows, nColumns_psf, nRows_psf, nFreeParams,
  										usingLevMar, usingCashTerms, options.saveResidualImage, 
  										options.saveModel, nColumns_psf_oversampled, nRows_psf_oversampled,
  										deltaX_oversampled, deltaY_oversampled, 
  										options.psfOversamplingScale);
  }
  else
    estimatedMemory = EstimateMemoryUse(nColumns, nRows, nColumns_psf, nRows_psf, nFreeParams,
  										usingLevMar, usingCashTerms, options.saveResidualImage, 
  										options.saveModel);
  nGBytes = (1.0*estimatedMemory) / GIGABYTE;
  if (nGBytes >= 1.0)
    printf("Estimated memory use: %ld bytes (%.1f GB)\n", estimatedMemory, nGBytes);
  else if (nGBytes >= 1.0e-3)
    printf("Estimated memory use: %ld bytes (%.1f MB)\n", estimatedMemory, nGBytes*1024.0);
  else
    printf("Estimated memory use: %ld bytes (%.1f KB)\n", estimatedMemory, nGBytes*1024.0*1024.0);
  if (estimatedMemory > MEMORY_WARNING_LIMT) {
    fprintf(stderr, "WARNING: Estimated memory needed by internal images =");
    fprintf(stderr, " %ld bytes (%g gigabytes)\n", estimatedMemory, nGBytes);
  }

  
  /* Copy initial parameter values into C array, correcting for X0,Y0 offsets */
  paramsVect = (double *) calloc(nParamsTot, sizeof(double));
  for (int i = 0; i < nParamsTot; i++) {
    if (theModel->GetParameterName(i) == X0_string) {
      paramsVect[i] = parameterList[i] - X0_offset;
    } else if (theModel->GetParameterName(i) == Y0_string) {
      paramsVect[i] = parameterList[i] - Y0_offset;
    } else
      paramsVect[i] = parameterList[i];
  }
  
  
  // OK, now we either print chi^2 value for the input parameters and quit, or
  // else call one of the solvers!
  if (options.printFitStatisticOnly) {
    printf("\n");
    PrintFitStatistic(paramsVect, theModel, nFreeParams);
    printf("\n");
    options.saveBestFitParams = false;
  }
  else {
    // DO THE FIT!
    printf("\nPerforming fit by minimizing ");
    if (options.useCashStatistic)
      printf("Cash statistic:\n");
    else if (options.usePoissonMLR)
      printf("Poisson MLR statistic:\n");
    else if (options.useModelForErrors)
      printf("chi^2 (model-based errors):\n");
    else
      printf("chi^2 (data-based errors):\n");
    gettimeofday(&timer_start_fit, NULL);
    fitStatus = DispatchToSolver(options.solver, nParamsTot, nFreeParams, nPixels_tot, 
    							paramsVect, parameterInfo, theModel, options.ftol, paramLimitsExist, 
    							options.verbose, &resultsFromSolver, options.nloptSolverName,
    							options.rngSeed);
    gettimeofday(&timer_end_fit, NULL);
    							
    PrintResults(paramsVect, theModel, nFreeParams, parameterInfo, fitStatus, resultsFromSolver);
  }


  // Optional bootstrap resampling
  if ((options.doBootstrap) && (options.bootstrapIterations > 0)) {
    if (options.saveBootstrap) {
      bootstrapSaveFile_ptr = fopen(options.outputBootstrapFileName.c_str(), "w");
      // write general info + best-fitting params as a commented-out header
      SaveParameters2(bootstrapSaveFile_ptr, paramsVect, theModel, parameterInfo, 
      				programHeader, "#");
    }
    
    printf("\nNow doing bootstrap resampling (%d iterations) to estimate errors...\n",
           options.bootstrapIterations);
    gettimeofday(&timer_start_bootstrap, NULL);
    nSucessfulIterations = BootstrapErrors(paramsVect, parameterInfo, paramLimitsExist, 
    									theModel, options.ftol, options.bootstrapIterations, 
    									nFreeParams, theModel->WhichFitStatistic(), 
    									bootstrapSaveFile_ptr, options.rngSeed);
    gettimeofday(&timer_end_bootstrap, NULL);
    if (options.saveBootstrap) {
      if (nSucessfulIterations > 0)
        printf("Bootstrap-resampling output saved to file \"%s\".\n", options.outputBootstrapFileName.c_str());
      // even if we didn't have any successful iterations, we need to close the file
      fclose(bootstrapSaveFile_ptr);
    }
    didBootstrap = true;
  }


  // Handle assorted output requests
  // Note that from this point on, we handle failures reported by SaveVectorAsImage as
  // "warnings" and don't immediately exit, since we're close to the end of the program
  // anyway, and the user might just have given us a bad path for one of the output images
  if (options.saveBestFitParams) {
    printf("Saving best-fit parameters in file \"%s\"\n", options.outputParameterFileName.c_str());
    SaveParameters(paramsVect, theModel, parameterInfo, options.outputParameterFileName,
    								programHeader, nFreeParams, options.solver, 
    								fitStatus, resultsFromSolver);
  }
  if (options.saveModel) {
    PrepareImageComments(&imageCommentsList, progNameVersion, options.outputParameterFileName,
    					options.psfImagePresent, options.psfFileName, HDR_MODELIMAGE,
    					options.imageFileName);
    printf("Saving model image in file \"%s\"\n", options.outputModelFileName.c_str());
    status = SaveVectorAsImage(theModel->GetModelImageVector(), options.outputModelFileName, 
                      nColumns, nRows, imageCommentsList);
    if (status != 0) {
      fprintf(stderr, "\n*** WARNING: Failure saving model-image file \"%s\"!\n\n",
      				options.outputModelFileName.c_str());
    }
  }
  if (options.saveResidualImage) {
    imageCommentsList.clear();
    PrepareImageComments(&imageCommentsList, progNameVersion, options.outputParameterFileName,
    					options.psfImagePresent, options.psfFileName, HDR_RESIDUALIMAGE,
    					options.imageFileName);
    printf("Saving residual (data - best-fit model) image in file \"%s\"\n", options.outputResidualFileName.c_str());
    status = SaveVectorAsImage(theModel->GetResidualImageVector(), options.outputResidualFileName, 
                      nColumns, nRows, imageCommentsList);
    if (status != 0) {
      fprintf(stderr, "\n*** WARNING: Failure saving residual-image file \"%s\"!\n\n",
      				options.outputResidualFileName.c_str());
    }
  }
  if (options.saveWeightImage) {
    imageCommentsList.clear();
    PrepareImageComments(&imageCommentsList, progNameVersion, options.outputParameterFileName,
    					options.psfImagePresent, options.psfFileName, HDR_WEIGHTIMAGE,
    					options.imageFileName);
    printf("Saving weight image in file \"%s\"\n", options.outputWeightFileName.c_str());
    status = SaveVectorAsImage(theModel->GetWeightImageVector(), options.outputWeightFileName, 
                      nColumns, nRows, imageCommentsList);
    if (status != 0) {
      fprintf(stderr, "\n*** WARNING: Failure saving weight-image file \"%s\"!\n\n",
      				options.outputWeightFileName.c_str());
    }
  }


  // Free up memory
  fftw_free(allPixels);                 // allocated externally, in ReadImageAsVector()
  if (errorPixels_allocated)
    fftw_free(allErrorPixels);          // allocated externally, in ReadImageAsVector()
  if (options.psfImagePresent)
    fftw_free(psfPixels);               // allocated externally, in ReadImageAsVector()
  if (options.psfOversampledImagePresent)
    fftw_free(psfOversampledPixels);    // allocated externally, in ReadImageAsVector()
  if (maskAllocated)
    fftw_free(allMaskPixels);           // allocated externally, in ReadImageAsVector()
  free(paramsVect);
  if (parameterInfo_allocated)
    free(parameterInfo);
  delete theModel;
  
  // Elapsed time reports
  if (options.verbose >= 0) {
    gettimeofday(&timer_end_all, NULL);
    double  microsecs, time_elapsed_all, time_elapsed_fit, time_elapsed_bootstrap;
    microsecs = timer_end_all.tv_usec - timer_start_all.tv_usec;
    time_elapsed_all = timer_end_all.tv_sec - timer_start_all.tv_sec + microsecs/1e6;
    if (options.printFitStatisticOnly)
      printf("\n(Elapsed time: %.6f sec)\n", time_elapsed_all);
    else {
      microsecs = timer_end_fit.tv_usec - timer_start_fit.tv_usec;
      time_elapsed_fit = timer_end_fit.tv_sec - timer_start_fit.tv_sec + microsecs/1e6;
      if (didBootstrap) {
        microsecs = timer_end_bootstrap.tv_usec - timer_start_bootstrap.tv_usec;
        time_elapsed_bootstrap = timer_end_bootstrap.tv_sec - timer_start_bootstrap.tv_sec + microsecs/1e6;
        printf("\n(Elapsed time: %.6f sec for fit, %.6f for bootstrap, %.6f sec total)\n", 
        		time_elapsed_fit, time_elapsed_bootstrap, time_elapsed_all);
      }
      else
        printf("\n(Elapsed time: %.6f sec for fit, %.6f sec total)\n", time_elapsed_fit, 
        		time_elapsed_all);
    }
  }
  
  printf("Done!\n\n");
  
  return 0;
}



void ProcessInput( int argc, char *argv[], imfitCommandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();
  string  tempString = "";

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   imfit [options] <imagefile.fits>");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine(" -v  --version                Prints version number");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --list-parameters        Prints list of parameter names for each available function");
  tempString = PrintToString("     --sample-config          Generates an example configuration file (%s)", configImfitFile.c_str());
  optParser->AddUsageLine(tempString);
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" -c  --config <config-file>   configuration file [REQUIRED!]");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --noise <noisemap.fits>  Noise/error/weight image to use");
  optParser->AddUsageLine("     --mask <mask.fits>       Mask image to use");
  optParser->AddUsageLine("     --psf <psf.fits>         PSF image to use");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --overpsf <psf.fits>      Oversampled PSF image to use");
  optParser->AddUsageLine("     --overpsf_scale <n>       Oversampling scale (integer)");
  optParser->AddUsageLine("     --overpsf_region <x1:x2,y1:y2>       Section of image to convolve with oversampled PSF");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --save-params <output-file>          Specify filename for best-fit parameters output [default = bestfit_parameters_imfit.dat]");
  optParser->AddUsageLine("     --save-model <outputname.fits>       Save best-fit model image");
  optParser->AddUsageLine("     --save-residual <outputname.fits>    Save residual (data - best-fit model) image");
  optParser->AddUsageLine("     --save-weights <outputname.fits>     Save weight image");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --sky <sky-level>        Original sky background (ADUs) which was subtracted from image");
  optParser->AddUsageLine("     --gain <value>           Image A/D gain (e-/ADU)");
  optParser->AddUsageLine("     --readnoise <value>      Image read noise (e-)");
  optParser->AddUsageLine("     --exptime <value>        Exposure time in sec (only if image counts are ADU/sec)");
  optParser->AddUsageLine("     --ncombined <value>      Number of images averaged to make final image (if counts are average or median)");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --errors-are-variances   Indicates that values in noise image = variances (instead of sigmas)");
  optParser->AddUsageLine("     --errors-are-weights     Indicates that values in noise image = weights (instead of sigmas)");
  optParser->AddUsageLine("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --model-errors           Use model values (instead of data) to estimate errors for chi^2 computation");
  optParser->AddUsageLine("     --cashstat               Use Cash statistic instead of chi^2");
  optParser->AddUsageLine("     --poisson-mlr            Use Poisson maximum-likelihood-ratio statistic instead of chi^2");
  optParser->AddUsageLine("     --mlr                    Same as --poisson-mlr");
  optParser->AddUsageLine("     --ftol                   Fractional tolerance in fit statistic for convergence [default = 1.0e-8]");
  optParser->AddUsageLine("");
#ifndef NO_NLOPT
  optParser->AddUsageLine("     --nm                     Use Nelder-Mead simplex solver (instead of Levenberg-Marquardt)");
  optParser->AddUsageLine("     --nlopt <name>           Select miscellaneous NLopt solver");
#endif
  optParser->AddUsageLine("     --de                     Use differential evolution solver");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --bootstrap <int>        Do this many iterations of bootstrap resampling to estimate errors");
  optParser->AddUsageLine("     --save-bootstrap <filename>        Save all bootstrap best-fit parameters to specified file");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --chisquare-only         Print fit statistic (e.g., chi^2) of input model and quit (no fitting done)");
  optParser->AddUsageLine("     --fitstat-only           Same as --chisquare-only");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --quiet                  Turn off printing of updates during the fit");
  optParser->AddUsageLine("     --silent                 Turn off ALL printouts (except fatal errors)");
  optParser->AddUsageLine("     --loud                   Print extra info during the fit");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --max-threads <int>      Maximum number of threads to use");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --seed <int>             RNG seed (for testing purposes)");
  optParser->AddUsageLine("     --nosubsampling          Turn off pixel subsampling near centers of functions");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("EXAMPLES:");
  optParser->AddUsageLine("   imfit -c model_config_n100a.dat ngc100.fits");
  optParser->AddUsageLine("   imfit -c model_config_n100b.dat ngc100.fits[405:700,844:1060] --mask ngc100_mask.fits[405:700,844:1060] --gain 4.5 --readnoise 0.7");
  optParser->AddUsageLine("");


  /* by default all options are checked on the command line and from option/resource file */
  optParser->AddFlag("help", "h");
  optParser->AddFlag("version", "v");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("list-parameters");
  optParser->AddFlag("sample-config");
  optParser->AddFlag("printimage");
  optParser->AddFlag("chisquare-only");
  optParser->AddFlag("fitstat-only");
  optParser->AddFlag("errors-are-variances");
  optParser->AddFlag("errors-are-weights");
  optParser->AddFlag("mask-zero-is-bad");
  optParser->AddFlag("nosubsampling");
  optParser->AddFlag("model-errors");
  optParser->AddFlag("cashstat");
  optParser->AddFlag("poisson-mlr");
  optParser->AddFlag("mlr");
#ifndef NO_NLOPT
  optParser->AddFlag("nm");
  optParser->AddOption("nlopt");
#endif
  optParser->AddFlag("de");
  optParser->AddFlag("quiet");
  optParser->AddFlag("silent");
  optParser->AddFlag("loud");
  optParser->AddOption("noise");
  optParser->AddOption("mask");
  optParser->AddOption("psf");
  optParser->AddOption("overpsf");
  optParser->AddOption("overpsf_scale");
  optParser->AddOption("overpsf_region");
  optParser->AddOption("save-params");
  optParser->AddOption("save-model");
  optParser->AddOption("save-residual");
  optParser->AddOption("save-weights");
  optParser->AddOption("sky");
  optParser->AddOption("gain");
  optParser->AddOption("readnoise");
  optParser->AddOption("exptime");
  optParser->AddOption("ncombined");
  optParser->AddOption("ftol");
  optParser->AddOption("bootstrap");
  optParser->AddOption("save-bootstrap");
  optParser->AddOption("config", "c");
  optParser->AddOption("max-threads");
  optParser->AddOption("seed");

  // Comment this out if you want unrecognized (e.g., mis-spelled) flags and options
  // to be ignored only, rather than causing program to exit
  optParser->UnrecognizedAreErrors();
  
  /* parse the command line:  */
  int status = optParser->ParseCommandLine( argc, argv );
  if (status < 0) {
    printf("\nError on command line... quitting...\n\n");
    delete optParser;
    exit(1);
  }


  /* Process the results: actual arguments, if any: */
  if (optParser->nArguments() > 0) {
    theOptions->imageFileName = optParser->GetArgument(0);
    theOptions->noImage = false;
    printf("\tImage file = %s\n", theOptions->imageFileName.c_str());
  }

  /* Process the results: options */
  // First four are options which print useful info and then exit the program
  if ( optParser->FlagSet("help") || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("version") ) {
    printf("imfit version %s\n\n", VERSION_STRING);
    delete optParser;
    exit(1);
  }
  if (optParser->FlagSet("list-functions")) {
    PrintAvailableFunctions();
    delete optParser;
    exit(1);
  }
  if (optParser->FlagSet("list-parameters")) {
    ListFunctionParameters();
    delete optParser;
    exit(1);
  }
  if (optParser->FlagSet("sample-config")) {
    int saveStatus = SaveExampleImfitConfig();
    if (saveStatus == 0)
      printf("Sample configuration file \"%s\" saved.\n", configImfitFile.c_str());
    delete optParser;
    exit(1);
  }

  if (optParser->FlagSet("printimage")) {
    theOptions->printImages = true;
  }
  if (optParser->FlagSet("chisquare-only")) {
    printf("\t* No fitting will be done!\n");
    theOptions->printFitStatisticOnly = true;
  }
  if (optParser->FlagSet("fitstat-only")) {
    printf("\t* No fitting will be done!\n");
    theOptions->printFitStatisticOnly = true;
  }
  if (optParser->FlagSet("model-errors")) {
  	printf("\t* Using model counts instead of data to compute errors for chi^2\n");
  	theOptions->useModelForErrors = true;
  }
  if (optParser->FlagSet("cashstat")) {
  	printf("\t* Using standard Cash statistic instead of chi^2 for minimization!\n");
  	theOptions->useCashStatistic = true;
  }
  if ( (optParser->FlagSet("poisson-mlr")) || (optParser->FlagSet("mlr")) ) {
  	printf("\t* Using Poisson maximum-likelihood-ratio statistic instead of chi^2 for minimization!\n");
  	theOptions->usePoissonMLR = true;
  }
#ifndef NO_NLOPT
  if (optParser->FlagSet("nm")) {
  	printf("\t* Nelder-Mead simplex solver selected!\n");
  	theOptions->solver = NMSIMPLEX_SOLVER;
  }
  if (optParser->OptionSet("nlopt")) {
    theOptions->solver = GENERIC_NLOPT_SOLVER;
    theOptions->nloptSolverName = optParser->GetTargetString("nlopt");
    if (! ValidNLOptSolverName(theOptions->nloptSolverName)) {
      fprintf(stderr, "*** ERROR: \"%s\" is not a valid NLOpt solver name!\n", 
      			theOptions->nloptSolverName.c_str());
      fprintf(stderr, "    (valid names for --nlopt: COBYLA, BOBYQA, NEWUOA, PRAXIS, NM, SBPLX)\n");
      delete optParser;
      exit(1);
    }
    printf("\tNLopt solver = %s\n", theOptions->nloptSolverName.c_str());
  }
#endif
  if (optParser->FlagSet("de")) {
  	printf("\t* Differential Evolution selected!\n");
  	theOptions->solver = DIFF_EVOLN_SOLVER;
  }
  if (optParser->FlagSet("nosubsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (optParser->FlagSet("silent")) {
    theOptions->verbose = -1;
  }
  if (optParser->FlagSet("quiet")) {
    theOptions->verbose = 0;
  }
  if (optParser->FlagSet("loud")) {
    theOptions->verbose = 2;
  }
  if (optParser->FlagSet("errors-are-variances")) {
    theOptions->errorType = WEIGHTS_ARE_VARIANCES;
  }
  if (optParser->FlagSet("errors-are-weights")) {
    theOptions->errorType = WEIGHTS_ARE_WEIGHTS;
  }
  if (optParser->FlagSet("mask-zero-is-bad")) {
    theOptions->maskFormat = MASK_ZERO_IS_BAD;
  }
  if (optParser->OptionSet("config")) {
    theOptions->configFileName = optParser->GetTargetString("config");
    printf("\tconfiguration file = %s\n", theOptions->configFileName.c_str());
  }
  if (optParser->OptionSet("noise")) {
    theOptions->noiseFileName = optParser->GetTargetString("noise");
    theOptions->noiseImagePresent = true;
    printf("\tnoise image = %s\n", theOptions->noiseFileName.c_str());
  }
  if (optParser->OptionSet("psf")) {
    theOptions->psfFileName = optParser->GetTargetString("psf");
    theOptions->psfImagePresent = true;
    printf("\tPSF image = %s\n", theOptions->psfFileName.c_str());
  }
  if (optParser->OptionSet("overpsf")) {
    theOptions->psfOversampledFileName = optParser->GetTargetString("overpsf");
    theOptions->psfOversampledImagePresent = true;
    printf("\tOversampled PSF image = %s\n", theOptions->psfOversampledFileName.c_str());
  }
  if (optParser->OptionSet("overpsf_scale")) {
    if (NotANumber(optParser->GetTargetString("overpsf_scale").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: overpsf_scale should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->psfOversamplingScale = atoi(optParser->GetTargetString("overpsf_scale").c_str());
    printf("\tPSF oversampling scale = %d\n", theOptions->psfOversamplingScale);
  }
  if (optParser->OptionSet("overpsf_region")) {
    theOptions->psfOversampleRegion = optParser->GetTargetString("overpsf_region");
    theOptions->oversampleRegionSet = true;
    theOptions->nOversampleRegions += 1;
    printf("\tPSF oversampling region = %s\n", theOptions->psfOversampleRegion.c_str());
  }
  if (optParser->OptionSet("mask")) {
    theOptions->maskFileName = optParser->GetTargetString("mask");
    theOptions->maskImagePresent = true;
    printf("\tmask image = %s\n", theOptions->maskFileName.c_str());
  }
  if (optParser->OptionSet("save-model")) {
    theOptions->outputModelFileName = optParser->GetTargetString("save-model");
    theOptions->saveModel = true;
    printf("\toutput best-fit model image = %s\n", theOptions->outputModelFileName.c_str());
  }
  if (optParser->OptionSet("save-residual")) {
    theOptions->outputResidualFileName = optParser->GetTargetString("save-residual");
    theOptions->saveResidualImage = true;
    printf("\toutput residual (data - best-fit model) image = %s\n", theOptions->outputResidualFileName.c_str());
  }
  if (optParser->OptionSet("save-weights")) {
    theOptions->outputWeightFileName = optParser->GetTargetString("save-weights");
    theOptions->saveWeightImage = true;
    printf("\toutput weight image = %s\n", theOptions->outputWeightFileName.c_str());
  }
  if (optParser->OptionSet("save-params")) {
    theOptions->outputParameterFileName = optParser->GetTargetString("save-params");
    theOptions->saveBestFitParams = true;
    printf("\toutput best-fit parameter file = %s\n", theOptions->outputParameterFileName.c_str());
  }
  if (optParser->OptionSet("sky")) {
    if (NotANumber(optParser->GetTargetString("sky").c_str(), 0, kAnyReal)) {
      fprintf(stderr, "*** ERROR: sky should be a real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->originalSky = atof(optParser->GetTargetString("sky").c_str());
    theOptions->originalSkySet = true;
    printf("\toriginal sky level = %g ADU\n", theOptions->originalSky);
  }
  if (optParser->OptionSet("gain")) {
    if (NotANumber(optParser->GetTargetString("gain").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** ERROR: gain should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->gain = atof(optParser->GetTargetString("gain").c_str());
    theOptions->gainSet = true;
    printf("\tgain = %g e-/ADU\n", theOptions->gain);
  }
  if (optParser->OptionSet("readnoise")) {
    if (NotANumber(optParser->GetTargetString("readnoise").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** ERROR: read noise should be a non-negative real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->readNoise = atof(optParser->GetTargetString("readnoise").c_str());
    theOptions->readNoiseSet = true;
    printf("\tread noise = %g e-\n", theOptions->readNoise);
  }
  if (optParser->OptionSet("exptime")) {
    if (NotANumber(optParser->GetTargetString("exptime").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** ERROR: exptime should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->expTime = atof(optParser->GetTargetString("exptime").c_str());
    theOptions->expTimeSet = true;
    printf("\texposure time = %g sec\n", theOptions->expTime);
  }
  if (optParser->OptionSet("ncombined")) {
    if (NotANumber(optParser->GetTargetString("ncombined").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: ncombined should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nCombined = atoi(optParser->GetTargetString("ncombined").c_str());
    theOptions->nCombinedSet = true;
    printf("\tn_combined = %d\n", theOptions->nCombined);
  }
  if (optParser->OptionSet("ftol")) {
    if (NotANumber(optParser->GetTargetString("ftol").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** ERROR: ftol should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->ftol = atof(optParser->GetTargetString("ftol").c_str());
    theOptions->ftolSet = true;
    printf("\tfractional tolerance ftol for fit-statistic convergence = %g\n", theOptions->ftol);
  }
  if (optParser->OptionSet("bootstrap")) {
    if (NotANumber(optParser->GetTargetString("bootstrap").c_str(), 0, kPosInt)) {
      printf("*** ERROR: number of bootstrap iterations should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->doBootstrap = true;
    theOptions->bootstrapIterations = atol(optParser->GetTargetString("bootstrap").c_str());
    printf("\tnumber of bootstrap iterations = %d\n", theOptions->bootstrapIterations);
  }
  if (optParser->OptionSet("save-bootstrap")) {
    theOptions->outputBootstrapFileName = optParser->GetTargetString("save-bootstrap");
    theOptions->saveBootstrap = true;
    printf("\tbootstrap best-fit parameters to be saved in %s\n", theOptions->outputBootstrapFileName.c_str());
  }
  if (optParser->OptionSet("max-threads")) {
    if (NotANumber(optParser->GetTargetString("max-threads").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: max-threads should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->maxThreads = atol(optParser->GetTargetString("max-threads").c_str());
    theOptions->maxThreadsSet = true;
  }
  if (optParser->OptionSet("seed")) {
    if (NotANumber(optParser->GetTargetString("seed").c_str(), 0, kPosInt)) {
      printf("*** WARNING: RNG seed should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->rngSeed = atol(optParser->GetTargetString("seed").c_str());
    printf("\tRNG seed = %ld\n", theOptions->rngSeed);
  }

  delete optParser;

}



/// Checks to see that all user-requested files are present; returns false if
/// any are missing, and prints appropriate error messages.
/// Files we check:
///    config file (options.configFileName)
///    data image (options.imageFileName)
/// and the following, if the user supplied names for them:
///    mask image (options.maskFileName)
///    noise image (options.noiseFileName)
///    PSF image (options.psfFileName)
bool RequestedFilesPresent( imfitCommandOptions &theOptions )
{
  bool  allFilesPresent = true;
  
  if (! FileExists(theOptions.configFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find configuration file \"%s\"!\n", 
           theOptions.configFileName.c_str());
    allFilesPresent = false;
  }
  if (! ImageFileExists(theOptions.imageFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find image file \"%s\"!\n", 
           theOptions.imageFileName.c_str());
    allFilesPresent = false;
  }
  if ( (theOptions.maskImagePresent) && (! ImageFileExists(theOptions.maskFileName.c_str())) ) {
    fprintf(stderr, "\n*** ERROR: Unable to find mask file \"%s\"!\n", 
           theOptions.maskFileName.c_str());
    allFilesPresent = false;
  }
  if ( (theOptions.noiseImagePresent) && (! ImageFileExists(theOptions.noiseFileName.c_str())) ) {
    fprintf(stderr, "\n*** ERROR: Unable to find noise-image file \"%s\"!\n", 
           theOptions.noiseFileName.c_str());
    allFilesPresent = false;
  }
  if ( (theOptions.psfImagePresent) && (! ImageFileExists(theOptions.psfFileName.c_str())) ) {
    fprintf(stderr, "\n*** ERROR: Unable to find PSF image file \"%s\"!\n", 
           theOptions.psfFileName.c_str());
    allFilesPresent = false;
  }

  return allFilesPresent;
}



// Note that we only use options from the config file if they have *not*
// already been set by the command line (i.e., command-line options override
// config-file values).
void HandleConfigFileOptions( configOptions *configFileOptions, imfitCommandOptions *mainOptions )
{
	double  newDblVal;
	int  newIntVal;
	
  if (configFileOptions->nOptions == 0)
    return;

  for (int i = 0; i < configFileOptions->nOptions; i++) {
    
    if (configFileOptions->optionNames[i] == kGainString) {
      if (mainOptions->gainSet) {
        printf("Gain value in config file ignored (using command-line value)\n");
      } else {
        newDblVal = strtod(configFileOptions->optionValues[i].c_str(), NULL);
        printf("Value from config file: gain = %f e-/ADU\n", newDblVal);
        mainOptions->gain = newDblVal;
      }
      continue;
    }
    if (configFileOptions->optionNames[i] == kReadNoiseString) {
      if (mainOptions->readNoiseSet) {
        printf("Read-noise value in config file ignored (using command-line value)\n");
      } else {
        newDblVal = strtod(configFileOptions->optionValues[i].c_str(), NULL);
        printf("Value from config file: read noise = %f e-\n", newDblVal);
        mainOptions->readNoise = newDblVal;
      }
      continue;
    }
    if (configFileOptions->optionNames[i] == kExpTimeString) {
      if (mainOptions->expTimeSet) {
        printf("Read-noise value in config file ignored (using command-line value)\n");
      } else {
        newDblVal = strtod(configFileOptions->optionValues[i].c_str(), NULL);
        printf("Value from config file: exposure time = %f sec\n", newDblVal);
        mainOptions->expTime = newDblVal;
      }
      continue;
    }
    if (configFileOptions->optionNames[i] == kOriginalSkyString) {
      if (mainOptions->originalSkySet) {
        printf("Original-sky value in config file ignored (using command-line value)\n");
      } else {
        newDblVal = strtod(configFileOptions->optionValues[i].c_str(), NULL);
        printf("Value from config file: original sky = %f\n", newDblVal);
        mainOptions->originalSky = newDblVal;
      }
      continue;
    }
    if (configFileOptions->optionNames[i] == kNCombinedString) {
      if (mainOptions->nCombinedSet) {
        printf("nCombined value in config file ignored (using command-line value)\n");
      } else {
        newIntVal = atoi(configFileOptions->optionValues[i].c_str());
        printf("Value from config file: nCombined = %d\n", newIntVal);
        mainOptions->nCombined = newIntVal;
      }
      continue;
    }
    // we only get here if we encounter an unknown option
    printf("Unknown keyword (\"%s\") in config file ignored\n", 
    				configFileOptions->optionNames[i].c_str());
    
  }
}




/// Takes user-supplied image filename and determines what, if any, x0 and y0
/// pixel offsets are implied by any section specification in the filename
/// (e.g. "somefile.fits[100:250,400:400]").
/// Note that offsets are always >= 0.
void DetermineImageOffset( const std::string &fullImageName, double *x_offset,
					double *y_offset)
{
  int  xStart, yStart;

  GetPixelStartCoords(fullImageName, &xStart, &yStart);
  *x_offset = xStart - 1;
  *y_offset = yStart - 1;
}



/* END OF FILE: imfit_main.cpp ------------------------------------------- */
