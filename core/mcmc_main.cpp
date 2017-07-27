/* FILE: mcmc_main.cpp --------------------------------------------------- */
/*
 * This is the main program file for imfit-mcmc.
 *
 * Useful reminder about FITS image sizes -- the proper translations are:
 * NAXIS1 = naxes[0] = nColumns = sizeX;
 * NAXIS2 = naxes[1] = nRows = sizeY.
 *
 *
 * HISTORY
 *    24 October 2016: Created as modification of imfit_main.cpp.
*/

// Copyright 2009--2017 by Peter Erwin.
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
#include "getimages.h"
#include "model_object.h"
#include "add_functions.h"
#include "param_struct.h"   // for mp_par structure
#include "options_base.h"
#include "options_mcmc.h"
#include "psf_oversampling_info.h"
#include "setup_model_object.h"

#include "commandline_parser.h"
#include "config_file_parser.h"
#include "estimate_memory.h"
#include "sample_configs.h"

#include "dream_params.h"
#include "dream.h"
#include "rng/GSLStream.h"

#ifdef USE_PLOG
#include "plog/Log.h"
#endif



/* ---------------- Definitions & Constants ----------------------------- */

// Option names for use in config files
static string  kGainString = "GAIN";
static string  kReadNoiseString = "READNOISE";
static string  kExpTimeString = "EXPTIME";
static string  kNCombinedString = "NCOMBINED";
static string  kOriginalSkyString = "ORIGINAL_SKY";

const char *  LOG_FILENAME = "log_imfit-mcmc.txt";


#ifdef USE_OPENMP
#define VERSION_STRING      "1.5.0a (OpenMP-enabled)"
#else
#define VERSION_STRING      "1.5.0a"
#endif



/* ------------------- Function Prototypes ----------------------------- */
void ProcessInput( int argc, char *argv[], MCMCOptions *theOptions );
bool RequestedFilesPresent( MCMCOptions *theOptions );
void HandleConfigFileOptions( configOptions *configFileOptions, 
								MCMCOptions *mainOptions );

double LikelihoodFuncForDREAM( int chain, int gen, const double* state, 
								const void* extraData, bool recalc );
void MakeMCMCOutputHeader( vector<string> *headerLines, const string& programName, 
						const int argc, char *argv[] );


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
  PsfOversamplingInfo  *psfOversampleInfo = NULL;
  vector<PsfOversamplingInfo *>  psfOversamplingInfoVect;
  double  *paramsVect;
  int  X0_offset = 0;
  int  Y0_offset = 0;
  int  x1_oversample = 0;
  int  x2_oversample = 0;
  int  y1_oversample = 0;
  int  y2_oversample = 0;
  vector<int> xyOsamplePos;
  std::string  noiseImage;
  ModelObject  *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  FunctionBlockIndices;
  vector< map<string, string> > optionalParamsMap;
  bool  paramLimitsExist = false;
  bool  parameterInfo_allocated = false;
  mp_par  *parameterInfo;
  int  status, fitStatus, nSucessfulIterations;
  vector<string>  imageCommentsList;
  OptionsBase *commandOpts;
  MCMCOptions *options;
  configOptions  userConfigOptions;
  const std::string  X0_string("X0");
  const std::string  Y0_string("Y0");
  string  progNameVersion = "imfit-mcmc ";
  vector<string> programHeader;
  // timing-related
  struct timeval  timer_start_all, timer_end_all;


  gettimeofday(&timer_start_all, NULL);

#ifdef USE_PLOG
  plog::init(plog::debug, LOG_FILENAME, 100000, 3); // Initialize the logger.
  LOG_INFO << "\n*** Starting up...";
#endif

  progNameVersion += VERSION_STRING;
  MakeMCMCOutputHeader(&programHeader, progNameVersion, argc, argv);

 
  // Define default options, then process the command line
  // Use a pointer to OptionsBase so we can use it in calls to SetupModelImage
  commandOpts = new MCMCOptions();
  // Use a pointer to MakeimageOptions so we can access all the extra, makeimage-specific
  // data members
  options = (MCMCOptions *)commandOpts;

  ProcessInput(argc, argv, options);

  // Check for presence of user-requested files; if any are missing, quit.
  // (Appropriate error messages regarding which files are missing will be printed
  // to stderr by RequestedFilesPresent)
  if (! RequestedFilesPresent(options)) {
    fprintf(stderr, "\n");
    exit(-1);
  }


  // Read configuration file, parse & process user-supplied (non-function-related) values
  status = ReadConfigFile(options->configFileName, true, functionList, parameterList, 
  							paramLimits, FunctionBlockIndices, paramLimitsExist, userConfigOptions);
  if (status != 0) {
    fprintf(stderr, "\n*** ERROR: Failure reading configuration file!\n\n");
    return -1;
  }
  HandleConfigFileOptions(&userConfigOptions, options);

  
  if (options->noImage) {
    fprintf(stderr, "*** ERROR: No image to fit!\n\n");
    return -1;
  }

  // Get image data and sizes
  printf("Reading data image (\"%s\") ...\n", options->imageFileName.c_str());
  allPixels = ReadImageAsVector(options->imageFileName, &nColumns, &nRows);
  if (allPixels == NULL) {
    fprintf(stderr,  "\n*** ERROR: Unable to read image file \"%s\"!\n\n", 
    			options->imageFileName.c_str());
    exit(-1);
  }
  // Reminder: nColumns = n_pixels_per_row = x-size; nRows = n_pixels_per_column = y-size
  nPixels_tot = (long)nColumns * (long)nRows;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns, nRows, nPixels_tot);
  // // Determine X0,Y0 pixel offset values if user specified an image section
  DetermineImageOffset(options->imageFileName, &X0_offset, &Y0_offset);


  // Get (and check) mask and/or error images
  std::tie(allMaskPixels, allErrorPixels, status) = GetMaskAndErrorImages(nColumns, nRows, 
  													options, maskAllocated, errorPixels_allocated);
  if (status < 0)
    exit(-1);

  // Read in PSF image, if supplied
  if (options->psfImagePresent) {
    std::tie(psfPixels, nColumns_psf, nRows_psf, status) = GetPsfImage(options);
    if (status < 0)
      exit(-1);
  }
  else
    printf("* No PSF image supplied -- no image convolution will be done!\n");

  // Read in oversampled PSF image(s), if supplied
  if ((options->psfOversampling) && (options->psfOversampledImagePresent)) {
    status = GetOversampledPsfInfo(options, X0_offset, Y0_offset, psfOversamplingInfoVect);
	if (status < 0)
	  exit(-1);
  }

  // Read in oversampled PSF image, if supplied
//   if (options->psfOversampledImagePresent) {
//     if (options->psfOversamplingScale < 1) {
//       fprintf(stderr, "\n*** ERROR: the oversampling scale for the oversampled PSF was not supplied!\n");
//       fprintf(stderr, "           (use --overpsf_scale to specify the scale)\n\n");
//       exit(-1);
//     }
//     if (! options->oversampleRegionSet) {
//       fprintf(stderr, "\n*** ERROR: the oversampling region was not defined!\n");
//       fprintf(stderr, "           (use --overpsf_region to specify the region)\n\n");
//       exit(-1);
//     }
//     printf("Reading oversampled PSF image (\"%s\") ...\n", options->psfOversampledFileName.c_str());
//     psfOversampledPixels = ReadImageAsVector(options->psfOversampledFileName, 
//     							&nColumns_psf_oversampled, &nRows_psf_oversampled);
//     if (psfOversampledPixels == NULL) {
//       fprintf(stderr, "\n*** ERROR: Unable to read oversampled PSF image file \"%s\"!\n\n", 
//     			options->psfOversampledFileName.c_str());
//       exit(-1);
//     }
//     nPixels_psf_oversampled = (long)nColumns_psf_oversampled * (long)nRows_psf_oversampled;
//     printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
//            nColumns_psf_oversampled, nRows_psf_oversampled, nPixels_psf_oversampled);
//     // Determine oversampling region; correct for X0,Y0 pixel offset values if user 
//     // specified an image section (if not, X0_offset and Y0_offset will be = 0)
//     GetAllCoordsFromBracket(options->psfOversampleRegion, &x1_oversample, &x2_oversample, 
//     						&y1_oversample, &y2_oversample);
//     x1_oversample -= (int)X0_offset;
//     x2_oversample -= (int)X0_offset;
//     y1_oversample -= (int)Y0_offset;
//     y2_oversample -= (int)Y0_offset;
//     string  msg, msg0;
//     if ((X0_offset != 0.0) || (Y0_offset != 0.0)) {
//       msg0 = "\tRegion to be oversampled within full image:";
//       msg = "\tRegion to be oversampled within fitted subset image:";
//       printf("%s [%d:%d,%d:%d]\n", msg0.c_str(),
//       		x1_oversample + (int)X0_offset,x2_oversample + (int)X0_offset,
//       		y1_oversample + (int)Y0_offset,y2_oversample + (int)Y0_offset);
//     }
//     else
//       msg = "\tRegion to be oversampled within image:";
//     printf("%s [%d:%d,%d:%d]\n", msg.c_str(), x1_oversample,x2_oversample,
//     		y1_oversample,y2_oversample);
//     xyOsamplePos.push_back(x1_oversample);
//     xyOsamplePos.push_back(x2_oversample);
//     xyOsamplePos.push_back(y1_oversample);
//     xyOsamplePos.push_back(y2_oversample);
//   }

  if (! options->subsamplingFlag)
    printf("* Pixel subsampling has been turned OFF.\n");


  // Set up the model object
  // Populate the column-and-row-numbers vector
  vector<int> nColumnsRowsVect;
  nColumnsRowsVect.push_back(nColumns);
  nColumnsRowsVect.push_back(nRows);
  nColumnsRowsVect.push_back(nColumns_psf);
  nColumnsRowsVect.push_back(nRows_psf);

  theModel = SetupModelObject(options, nColumnsRowsVect, allPixels, psfPixels, allMaskPixels,
  								allErrorPixels, psfOversamplingInfoVect);



  // Add functions to the model object
  status = AddFunctions(theModel, functionList, FunctionBlockIndices, 
  						options->subsamplingFlag, options->verbose, optionalParamsMap);
  if (status < 0) {
  	fprintf(stderr, "*** ERROR: Failure in AddFunctions!\n\n");
  	exit(-1);
  }
  
  
  // Determine nParamsTot, nFreeParams and nDegFreedom
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  if (nParamsTot != (int)parameterList.size()) {
  	fprintf(stderr, "*** ERROR: number of input parameters (%d) does not equal", 
  	       (int)parameterList.size());
  	fprintf(stderr, " required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	exit(-1);
  }
  for (int i = 0; i < nParamsTot; i++) {
    if (paramLimits[i].fixed == 1)
      nFreeParams--;
  }
  nDegFreedom = (long)(theModel->GetNValidPixels() - nFreeParams);
  printf("%d free parameters (%ld degrees of freedom)\n", nFreeParams, nDegFreedom);

  theModel->PrintDescription();

  // Final fitting-oriented setup for ModelObject instance (generates data-based error
  // vector if needed, created final weight vector from mask and optionally from
  // error vector)
  status = theModel->FinalSetupForFitting();
  if (status < 0) {
    fprintf(stderr, "*** ERROR: Failure in ModelObject::FinalSetupForFitting!\n\n");
    exit(-1);
  }


  // Set up parameter-info/limits array (mostly useful for storing
  // X0_offset, Y0_offset, which will be needed for printing output
  // parameter values)
  if (nParamsTot <= 0) {
    fprintf(stderr, "*** ERROR: nParamsTot was not set correctly!\n\n");
    exit(-1);
  }
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

  // tell ModelObject about parameterInfo (mainly useful for printing-related methods)
  theModel->AddParameterInfo(parameterInfo);


 
  // Get estimate of memory use (do this *after* we know number of free parameters); 
  // warn if it will be large
  long  estimatedMemory;
  double  nGBytes;
  bool  usingLevMar, usingCashTerms;
  usingLevMar = false;
  if ((options->useCashStatistic) || (options->usePoissonMLR))
    usingCashTerms = true;
  else
    usingCashTerms = false;

  estimatedMemory = EstimateMemoryUse(nColumns, nRows, nColumns_psf, nRows_psf, nFreeParams,
										usingLevMar, usingCashTerms, options->saveResidualImage, 
  										options->saveModel);
  if (options->psfOversampledImagePresent)
    estimatedMemory += EstimatePsfOversamplingMemoryUse(psfOversamplingInfoVect);

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

  
  // Copy initial parameter values into C array, correcting for X0,Y0 offsets
  paramsVect = (double *) calloc(nParamsTot, sizeof(double));
  for (int i = 0; i < nParamsTot; i++) {
    if (theModel->GetParameterName(i) == X0_string) {
      paramsVect[i] = parameterList[i] - X0_offset;
    } else if (theModel->GetParameterName(i) == Y0_string) {
      paramsVect[i] = parameterList[i] - Y0_offset;
    } else
      paramsVect[i] = parameterList[i];
  }
  
  
  // Stuff for dream_pars
  string *paramNames = new string[nParamsTot];
  
  printf("Setting up MCMC-related arrays ...\n");
  int lockFlags[nParamsTot];
  double lowVals[nParamsTot];
  double highVals[nParamsTot];
  for (int i = 0; i < nParamsTot; i++) {
    if (paramLimits[i].fixed == 1) {
      printf("Fixed parameter detected (i = %d)\n", i);
      lockFlags[i] = 1;
      // cdream stores parameter value for fixed params in lowVals; we'll also
      // store it in highVals, just for consistency
      lowVals[i] = paramsVect[i];
      highVals[i] = paramsVect[i];
    }
    else {
      lockFlags[i] = 0;
      lowVals[i] = paramLimits[i].limits[0];
      highVals[i] = paramLimits[i].limits[1];
    }
    paramNames[i] = theModel->GetParameterName(i);
  }

  // Default is nChains = N_free_params (dimensionality of parameter space)
  if (options->nChains <= 0) {
    printf("Setting nChains = %d (nFreeParams)\n", nFreeParams);
    options->nChains = nFreeParams;
  }

  dream_pars dreamPars;
  SetupDreamParams(&dreamPars, nParamsTot, paramsVect, paramNames, lockFlags, lowVals,
  					highVals);

  // Set up various things in dream_pars struct
  dreamPars.outputRootname = options->outputFileRoot;
  if (options->appendToOutput)
    dreamPars.appendFile = 1;
  dreamPars.numChains = options->nChains;
  dreamPars.maxEvals = options->maxEvals;
  dreamPars.burnIn = options->nBurnIn;
  // the following is tricky, since cdream actually runs a Gelman-Rubin convergence
  // test every dreamPars.loopSteps * dreamPars.gelmanEvals generations
  dreamPars.gelmanEvals = (int)(options->nGelmanEvals / 10);
  dreamPars.scaleReductionCrit = options->GRScaleReductionLimit;
  dreamPars.noise = options->mcmcNoise;
  dreamPars.bstar_zero = options->mcmc_bstar;
  dreamPars.verboseLevel = options->verbose;
  for (int i = 0; i < nParamsTot; i++)
    dreamPars.parameterNames.push_back(paramNames[i]);

  // Construct header using parameter names, etc.
  programHeader.push_back("#\n# Model definition:\n");
//   theModel->PrintModelParamsToStrings(programHeader, paramsVect, paramLimits, NULL,
//   									"#", true);
  theModel->PrintModelParamsToStrings(programHeader, paramsVect, NULL, "#", true);
  programHeader.push_back("#\n# Column Headers\n");
  string headerString = theModel->GetParamHeader();
  headerString += "  likelihood  burn-in    ";
  for (int i = 0; i < dreamPars.nCR - 1; i++)
    headerString += PrintToString("CR%d               ", i + 1);
  headerString += PrintToString("CR%d        ", dreamPars.nCR);
  headerString += "accept\n";
  programHeader.push_back(headerString);
  SetHeaderDreamParams(&dreamPars, programHeader);
  
  dreamPars.fun = &LikelihoodFuncForDREAM;
  // Assign extra "data" that will be passed to likelihood function
  dreamPars.extraData = theModel;

  rng::GSLStream rng;
  if (options->rngSeed > 0)
    rng.alloc(options->rngSeed);
  else
    rng.alloc();   // void alloc(unsigned long seed = time(NULL))


  // OK, now we execute the MCMC process
  // DO THE MCMC!
  printf("\nStart of MCMC processing...\n");
  dream(&dreamPars, &rng);
  printf("\nMCMC chains written to output files %s.0.txt through %s.%d.txt", 
  		options->outputFileRoot.c_str(), options->outputFileRoot.c_str(), options->nChains - 1);


#ifdef USE_PLOG
  LOG_INFO << "*** Shutting down...";
#endif

  // Free up memory
  fftw_free(allPixels);                 // allocated externally, in ReadImageAsVector()
  if (errorPixels_allocated)
    fftw_free(allErrorPixels);          // allocated externally, in ReadImageAsVector()
  if (options->psfImagePresent)
    fftw_free(psfPixels);               // allocated externally, in ReadImageAsVector()
  if (options->psfOversampledImagePresent)
    fftw_free(psfOversampledPixels);    // allocated externally, in ReadImageAsVector()
  if (maskAllocated)
    fftw_free(allMaskPixels);           // allocated externally, in ReadImageAsVector()
  free(paramsVect);
  if (parameterInfo_allocated)
    free(parameterInfo);
  delete theModel;

  FreeVarsDreamParams(&dreamPars);
  delete[] paramNames;

  // Elapsed time reports
  if (options->verbose >= 0) {
    gettimeofday(&timer_end_all, NULL);
    double  microsecs, time_elapsed_all, time_elapsed_fit;
    microsecs = timer_end_all.tv_usec - timer_start_all.tv_usec;
    time_elapsed_all = timer_end_all.tv_sec - timer_start_all.tv_sec + microsecs/1e6;
    printf("\n(Elapsed time: %.6f sec)\n", time_elapsed_all);
  }
  
  printf("\nDone!\n\n");
  
  return 0;
}



void ProcessInput( int argc, char *argv[], MCMCOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();
  string  tempString = "";

  // SET THE USAGE/HELP
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   imfit-mcmc [options] <imagefile.fits>");
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
  optParser->AddUsageLine("     --no-normalize           Do *not* normalize input PSF image");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     (Note that the following 3 options can be specified multiple times)");
  optParser->AddUsageLine("     --overpsf <psf.fits>      Oversampled PSF image to use");
  optParser->AddUsageLine("     --overpsf_scale <n>       Oversampling scale (integer)");
  optParser->AddUsageLine("     --overpsf_region <x1:x2,y1:y2>       Section of image to convolve with oversampled PSF");
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
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" -o  --output <output-root>       root name for output MCMC chain files [default = mcmc_out]");
  optParser->AddUsageLine("     --append                     load state from existing output files and continue from there");
  optParser->AddUsageLine("     --nchains <int>              Number of separate MCMC chains [default = # free parameters in model]");
  optParser->AddUsageLine("     --max-chain-length <int>     Maximum number of likelihood evaluations per chain [default = 100000]");
  optParser->AddUsageLine("     --burnin-length <int>        Number of generations in burn-in phase [default = 5000]");
  optParser->AddUsageLine("     --gelman-evals <int>         Perform Gelman-Rubin convergence check every N generations [default = 5000]");
  optParser->AddUsageLine("     --gelman-rubin-limit <float> Gelman-Rubin scale reduction factor limit [default = 1.01])");
  optParser->AddUsageLine("     --uniform-offset <float>     MCMC uniform-offset term [boundary for uniform offsets of scaling; default = 0.01]");
  optParser->AddUsageLine("     --gaussian-offset <float>    MCMC b^star term [sigma for absolute Gaussian offsets; default = 1.0e-6]");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --quiet                  Turn off printing of updates during the fit");
  optParser->AddUsageLine("     --silent                 Turn off ALL printouts (except fatal errors)");
  optParser->AddUsageLine("     --loud                   Print extra info during the fit");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --max-threads <int>      Maximum number of threads to use");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --seed <int>             RNG seed (for testing purposes)");
  optParser->AddUsageLine("     --no-subsampling         Turn off pixel subsampling near centers of functions");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("EXAMPLES:");
  optParser->AddUsageLine("   imfit-mcmc -c model_config_n100a.dat ngc100.fits -o n100a_mcmc_chain");
  optParser->AddUsageLine("   imfit-mcmc -c model_config_n100b.dat ngc100.fits[405:700,844:1060] --mask ngc100_mask.fits[405:700,844:1060] --gain 4.5 --readnoise 0.7");
  optParser->AddUsageLine("");


  // by default all options are checked on the command line and from option/resource file
  optParser->AddFlag("help", "h");
  optParser->AddFlag("version", "v");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("list-parameters");
  optParser->AddFlag("sample-config");
  optParser->AddFlag("errors-are-variances");
  optParser->AddFlag("errors-are-weights");
  optParser->AddFlag("mask-zero-is-bad");
  optParser->AddFlag("no-normalize");
  optParser->AddFlag("no-subsampling");
  optParser->AddFlag("model-errors");
  optParser->AddFlag("cashstat");
  optParser->AddFlag("poisson-mlr");
  optParser->AddFlag("mlr");
  optParser->AddFlag("quiet");
  optParser->AddFlag("silent");
  optParser->AddFlag("loud");
  optParser->AddOption("noise");
  optParser->AddOption("mask");
  optParser->AddOption("psf");
  optParser->AddOption("overpsf");
  optParser->AddOption("overpsf_scale");
  optParser->AddOption("overpsf_region");
  optParser->AddOption("sky");
  optParser->AddOption("gain");
  optParser->AddOption("readnoise");
  optParser->AddOption("exptime");
  optParser->AddOption("ncombined");
  optParser->AddOption("config", "c");
  optParser->AddOption("output", "o");
  optParser->AddFlag("append");
  optParser->AddOption("nchains");
  optParser->AddOption("max-chain-length");
  optParser->AddOption("burnin-length");
  optParser->AddOption("gelman-evals");
  optParser->AddOption("gelman-rubin-limit");
  optParser->AddOption("uniform-offset");
  optParser->AddOption("gaussian-offset");
  optParser->AddOption("max-threads");
  optParser->AddOption("seed");

  // Comment this out if you want unrecognized (e.g., mis-spelled) flags and options
  // to be ignored only, rather than causing program to exit
  optParser->UnrecognizedAreErrors();
  
  // parse the command line:
  int status = optParser->ParseCommandLine( argc, argv );
  if (status < 0) {
    printf("\nError on command line... quitting...\n\n");
    delete optParser;
    exit(1);
  }


  // Process the results: actual arguments, if any:
  if (optParser->nArguments() > 0) {
    theOptions->imageFileName = optParser->GetArgument(0);
    theOptions->noImage = false;
    printf("\tImage file = %s\n", theOptions->imageFileName.c_str());
  }

  // Process the results: options
  // First five are options which print useful info and then exit the program
  if ( optParser->FlagSet("help") || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("version") ) {
    printf("imfit-mcmc version %s\n\n", VERSION_STRING);
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
  if (optParser->FlagSet("no-normalize")) {
    theOptions->normalizePSF = false;
  }
  if (optParser->FlagSet("no-subsampling")) {
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
  
  // oversampled PSF(s) and region(s)
  if (optParser->OptionSet("overpsf")) {
    string fileName = optParser->GetTargetString("overpsf");
    theOptions->psfOversampledFileNames.push_back(fileName);
    theOptions->psfOversampling = true;
    theOptions->psfOversampledImagePresent = true;
    printf("\tOversampled PSF image = %s\n", fileName.c_str());
  }
  if (optParser->OptionSet("overpsf_scale")) {
    if (NotANumber(optParser->GetTargetString("overpsf_scale").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: overpsf_scale should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    int scale = atoi(optParser->GetTargetString("overpsf_scale").c_str());
    theOptions->psfOversamplingScales.push_back(scale);
    printf("\tPSF oversampling scale = %d\n", scale);
  }
  if (optParser->OptionSet("overpsf_region")) {
  	string psfRegion = optParser->GetTargetString("overpsf_region");
    theOptions->psfOversampleRegions.push_back(psfRegion);
    theOptions->oversampleRegionSet = true;
    theOptions->nOversampleRegions += 1;
    printf("\tPSF oversampling region = %s\n", psfRegion.c_str());
  }

  if (optParser->OptionSet("mask")) {
    theOptions->maskFileName = optParser->GetTargetString("mask");
    theOptions->maskImagePresent = true;
    printf("\tmask image = %s\n", theOptions->maskFileName.c_str());
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
  if (optParser->OptionSet("output")) {
    theOptions->outputFileRoot = optParser->GetTargetString("output");
  }
  if (optParser->FlagSet("append")) {
    printf("\t Current state will be loaded from output files; extended chains will be appended\n");
    theOptions->appendToOutput = true;
  }
  if (optParser->OptionSet("nchains")) {
    if (NotANumber(optParser->GetTargetString("nchains").c_str(), 0, kPosInt)) {
      printf("*** WARNING: number of chains should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nChains = atol(optParser->GetTargetString("nchains").c_str());
    printf("\tNumber of chains = %d\n", theOptions->nChains);
  }
  if (optParser->OptionSet("max-chain-length")) {
    if (NotANumber(optParser->GetTargetString("max-chain-length").c_str(), 0, kPosInt)) {
      printf("*** WARNING: maximum number of evaluations per chain should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->maxEvals = atol(optParser->GetTargetString("max-chain-length").c_str());
    printf("\tMaximum number of likelihood evaluations per chain = %d\n", theOptions->maxEvals);
  }
  if (optParser->OptionSet("burnin-length")) {
    if (NotANumber(optParser->GetTargetString("burnin-length").c_str(), 0, kPosInt)) {
      printf("*** WARNING: number of burn-in evaluations should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nBurnIn = atol(optParser->GetTargetString("burnin-length").c_str());
    printf("\tBurn-in evaluations = %d\n", theOptions->nBurnIn);
  }
  if (optParser->OptionSet("gelman-evals")) {
    if (NotANumber(optParser->GetTargetString("gelman-evals").c_str(), 0, kPosInt)) {
      printf("*** WARNING: Gelman-Rubin generation timing should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nGelmanEvals = atol(optParser->GetTargetString("gelman-evals").c_str());
    printf("\tGelman-Rubin convergence evaluations every %d generations\n", theOptions->nGelmanEvals);
  }
  if (optParser->OptionSet("gelman-rubin-limit")) {
    if (NotANumber(optParser->GetTargetString("gelman-rubin-limit").c_str(), 0, kPosReal)) {
      printf("*** WARNING: Gelman-Rubin scale reduction limit should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->GRScaleReductionLimit = atof(optParser->GetTargetString("gelman-rubin-limit").c_str());
    printf("\tGelman-Rubin scale reduction limit = %f\n", theOptions->GRScaleReductionLimit);
  }
  if (optParser->OptionSet("uniform-offset")) {
    if (NotANumber(optParser->GetTargetString("uniform-offset").c_str(), 0, kPosReal)) {
      printf("*** WARNING: MCMC uniform-offset scale parameter should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->mcmcNoise = atof(optParser->GetTargetString("uniform-offset").c_str());
    printf("\tMCMC uniform-offset parameter = %f\n", theOptions->mcmcNoise);
  }
  if (optParser->OptionSet("gaussian-offset")) {
    if (NotANumber(optParser->GetTargetString("gaussian-offset").c_str(), 0, kPosReal)) {
      printf("*** WARNING: MCMC Gaussian-offset sigma should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->mcmc_bstar = atof(optParser->GetTargetString("gaussian-offset").c_str());
    printf("\tMCMC Gaussian-offset sigma = %f\n", theOptions->mcmc_bstar);
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
bool RequestedFilesPresent( MCMCOptions *theOptions )
{
  bool  allFilesPresent = true;
  
  if (! FileExists(theOptions->configFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find configuration file \"%s\"!\n", 
           theOptions->configFileName.c_str());
    allFilesPresent = false;
  }
  if (! ImageFileExists(theOptions->imageFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find image file \"%s\"!\n", 
           theOptions->imageFileName.c_str());
    allFilesPresent = false;
  }
  if ( (theOptions->maskImagePresent) && (! ImageFileExists(theOptions->maskFileName.c_str())) ) {
    fprintf(stderr, "\n*** ERROR: Unable to find mask file \"%s\"!\n", 
           theOptions->maskFileName.c_str());
    allFilesPresent = false;
  }
  if ( (theOptions->noiseImagePresent) && (! ImageFileExists(theOptions->noiseFileName.c_str())) ) {
    fprintf(stderr, "\n*** ERROR: Unable to find noise-image file \"%s\"!\n", 
           theOptions->noiseFileName.c_str());
    allFilesPresent = false;
  }
  if ( (theOptions->psfImagePresent) && (! ImageFileExists(theOptions->psfFileName.c_str())) ) {
    fprintf(stderr, "\n*** ERROR: Unable to find PSF image file \"%s\"!\n", 
           theOptions->psfFileName.c_str());
    allFilesPresent = false;
  }

  return allFilesPresent;
}



// Note that we only use options from the config file if they have *not*
// already been set by the command line (i.e., command-line options override
// config-file values).
void HandleConfigFileOptions( configOptions *configFileOptions, MCMCOptions *mainOptions )
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




/* ---------------- LikelihoodFuncForDREAM ------------------------------- */

// Function which returns log likelihood given a set of parameters in the
// input parameter state.
// Assumes that extraData is pointer to an instance of ModelObject.
//    chain = which chain this is being called for (ignored)
//    gen = current generation of the chain (ignored)
//    state = current parameter vector
//    extraData = pointer to ModelObject instance
//    recalc -- ignored
double LikelihoodFuncForDREAM( int chain, int gen, const double* state, 
								const void* extraData, bool recalc )
{
  // following is a necessary kludge bcs theModel->GetFitStatistic() won't accept const double*
  double  *params = (double *)state;
  ModelObject *theModel = (ModelObject *)extraData;
  double  chi2;
  
  chi2 = theModel->GetFitStatistic(params);
  return -chi2/2.0;
}




/* ---------------- FUNCTION: MakeMCMCOutputHeader() --------------------- */

void MakeMCMCOutputHeader( vector<string> *headerLines, const string& programName, 
						const int argc, char *argv[] )
{  
  char  *timeStamp = TimeStamp();
  string  tempString;

  tempString = PrintToString("# Markov chain Monte Carlo output from %s\n", programName.c_str());
  headerLines->push_back(tempString);
  tempString = PrintToString("# Generated on %s by the following command:", 
          timeStamp);
  headerLines->push_back(tempString + "\n");
  tempString = "#   ";
  for (int i = 0; i < argc; i++) {
    tempString += PrintToString(" %s", argv[i]);
  }
  headerLines->push_back(tempString + "\n");
}



/* END OF FILE: mcmc_main.cpp -------------------------------------------- */
