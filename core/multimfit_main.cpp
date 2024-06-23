/* FILE: multimfit_main.cpp ---------------------------------------------- */
/*
 * This is the main program file for multimfit.
 *
 *
*/

// Copyright 2009--2024 by Peter Erwin.
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
#include <tuple>
#include <sys/time.h>
#include "fftw3.h"

// logging
//#include "spdlog/spdlog.h"
#include <loguru/loguru.cpp>

#include "definitions.h"
#include "utilities_pub.h"
#include "utilities_multimfit.h"
#include "image_io.h"
#include "getimages.h"
#include "model_object.h"
#include "model_object_multimage.h"
#include "add_functions.h"
#include "param_struct.h"   // for mp_par structure
#include "solver_results.h"
#include "bootstrap_errors.h"
#include "options_base.h"
#include "options_multimfit.h"
#include "psf_oversampling_info.h"
#include "store_psf_oversampling.h"
#include "setup_model_object.h"
#include "imageparams_file_parser.h"
#include "statistics.h"

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
#include "print_results_multi.h"
#include "estimate_memory.h"
#include "sample_configs.h"

using namespace std;


/* ---------------- Quasi-Global Variable Definitions ------------------- */

#ifndef NO_SIGNALS
volatile sig_atomic_t  stopSignal_flag = 0;
#endif


/* ---------------- Definitions & Constants ----------------------------- */

// Option names for use in config files
static string  kGainString = "GAIN";
static string  kReadNoiseString = "READNOISE";
static string  kExpTimeString = "EXPTIME";
static string  kNCombinedString = "NCOMBINED";
static string  kOriginalSkyString = "ORIGINAL_SKY";

const char *  LOG_FILENAME = "log_multimfit.txt";
// const char *  LOG_FILENAME_BASE = "log_multimfit";

#ifdef USE_OPENMP
#define VERSION_STRING      "0.8 (OpenMP-enabled)"
#else
#define VERSION_STRING      "0.8"
#endif



/* ------------------- Function Prototypes ----------------------------- */
void ProcessInput( int argc, char *argv[], shared_ptr<MultimfitOptions> theOptions );
bool RequestedFilesPresent( vector<ImageInfo> &imageInfoVect, int nDataImages );
void signal_handler( int signal );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */





/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  long  nPixels_tot;
  int  nRows_psf = 0;
  int  nColumns_psf = 0;
  long  nDegFreedom;
  int  nDataImages;
  int  nParamsTot, nFreeParams;
  int  nGlobalParams = 0;
  vector<double>  localParamsList;
  int  nLocalParams = 0;
  double  *dataPixels;
  double  *psfPixels;
  double  *errorPixels;
  double  *maskPixels;
  vector<PsfOversamplingInfo *>  psfOversamplingInfoVect;
  vector< vector<PsfOversamplingInfo *> > allPsfOversamplingInfoVectors;

  vector<double *> dataPixelsUsed;   // so we know how many allocated pixel-vectors to free
  vector<double *> maskPixelsUsed;   // so we know how many allocated pixel-vectors to free
  vector<double *> errorPixelsUsed;  // so we know how many allocated pixel-vectors to free
  vector<double *> psfPixelsUsed;    // so we know how many allocated pixel-vectors to free
  vector<double *> psfOversampledPixelsUsed;  // so we know how many allocated pixel-vectors to free

  double  *paramsVect;
  string  noiseImage;
  vector<string>  functionList;
  vector<string>  functionLabelList;
  vector<double>  parameterList;
  mp_par  newParamLimit;
  vector<mp_par>  paramLimits;
  vector<mp_par>  imageDescParamLimits;
  vector<mp_par>  paramLimits_local;
  vector<int>  functionSetIndices;
  vector< map<string, string> > optionalParamsMap;
  bool  paramLimitsExist_func = false;
  bool  paramLimitsExist_imageDesc = false;
  bool  paramLimitsExist_local = false;
  bool  paramLimitsExist = false;
  vector<mp_par>  parameterInfo;
  int  status, fitStatus, nSucessfulIterations;
  SolverResults  resultsFromSolver;
  vector<string>  imageCommentsList;
  vector<int> nColumnsRowsVect;
  vector<int> nColumnsVect, nRowsVect, nPsfColumnsVect, nPsfRowsVect;

  shared_ptr<MultimfitOptions> options;
  ModelObjectMultImage  *theMultImageModel;
  int  nColumns_thisImage, nRows_thisImage;
  vector<ImageInfo> imageInfoVect;
  int  nGlobalFuncs = 0;
  double  *imageDescParamsVect;
  int  nImageDescParams;
  configOptions  userConfigOptions;
  string  progNameVersion = "multimfit ";
  vector<string> programHeader;
  bool  userInterrupted = false;
  FILE  *bootstrapSaveFile_ptr = nullptr;
  bool  didBootstrap = false;
  // timing-related
  struct timeval  timer_start_all, timer_end_all;
  struct timeval  timer_start_fit, timer_end_fit;
  struct timeval  timer_start_bootstrap, timer_end_bootstrap;

  gettimeofday(&timer_start_all, NULL);

  progNameVersion += VERSION_STRING;
  MakeOutputHeader(&programHeader, progNameVersion, argc, argv);

 
  // ** Define default options, then process the command line
  options = make_shared<MultimfitOptions>();
  ProcessInput(argc, argv, options);

  if (options->loggingOn) {
    // turn off writing log output to stderr
    loguru::g_stderr_verbosity = loguru::Verbosity_OFF;
    // initialize logging
    loguru::init(argc, argv);
    loguru::add_file(LOG_FILENAME, loguru::Append, loguru::Verbosity_MAX);
    LOG_F(INFO, "*** Starting up...");
  }


  // MULTIMFIT: Read configuration file to get model
  if (! FileExists(options->configFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find configuration file \"%s\"!\n\n", 
           options->configFileName.c_str());
    if (options->loggingOn)
      LOG_F(ERROR, "Unable to find configuration file \"%s\"", options->configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options->configFileName, true, functionList, functionLabelList,
  							parameterList, paramLimits, functionSetIndices, 
  							paramLimitsExist_func, userConfigOptions);
  if (status != 0) {
    fprintf(stderr, "\n*** ERROR: Failure reading configuration file \"%s\"!\n\n", 
    			options->configFileName.c_str());
    if (options->loggingOn)
      LOG_F(ERROR, "Problem reading configuration file \"%s\"", options->configFileName.c_str());
    return -1;
  }
  nGlobalFuncs = functionList.size();


  // MULTIMFIT: Read image-info file to get info about individual images
  if (options->imageInfoFile_exists) {
    if (! FileExists(options->imageInfoFile.c_str())) {
      fprintf(stderr, "\n*** ERROR: Unable to find image-info file \"%s\"!\n\n", 
             options->imageInfoFile.c_str());
      if (options->loggingOn)
        LOG_F(ERROR, "Unable to find image-info file \"%s\"", options->imageInfoFile.c_str());
      return -1;
    }
    status = ReadImageParamsFile(options->imageInfoFile, nDataImages, imageInfoVect);
    if (status < 0) {
      fprintf(stderr, "\n*** ERROR: Problems parsing image-info file \"%s\"!\n\n", 
             options->imageInfoFile.c_str());
      if (options->loggingOn)
        LOG_F(ERROR, "Problem parsing image-info file \"%s\"", options->imageInfoFile.c_str());
      return -1;
    }
    // FIXME: Determine whether we have parameter limits in the image-description
    // parameters, to make sure that paramLimitsExist is = true even if the *model* parameters
    // had no limits
    std::tie(imageDescParamsVect, nImageDescParams) = GetImageDescriptionParams(imageInfoVect, 
    												imageDescParamLimits, paramLimitsExist_imageDesc);
    printf("Image-description params from image-info file \"%s\": ", options->imageInfoFile.c_str());
    for (int i = 0; i < nImageDescParams; i++)
      printf("%.1f ", imageDescParamsVect[i]);
    printf("\n");
  }
  else {
    fprintf(stderr, "\n*** ERROR: No image-info file was supplied!\n\n");
    if (options->loggingOn)
      LOG_F(ERROR, "No image-info file supplied!");
    return -1;
  }


  // MULTIMFIT: Check for valid info from image-params file
  //    Check for whether data-image filenames were supplied, etc.
  //    Check whether # data-image filenames = # image specifications
  //    Check for existence of files
  for (int i = 0; i < nDataImages; i++) {
    if (! imageInfoVect[i].dataImage_present) {
      fprintf(stderr, "\n*** ERROR: No data image filename was specified for image #%d!\n\n", 
             i + 1);
      if (options->loggingOn)
        LOG_F(ERROR, "No data image filename for image %d", i + 1);
      return -1;
    }
  }
  
  // Check for presence of user-requested files; if any are missing, quit.
  // (Appropriate error messages regarding which files are missing will be printed
  // to stderr by RequestedFilesPresent)
  if (! RequestedFilesPresent(imageInfoVect, nDataImages)) {
    fprintf(stderr, "\n");
    return -1;
  }
  

  // MULTIMFIT: Creation and initial setup of ModelObjectMultImage
  if (options->loggingOn)
    LOG_F(INFO, "creating ModelObjectMultImage instance...");
  theMultImageModel = new ModelObjectMultImage();
  theMultImageModel->SetupGeneralOptions(options);
  
  // MULTIMFIT: loop over creation of ModelObject instances, one per image
  //    Create new ModelObject instance via SetupModelObject
  //    Call AddFunctions using new ModelObject instance
  //    Add new ModelObject instance to ModelObjectMultImage
  // After loop finishes, call AddFunctions using ModelObjectMultImage
  
  ModelObject * newModel;
  if (options->loggingOn)
    LOG_F(INFO, "Starting loop setting up and adding individual ModelObject instances...");
  for (int i = 0; i < nDataImages; i++) {
    // reset things to default values
    options->maskImagePresent = false;
    options->noiseImagePresent = false;
    options->psfImagePresent = false;
    options->psfOversampledImagePresent = false;
    options->oversampleRegionSet = false;
    nColumnsRowsVect.clear();
    psfOversamplingInfoVect.clear();
    // reset these to point to nullptr, otherwise we'll get e.g. psf pixels from
    // previous rounds passed in to SetupModelObject
    maskPixels = nullptr;
    errorPixels = nullptr;
    psfPixels = nullptr;

    // image characteristics
    // options->gain, options->readNoise, options->expTime, options->nCombined, options->originalSky
    options->readNoise = imageInfoVect[i].readNoise;
    options->gain = imageInfoVect[i].gain;
    options->originalSky = imageInfoVect[i].originalSky;
    options->expTime = imageInfoVect[i].expTime;
    options->nCombined = imageInfoVect[i].nCombined;

    // Get data image and sizes
    printf("Reading data image (\"%s\") ...\n", imageInfoVect[i].dataImageFileName.c_str());
    if (options->loggingOn)
      LOG_F(INFO, "   i = %d: reading data image \"%s\"", i, imageInfoVect[i].dataImageFileName.c_str());
    dataPixels = ReadImageAsVector(imageInfoVect[i].dataImageFileName, 
    								&nColumns_thisImage, &nRows_thisImage);
    if (dataPixels == nullptr) {
      fprintf(stderr,  "\n*** ERROR: Unable to read image file \"%s\"!\n\n", 
    			imageInfoVect[i].dataImageFileName.c_str());
      if (options->loggingOn)
        LOG_F(ERROR, "   Unable to read data image file \"%s\"", options->imageInfoFile.c_str());
      return -1;
    }
    nPixels_tot = (long)nColumns_thisImage * (long)nRows_thisImage;
    printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns_thisImage, nRows_thisImage, nPixels_tot);
    // Determine X0,Y0 pixel offset values if user specified an image section
    int  X0_offset_current, Y0_offset_current;
    std::tie(X0_offset_current, Y0_offset_current) = DetermineImageOffset(imageInfoVect[i].dataImageFileName);

    // Get (and check) mask and/or error images
    if ((imageInfoVect[i].maskImage_present) || (imageInfoVect[i].errorImage_present)) {
      bool maskAllocated = false, errorPixels_allocated = false;
      std::tie(maskPixels, errorPixels, status) = GetMaskAndErrorImages(nColumns_thisImage, 
      							nRows_thisImage, imageInfoVect[i].maskImageFileName, 
      							imageInfoVect[i].errorImageFileName, maskAllocated, errorPixels_allocated);
      if (status < 0)
        return -1;
      if (maskAllocated) {
        options->maskImagePresent = true;
        maskPixelsUsed.push_back(maskPixels);
        options->maskFormat = imageInfoVect[i].maskFormat;
      }
      if (errorPixels_allocated) {
        options->noiseImagePresent = true;
        errorPixelsUsed.push_back(errorPixels);
        options->errorType = imageInfoVect[i].errorType;
      }
    }
    
    // Read in PSF image, if supplied
    nColumns_psf = nRows_psf = 0;
    if (imageInfoVect[i].psfImage_present) {
      if (options->loggingOn)
        LOG_F(INFO, "          reading PSF image \"%s\"", imageInfoVect[i].psfImageFileName.c_str());
      std::tie(psfPixels, nColumns_psf, nRows_psf, status) = GetPsfImage(imageInfoVect[i].psfImageFileName);
      if (status < 0)
        return -1;
      options->psfImagePresent = true;
      psfPixelsUsed.push_back(psfPixels);
    }

    // Read in oversampled PSF image, if supplied
    if (imageInfoVect[i].psfOversampledImage_present) {
      options->psfOversampling = true;
      status = ExtractAndStorePsfOversampling(imageInfoVect[i], i, psfOversamplingInfoVect);
      if (status < 0)
        exit(status);
    }
    // store copy of current psfOversamplingInfoVect so we can properly de-allocate 
    // its memory later (also useful for estimating memory usage)
    allPsfOversamplingInfoVectors.push_back(psfOversamplingInfoVect);

    nColumnsRowsVect.push_back(nColumns_thisImage);
    nColumnsRowsVect.push_back(nRows_thisImage);
    nColumnsRowsVect.push_back(nColumns_psf);
    nColumnsRowsVect.push_back(nRows_psf);
    // store copies of # columns, rows, etc. for later use in estimating memory usage
    nColumnsVect.push_back(nColumns_thisImage);
    nRowsVect.push_back(nRows_thisImage);
    nPsfColumnsVect.push_back(nColumns_psf);
    nPsfRowsVect.push_back(nRows_psf);
    // generate and set up a new ModelObject instance!
    newModel = SetupModelObject(options, nColumnsRowsVect, dataPixels, psfPixels, 
    						maskPixels, errorPixels, psfOversamplingInfoVect);
    newModel->AddImageOffsets(X0_offset_current, Y0_offset_current);
    newModel->AddDataFilename(imageInfoVect[i].dataImageFileName);

    // Handle possible local, per-image functions
    vector<string> currentFunctionList = functionList;
    vector<string> currentFunctionLabelList = functionLabelList;
    vector<int> currentFunctionSetIndices = functionSetIndices;
    vector< map<string, string> > currentOptionalParamsMap = optionalParamsMap;
    vector<bool> globalFuncFlags(currentFunctionList.size(), true);
    if (imageInfoVect[i].perImageFunctionsExist) {
      int  nLocalFuncs = imageInfoVect[i].perImageFuncNames.size();
      printf("Image %d: %d local functions found!\n", i, nLocalFuncs);
      for (int j = 0; j < nLocalFuncs; j++) {
        currentFunctionList.push_back(imageInfoVect[i].perImageFuncNames[j]);
        currentFunctionLabelList.push_back(imageInfoVect[i].perImageFuncLabels[j]);
        printf("newFunctionSetIndex = nGlobalFuncs + imageInfoVect[i].perImageFuncSetIndices[j]");
        printf("  %d + %d\n", nGlobalFuncs, imageInfoVect[i].perImageFuncSetIndices[j]);
        int  newFunctionSetIndex = nGlobalFuncs + imageInfoVect[i].perImageFuncSetIndices[j];
        currentFunctionSetIndices.push_back(newFunctionSetIndex);
        //FIXME: handle any extra, local-function optional-parameter-map values
        globalFuncFlags.push_back(false);
      }
    }
#ifdef USE_LOGGING
    if (options->loggingOn) {
      LOG_F(INFO, "Calling AddFunctions...");
      LOG_F(INFO, "   function (label): globalFuncFlags");
      for (int j = 0; j < (int)currentFunctionList.size(); j++) {
        string  bool_str;
        if (globalFuncFlags[j])
          bool_str = "true";
        else
          bool_str = "false";
        LOG_F(INFO, "   %s (%s): %s", currentFunctionList[j].c_str(),
        		currentFunctionLabelList[j].c_str(), bool_str.c_str());
      }
      for (int j = 0; j < (int)currentFunctionSetIndices.size(); j++) {
        LOG_F(INFO, "  currentFunctionSetIndices[%d] = %d", j, currentFunctionSetIndices[j]);
      }
    }
#endif

    // Add functions to the model object; also tells model object where function sets start
    status = AddFunctions(newModel, currentFunctionList, currentFunctionLabelList, 
    					currentFunctionSetIndices, options->subsamplingFlag, 0, 
    					currentOptionalParamsMap, globalFuncFlags);
    if (status < 0) {
  	  fprintf(stderr, "*** ERROR: Failure in AddFunctions!\n\n");
  	  exit(-1);
    }
    newModel->FinalModelSetup();
    // Note that ModelObjectMultImage will handle de-allocation of ModelObject instance
    theMultImageModel->AddModelObject(newModel);
  }
  if (options->loggingOn)
    LOG_F(INFO, "Done with adding ModelObject instances.");

  printf("main: theMultImageModel has %d data images (ModelObject instances)\n",
  		theMultImageModel->GetNImages());
  theMultImageModel->PrintDescription();
  
  // Finally, add functions to ModelObjectMultImage object (important!)
  if (options->loggingOn)
    LOG_F(INFO, "Adding functions to theMultImageModel...");
  // NOTE: this does *not* add any of the per-image functions!
  status = AddFunctions(theMultImageModel, functionList, functionLabelList, functionSetIndices, 
  						options->subsamplingFlag, 0, optionalParamsMap);
  if (status < 0) {
  	fprintf(stderr, "*** ERROR: Failure in AddFunctions!\n\n");
    if (options->loggingOn)
  	  LOG_F(ERROR, "Failure in AddFunctions!");
  	return -1;
  }
  theMultImageModel->FinalModelSetup();
  
  // Note that SetupModelObject has already told the individual ModelObject instances
  // which fit statistic we're using. The following is to alert the ModelObjectMultImage
  // instance. (If the default of chi^2 is being used, there's nothing to do here.)
  if (options->useCashStatistic) {
    if ((options->solver == MPFIT_SOLVER) && (! options->printFitStatisticOnly)) {
      fprintf(stderr, "*** ERROR -- Cash statistic cannot be used with L-M solver!\n\n");
      if (options->loggingOn)
        LOG_F(ERROR, "Cash statistic cannot be used with L-M solver!");
      return -1;
    }
    status = theMultImageModel->UseCashStatistic();
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObjectMultImage::UseCashStatistic!\n\n");
      if (options->loggingOn)
        LOG_F(ERROR, "Failure in ModelObjectMultImage::UseCashStatistic!");
      return -1;
    }
  } 
  else if (options->usePoissonMLR) {
    theMultImageModel->UsePoissonMLR();
  }


  // Set up parameter vector(s), now that we know how many total parameters there are
  nParamsTot = nFreeParams = theMultImageModel->GetNParams();
  nGlobalParams = (int)parameterList.size();
  for (int i = 0; i < nDataImages; i++) {
    int  nLocalParams_thisImage = (int)imageInfoVect[i].perImageParamVals.size();
    nLocalParams += nLocalParams_thisImage;
    if (nLocalParams_thisImage > 0)
      localParamsList.insert(localParamsList.end(), imageInfoVect[i].perImageParamVals.begin(), 
      							imageInfoVect[i].perImageParamVals.end());
  }

  printf("(global) parameterList: ");
  for (int i = 0; i < nGlobalParams; i++)
    printf(" %.2f", parameterList[i]);
  int  nInputParams_nonLocal = nImageDescParams + nGlobalParams;
  int  nInputParamsTot = nInputParams_nonLocal + nLocalParams;
  printf("\n%d total parameters in ModelObjectMultImage\n", nParamsTot);
  if (nParamsTot != nInputParamsTot) {
  	fprintf(stderr, "*** ERROR: number of (image-desc. + global) input parameters (%d) does not equal", 
  	       nInputParamsTot);
  	fprintf(stderr, " required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	return -1;
  }

  theMultImageModel->PrintDescription();


  // Final fitting-oriented setup for ModelObject instance (generates data-based error
  // vector if needed, created final weight vector from mask and optionally from
  // error vector)
  if (options->loggingOn)
    LOG_F(INFO, "Calling theMultImageModel->FinalSetupForFitting()...");
  status = theMultImageModel->FinalSetupForFitting();
  if (status < 0) {
    fprintf(stderr, "*** ERROR: Failure in ModelObjectMultImage::FinalSetupForFitting!\n\n");
    if (options->loggingOn)
      LOG_F(ERROR, "Failure in ModelObjectMultImage::FinalSetupForFitting!");
    return -1;
  }
  
  
  /* START OF MINIMIZATION-ROUTINE-RELATED CODE */
  // Parameter limits and other info:
  // First we create a vector of mp_par structures, containing parameter constraints
  // (if any) *and* any other useful info (like X0,Y0 offset values). We populate
  // this with copies of the mp_par structures from imageDescParamLimits and paramLimits.
  
  printf("Setting up parameter information vector ...\n");
  if (options->loggingOn) {
    LOG_F(INFO, "Setting up parameterInfo vector...");
    LOG_F(INFO, "Populating parameterInfo from imageDescParamLimits...");
  }
  parameterInfo.clear();
  
  // Reminder
  //    nImageDescParams = number of image-description parameters
  //    nGlobalParams = number of X0,Y0 and function parameters for global model
  //    nLocalParams = number of parameters for local, per-image functions (can be 0)
  //
  //    nInputParams_nonLocal = nImageDescParams + nGlobalParams;
  //    nParamsTot = nInputParams_tot = sum of all above
  
  // Parameter limits
  // Add parameter limits, etc., for image-description parameters
  for (int i = 0; i < nImageDescParams; i++) {
  	bzero(&newParamLimit, sizeof(mp_par));
    parameterInfo.push_back(newParamLimit);
    parameterInfo[i].fixed = imageDescParamLimits[i].fixed;
    parameterInfo[i].limited[0] = imageDescParamLimits[i].limited[0];
    parameterInfo[i].limited[1] = imageDescParamLimits[i].limited[1];
    parameterInfo[i].limits[0] = imageDescParamLimits[i].limits[0];
    parameterInfo[i].limits[1] = imageDescParamLimits[i].limits[1];
  }
  if (options->loggingOn)
    LOG_F(INFO, "Populating parameterInfo from paramLimits...");
  // Add parameter limits, etc., for model parameters
  for (int i = nImageDescParams; i < nInputParams_nonLocal; i++) {
    // counting from i = nImageDescParams to nParamsTot, but paramLimits only
    // runs from 0 to (nParamsTot - nImageDescParams)
  	bzero(&newParamLimit, sizeof(mp_par));
    parameterInfo.push_back(newParamLimit);
    parameterInfo[i].fixed = paramLimits[i - nImageDescParams].fixed;
    parameterInfo[i].limited[0] = paramLimits[i - nImageDescParams].limited[0];
    parameterInfo[i].limited[1] = paramLimits[i - nImageDescParams].limited[1];
    parameterInfo[i].limits[0] = paramLimits[i - nImageDescParams].limits[0];
    parameterInfo[i].limits[1] = paramLimits[i - nImageDescParams].limits[1];
  }
  // Finally, add parameter limits, etc., for local, per-image function parameters
  if (nLocalParams > 0) {
    // First, merge together all the local, per-image limits
    for (int n = 0; n < nDataImages; n++) {
      int nLocalParams_thisImage = (int)imageInfoVect[n].perImageParamVals.size();
      printf("   image %d per-image parameter limits:\n", n);
      for (int j = 0; j < nLocalParams_thisImage; j++) {
	    bzero(&newParamLimit, sizeof(mp_par));
	    int  fixed = imageInfoVect[n].perImageParamLimits[j].fixed;
	    if (fixed)
	      paramLimitsExist_local = true;
	    newParamLimit.fixed = fixed;
	    int  limited_0 = imageInfoVect[n].perImageParamLimits[j].limited[0];
	    int  limited_1 = imageInfoVect[n].perImageParamLimits[j].limited[1];
	    if ( (limited_0) || (limited_1) )
	      paramLimitsExist_local = true;
	    newParamLimit.limited[0] = limited_0;
	    newParamLimit.limited[1] = limited_1;
	    newParamLimit.limits[0] = imageInfoVect[n].perImageParamLimits[j].limits[0];
	    newParamLimit.limits[1] = imageInfoVect[n].perImageParamLimits[j].limits[1];
        parameterInfo.push_back(newParamLimit);
      }
    }
  }

  for (int i = 0; i < nParamsTot; i++)
    if (parameterInfo[i].fixed == 1)
      nFreeParams--;
  nDegFreedom = theMultImageModel->GetNValidPixels() - nFreeParams;
  printf("%d free parameters (%ld degrees of freedom)\n", nFreeParams, nDegFreedom);



  // Now that we know all about the model (including nFreeParams), estimate the
  // memory usage and warn if it will be large
  long  estimatedMemory = 0;
  double  nGBytes;
  bool  usingLevMar, usingCashTerms;
  if (options->solver == MPFIT_SOLVER)
    usingLevMar = true;
  else
    usingLevMar = false;
  if ((options->useCashStatistic) || (options->usePoissonMLR))
    usingCashTerms = true;
  else
    usingCashTerms = false;
  for (int i = 0; i < nDataImages; i++) {
    estimatedMemory += EstimateMemoryUse(nColumnsVect[i], nRowsVect[i], nPsfColumnsVect[i], 
    									nPsfRowsVect[i], nFreeParams, usingLevMar, 
    									usingCashTerms, options->saveResidualImage, 
  										options->saveModel);
  	psfOversamplingInfoVect = allPsfOversamplingInfoVectors[i];
  	if (psfOversamplingInfoVect.size() > 0)
      estimatedMemory += EstimatePsfOversamplingMemoryUse(psfOversamplingInfoVect);
  }

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
  if (options->loggingOn)
    LOG_F(INFO, "Copying initial image-description parameter values into paramsVect...");
  for (int i = 0; i < nImageDescParams; i++)
    paramsVect[i] = imageDescParamsVect[i];
  if (options->loggingOn)
    LOG_F(INFO, "Copying initial model parameter values into paramsVect...");
  for (int i = nImageDescParams; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i - nImageDescParams];
  for (int i = nInputParams_nonLocal; i < nInputParamsTot; i++)
    paramsVect[i] = localParamsList[i - nInputParams_nonLocal];
  // Correct for image subsection (X0, Y0 offsets), if any
  CorrectForImageOffsets(paramsVect, parameterInfo, theMultImageModel);

  // tell ModelObject about parameterInfo (mainly useful for printing-related methods)
  theMultImageModel->AddParameterInfo(parameterInfo);


  // Debugging printout code to check that parameter & parameter-limit vectors
  // are set up OK
  if (options->loggingOn) {
    LOG_F(INFO, "\nInitial parameter values and limits:\n");
    LOG_F(INFO, "   i     param   fixed  limited    lower,upper\n");
    for (int i = 0; i < nParamsTot; i++) {
      LOG_F(INFO, "  %2d: %8.3f     %d  %d,%d  %8.3f,%8.3f\n",
    		i, paramsVect[i], parameterInfo[i].fixed, parameterInfo[i].limited[0], 
    		parameterInfo[i].limited[1], parameterInfo[i].limits[0], parameterInfo[i].limits[1]);
    }
  }

  
  // OK, now we either print chi^2 value for the input parameters and quit, or
  // else call one of the solvers!
  if (options->printFitStatisticOnly) {
    if (options->loggingOn)
      LOG_F(INFO, "Calling PrintFitStatistic...");
    printf("\n");
    PrintFitStatistic(paramsVect, theMultImageModel, nFreeParams);
    printf("\n");
    options->saveBestFitParams = false;
  }
  else {
    // DO THE FIT!
    printf("\nPerforming fit by minimizing ");
    if (options->useCashStatistic)
      printf("Cash statistic:\n");
    else if (options->usePoissonMLR)
      printf("Poisson MLR statistic:\n");
    else
      printf("chi^2:\n");
    gettimeofday(&timer_start_fit, NULL);
    // Set signal-handling so Ctrl-C (SIGINT) is intercepted
    signal(SIGINT, signal_handler);
    if (options->loggingOn)
      LOG_F(INFO, "* Starting fit...");
    nPixels_tot = theMultImageModel->GetNDataValues();
    if ((paramLimitsExist_func) || (paramLimitsExist_imageDesc) || (paramLimitsExist_local))
      paramLimitsExist = true;
    if (options->loggingOn) {
      LOG_F(5, "paramLimitsExist_func = %d", paramLimitsExist_func);
      LOG_F(5, "paramLimitsExist_imageDesc = %d", paramLimitsExist_imageDesc);
      LOG_F(5, "paramLimitsExist_local = %d", paramLimitsExist_local);
      LOG_F(5, "paramLimitsExist = %d", paramLimitsExist);
    }
    fitStatus = DispatchToSolver(options->solver, nParamsTot, nFreeParams, nPixels_tot, 
    						paramsVect, parameterInfo, theMultImageModel, options->ftol, 
    						paramLimitsExist, options->verbose, &resultsFromSolver, 
    						options->nloptSolverName, options->rngSeed);
    gettimeofday(&timer_end_fit, NULL);
#ifndef NO_SIGNALS
    if (stopSignal_flag == 1)
      userInterrupted = true;
#endif
    if (options->loggingOn)
      LOG_F(INFO, "* Done with fit.");
    							
    PrintResults(paramsVect, theMultImageModel, nFreeParams, fitStatus, resultsFromSolver);
    string  outputSummaryFile = options->outputParameterRootName + "_summary.dat";
    printf("\nSaving summary of fit in %s...\n", outputSummaryFile.c_str());
    SaveParameters(paramsVect, theMultImageModel, outputSummaryFile, programHeader,
    			nFreeParams, options->solver, fitStatus, resultsFromSolver);
  }



  printf("Saving single-image best-fit parameter files (root name = \"%s\"):\n",
  		options->outputParameterRootName.c_str());
  SaveMultImageParameters(paramsVect, theMultImageModel, options->outputParameterRootName,
  							programHeader);

  string  outputImageInfoFilename = options->outputParameterRootName + "_imageinfo.dat";
  printf("Saving best-fit image-info file \"%s\"...\n", outputImageInfoFilename.c_str());
  SaveImageInfoParameters(paramsVect, theMultImageModel, imageInfoVect, resultsFromSolver,
  							outputImageInfoFilename, programHeader);


  // Optional bootstrap resampling
  if ((options->doBootstrap) && (options->bootstrapIterations > 0) && (! userInterrupted)) {
    if (options->saveBootstrap) {
      bootstrapSaveFile_ptr = fopen(options->outputBootstrapFileName.c_str(), "w");
      // write general info + best-fitting params as a commented-out header
      // SaveParameters2(bootstrapSaveFile_ptr, paramsVect, theModel, programHeader, "#");
    }
    
    printf("\nNow doing bootstrap resampling (%d iterations) to estimate errors...\n",
           options->bootstrapIterations);
    gettimeofday(&timer_start_bootstrap, NULL);
    nSucessfulIterations = BootstrapErrors(paramsVect, parameterInfo, paramLimitsExist, 
    									theMultImageModel, options->ftol, options->bootstrapIterations, 
    									nFreeParams, theMultImageModel->WhichFitStatistic(), 
    									bootstrapSaveFile_ptr, options->rngSeed, true);
    gettimeofday(&timer_end_bootstrap, NULL);
    if (options->saveBootstrap) {
      if (nSucessfulIterations > 0)
        printf("\nBootstrap-resampling output saved to file \"%s\".\n", options->outputBootstrapFileName.c_str());
      // even if we didn't have any successful iterations, we need to close the file
      fclose(bootstrapSaveFile_ptr);
    }
    didBootstrap = true;
  }


  if (options->loggingOn) {
    LOG_F(INFO, "*** Shutting down...");
    printf("\nInternal logging output saved to %s\n\n", LOG_FILENAME);
  }


  // Free up memory
  if (dataPixelsUsed.size() > 0) {
    for (int i = 0; i < (int)dataPixelsUsed.size(); i++)
      fftw_free(dataPixelsUsed[i]);      // allocated in ReadImageAsVector()
    dataPixelsUsed.clear();
  }
  if (maskPixelsUsed.size() > 0) {
    for (int i = 0; i < (int)maskPixelsUsed.size(); i++)
      fftw_free(maskPixelsUsed[i]);      // allocated in ReadImageAsVector()
    maskPixelsUsed.clear();
  }
  if (errorPixelsUsed.size() > 0) {
    for (int i = 0; i < (int)errorPixelsUsed.size(); i++)
      fftw_free(errorPixelsUsed[i]);      // allocated in ReadImageAsVector()
    errorPixelsUsed.clear();
  }
  if (psfPixelsUsed.size() > 0) {
    for (int i = 0; i < (int)psfPixelsUsed.size(); i++)
      fftw_free(psfPixelsUsed[i]);      // allocated in ReadImageAsVector()
    psfPixelsUsed.clear();
  }
  if (allPsfOversamplingInfoVectors.size() > 0) {
    for (int i = 0; i < (int)allPsfOversamplingInfoVectors.size(); i++) {
      vector<PsfOversamplingInfo *>  oversamplingInfoVect = allPsfOversamplingInfoVectors[i];
      if (oversamplingInfoVect.size() > 0) {
        for (int nn = 0; nn < (int)oversamplingInfoVect.size(); nn++)
          free(oversamplingInfoVect[nn]);
        oversamplingInfoVect.clear();
      }
    }
  }


  // Elapsed time reports
  if (options->verbose >= 0) {
    gettimeofday(&timer_end_all, NULL);
    double  microsecs, time_elapsed_all, time_elapsed_fit, time_elapsed_bootstrap;
    microsecs = timer_end_all.tv_usec - timer_start_all.tv_usec;
    time_elapsed_all = timer_end_all.tv_sec - timer_start_all.tv_sec + microsecs/1e6;
    if (options->printFitStatisticOnly)
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



  // FIXME: handle possible outputs of model and residual images
  // Handle assorted output requests
  // Note that from this point on, we handle failures reported by SaveVectorAsImage as
  // "warnings" and don't immediately exit, since we're close to the end of the program
  // anyway, and the user might just have given us a bad path for one of the output images
  if (options->saveModel) {
    theMultImageModel->CreateAllModelImages(paramsVect);
    for (int i = 0; i < nDataImages; i++) {
      std::tie(nColumns_thisImage, nRows_thisImage) = theMultImageModel->GetImageDimensions(i);
      string multImageOutputName = PrintToString("%s_%d.fits", options->outputModelFileName.c_str(), i);
      PrepareImageComments(&imageCommentsList, progNameVersion, multImageOutputName,
    					  imageInfoVect[i].psfImage_present, imageInfoVect[i].psfImageFileName, 
    					  HDR_MODELIMAGE, imageInfoVect[i].dataImageFileName);
      printf("Saving model image in file \"%s\"\n", multImageOutputName.c_str());
      status = SaveVectorAsImage(theMultImageModel->GetModelImageVector(i), multImageOutputName, 
                        nColumns_thisImage, nRows_thisImage, imageCommentsList);
      if (status != 0) {
        fprintf(stderr, "\n*** WARNING: Failure saving model-image file \"%s\"!\n\n",
        				multImageOutputName.c_str());
      }
    }
  }


  free(paramsVect);
  free(imageDescParamsVect);               // allocated externally, in GetImageDescriptionParams()
  delete theMultImageModel;
//   delete options;
 
  printf("Done!\n\n");

  return 0;
}



void ProcessInput( int argc, char *argv[], shared_ptr<MultimfitOptions> theOptions )
{

  CLineParser *optParser = new CLineParser();
  string  tempString = "";

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   multimfit [options]");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine(" -v  --version                Prints version number");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --list-parameters        Prints list of parameter names for each available function");
  tempString = PrintToString("     --sample-config          Generates an example configuration file (%s)", configImfitFile.c_str());
  optParser->AddUsageLine(tempString);
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" -c  --config <config-file>   configuration file [REQUIRED!]");
  optParser->AddUsageLine(" -i  --image-info <image-info-file> image-description file [REQUIRED!]");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --save-params-root <root-name>       Root name for best-fit parameters outputs [default = bestfit_parameters_multimfit]");
  optParser->AddUsageLine("     --save-model <root-name>             Save best-fit model image");
  optParser->AddUsageLine("     --save-residual <outputname.fits>    Save residual (data - best-fit model) image");
  optParser->AddUsageLine("     --save-weights <outputname.fits>     Save weight image");
  optParser->AddUsageLine("");
  // the following options *may* need to be removed
  optParser->AddUsageLine("     --errors-are-variances   Indicates that values in noise image = variances (instead of sigmas)");
  optParser->AddUsageLine("     --errors-are-weights     Indicates that values in noise image = weights (instead of sigmas)");
  optParser->AddUsageLine("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --model-errors           Use model values (instead of data) to estimate errors for chi^2 computation");
  optParser->AddUsageLine("     --cashstat               Use Cash statistic instead of chi^2");
  optParser->AddUsageLine("     --poisson-mlr            Use Poisson maximum-likelihood-ratio statistic instead of chi^2");
  optParser->AddUsageLine("     --mlr                    Same as --poisson-mlr");
  // end of possible elimination zone
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
  optParser->AddUsageLine("     --logging                Save logging outputs to file");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --max-threads <int>      Maximum number of threads to use");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --seed <int>             RNG seed (for testing purposes)");
  optParser->AddUsageLine("     --nosubsampling          Turn off pixel subsampling near centers of functions");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("EXAMPLES:");
  optParser->AddUsageLine("   multimfit -c model_config_n100a.dat -i imageinfo_n100a.dat");
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
  optParser->AddFlag("logging");
  optParser->AddOption("save-params-root");
  optParser->AddOption("save-model");
  optParser->AddOption("save-residual");
  optParser->AddOption("save-weights");
  optParser->AddOption("ftol");
  optParser->AddOption("bootstrap");
  optParser->AddOption("save-bootstrap");
  optParser->AddOption("config", "c");
  optParser->AddOption("image-info", "i");
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
  // For multimfit, we currently don't have any arguments


  /* Process the results: options */
  // First four are options which print useful info and then exit the program
  if ( optParser->FlagSet("help") || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("version") ) {
    printf("multimfit version %s\n\n", VERSION_STRING);
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
  if (optParser->FlagSet("logging")) {
    theOptions->loggingOn = true;
    printf("\tLogging to file %s ...\n", LOG_FILENAME);
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
  if (optParser->OptionSet("image-info")) {
    theOptions->imageInfoFile = optParser->GetTargetString("image-info");
    theOptions->imageInfoFile_exists = true;
  }
  if (optParser->OptionSet("save-model")) {
    theOptions->outputModelFileName = optParser->GetTargetString("save-model");
    theOptions->saveModel = true;
    printf("\toutput best-fit model image root name = %s\n", theOptions->outputModelFileName.c_str());
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
  if (optParser->OptionSet("save-params-root")) {
    theOptions->outputParameterRootName = optParser->GetTargetString("save-params-root");
    theOptions->saveBestFitParams = true;
    printf("\troot name for output best-fit parameter files = %s\n", theOptions->outputParameterRootName.c_str());
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
//     string  dataImageFileName;
//     string  maskImageFileName;
//     string  errorImageFileName;
//     string  psfImageFileName;
//     string  psfOversampledFileName;
//     double  weight, gain, readNoise, originalSky;
//     double  x0, y0;
//     int  nColumns, nRows;
//     double  pixelScale, intensityScale, rotation;
// 
//     bool  dataImage_present, maskImage_present, errorImage_present;
//     bool  psfImage_present, psfOversampledImage_present;
bool RequestedFilesPresent( vector<ImageInfo> &imageInfoVect, int nImages )
{
  bool  allFilesPresent = true;
  
  for (int i = 0; i < nImages; i++) {
    if (imageInfoVect[i].dataImage_present)
      if (! ImageFileExists(imageInfoVect[i].dataImageFileName.c_str())) {
        fprintf(stderr, "\n*** ERROR: Unable to find data image file \"%s\"!\n", 
                imageInfoVect[i].dataImageFileName.c_str());
        allFilesPresent = false;
      }
    if (imageInfoVect[i].maskImage_present)
      if (! ImageFileExists(imageInfoVect[i].maskImageFileName.c_str())) {
        fprintf(stderr, "\n*** ERROR: Unable to find mask image file \"%s\"!\n", 
                imageInfoVect[i].maskImageFileName.c_str());
        allFilesPresent = false;
      }
    if (imageInfoVect[i].errorImage_present)
      if (! ImageFileExists(imageInfoVect[i].errorImageFileName.c_str())) {
        fprintf(stderr, "\n*** ERROR: Unable to find error image file \"%s\"!\n", 
                imageInfoVect[i].errorImageFileName.c_str());
        allFilesPresent = false;
      }
    if (imageInfoVect[i].psfImage_present)
      if (! ImageFileExists(imageInfoVect[i].psfImageFileName.c_str())) {
        fprintf(stderr, "\n*** ERROR: Unable to find PSF image file \"%s\"!\n", 
                imageInfoVect[i].psfImageFileName.c_str());
        allFilesPresent = false;
      }
      
    if (imageInfoVect[i].psfOversampledImage_present)
      for (int n = 0; n < imageInfoVect[i].nOversampledFileNames; n++)
        if (! ImageFileExists(imageInfoVect[i].psfOversampledFileNames[n].c_str())) {
          fprintf(stderr, "\n*** ERROR: Unable to find oversampled PSF image file \"%s\"!\n", 
                imageInfoVect[i].psfOversampledFileNames[n].c_str());
          allFilesPresent = false;
        }
  }

  return allFilesPresent;
}



#ifndef NO_SIGNALS
void signal_handler( int signal )
{
  if (signal == SIGINT)
    stopSignal_flag = 1;
}
#endif


/* END OF FILE: multimfit_main.cpp --------------------------------------- */
