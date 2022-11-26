/* FILE: makeimage_main.cpp ---------------------------------------------- */
/*
 * Program for generating model images, using same code and approaches as
 * imfit does (but without the image-fitting parts).
 * 
 * The proper translations are:
 * NAXIS1 = naxes[0] = nColumns = sizeX;
 * NAXIS2 = naxes[1] = nRows = sizeY.
*/

// Copyright 2010--2022 by Peter Erwin.
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
#include <string.h>
#include <math.h>
#include <string>
#include <memory>
#include <sys/time.h>
#include "fftw3.h"

// logging
#ifdef USE_LOGGING
#include "loguru/loguru.hpp"
#endif

#include "definitions.h"
#include "image_io.h"
#include "getimages.h"
#include "model_object.h"
#include "add_functions.h"
#include "options_base.h"
#include "options_makeimage.h"
#include "commandline_parser.h"
#include "config_file_parser.h"
#include "utilities_pub.h"
#include "sample_configs.h"
#include "psf_oversampling_info.h"
#include "setup_model_object.h"
#include "count_cpu_cores.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
#define EST_SIZE_HELP_STRING "     --estimation-size <int>  Size of square image to use for estimating fluxes [default = 5000]"

const string  LOG_FILENAME = "log_makeimage.txt";


// Option names for use in config files
static string  kNCols1 = "NCOLS";
static string  kNCols2 = "NCOLUMNS";
static string  kNRows = "NROWS";


#ifdef USE_OPENMP
#define VERSION_STRING      "1.9.0 (OpenMP-enabled)"
#else
#define VERSION_STRING      "1.9.0"
#endif



/* ------------------- Function Prototypes ----------------------------- */

void ProcessInput( int argc, char *argv[], shared_ptr<MakeimageOptions> theOptions );
void HandleConfigFileOptions( configOptions *configFileOptions, 
								shared_ptr<MakeimageOptions> mainOptions );
void DetermineFluxes( ModelObject *theModel, double *parameters, 
					shared_ptr<MakeimageOptions> options, const string &programName );





/* ---------------- MAIN ----------------------------------------------- */

int main( int argc, char *argv[] )
{
  int  nColumns, nRows;
  int  nColumns_psf, nRows_psf;
  int  nParamsTot;
  int  status;
  double  *psfPixels = NULL;
  long  nPixels_psf_oversampled;
  vector<PsfOversamplingInfo *>  psfOversamplingInfoVect;
  double  *paramsVect;
  ModelObject  *theModel;
  vector<string>  functionList;
  vector<string>  functionLabelList;
  vector<double>  parameterList;
  vector<int>  functionSetIndices;
  vector< map<string, string> > optionalParams;
  vector<string>  imageCommentsList;
  double  *singleFunctionImage;
  shared_ptr<MakeimageOptions> options;
  configOptions  userConfigOptions;
  bool  printFluxesOnly = false;
  
  string  progName = "makeimage ";
  progName += VERSION_STRING;

  
  
  options = make_shared<MakeimageOptions>();    
  // Set maximum number of threads = number of hardware cores by default
  // (user can still override this with --max-threads option)
  options->maxThreads = GetPhysicalCoreCount();
  options->maxThreadsSet = true;
  /* Process command line and parse config file: */
  ProcessInput(argc, argv, options);

#ifdef USE_LOGGING
  if (options->loggingOn) {
    // turn off writing log output to stderr and remove "thread name" from outputs
    loguru::g_stderr_verbosity = loguru::Verbosity_OFF;
    loguru::g_preamble_thread = false;
    // initialize logging
    loguru::init(argc, argv);
    loguru::add_file(LOG_FILENAME.c_str(), loguru::Append, loguru::Verbosity_MAX);
    LOG_F(INFO, "*** Starting up...");
  }
#endif
  
  if ( ((options->printFluxes) || options->saveFluxes) && 
  		(! options->saveImage) && (! options->printImages) )
    printFluxesOnly = true;


  // Read configuration file
  if (! FileExists(options->configFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find configuration file \"%s\"!\n\n", 
           options->configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options->configFileName, true, functionList, 
  							functionLabelList, parameterList, functionSetIndices, 
  							userConfigOptions, optionalParams);
  if (status != 0) {
    fprintf(stderr, "\n*** ERROR: Failure reading configuration file \"%s\"!\n\n", 
    			options->configFileName.c_str());
    return -1;
  }

  // Parse and process user-supplied (non-function) values from config file, if any
  HandleConfigFileOptions(&userConfigOptions, options);


  if (! options->saveImage) {
    printf("\nUser requested that no images be saved!\n\n");
  }
  
  
  // Figure out size of model image
  // First, figure out if we have enough information, given what user wants us to do
  // (e.g., if user wants image saved, we need nColumns and nRows, or else a reference image)
  if ((options->nColumns > 0) && (options->nRows > 0))
    options->noImageDimensions = false;
  if ( (options->noRefImage) && (options->noImageDimensions)) {
    if (options->saveImage) {
      fprintf(stderr, "\n*** ERROR: Insufficient image dimensions (or no reference image) supplied!\n");
      fprintf(stderr, "           (Use --nrows and --ncols, or else --refimage)\n\n");
      return -1;
    }
    else {
      // minimal image size for the case where we don't actually generate an image
      // (except for --print-fluxes purposes)
      options->nColumns = 2;
      options->nRows = 2;
    }
  }
  // Get image size from reference image, if necessary
  if ((! printFluxesOnly) && (options->noImageDimensions)) {
    //status = GetImageSize(options->referenceImageName, &nColumns, &nRows);
    std::tie(nColumns, nRows, status) = GetImageSize(options->referenceImageName);
    if (status != 0) {
      fprintf(stderr,  "\n*** ERROR: Failure determining size of image file \"%s\"!\n\n", 
      			options->referenceImageName.c_str());
      exit(-1);
    }
    // Reminder: nColumns = n_pixels_per_row
    // Reminder: nRows = n_pixels_per_column
    printf("Reference image read: naxis1 [# rows] = %d, naxis2 [# columns] = %d\n",
           nRows, nColumns);
  }
  else {
    nColumns = options->nColumns;
    nRows = options->nRows;
  }
  
  // Read in PSF image, if supplied
  if (options->psfImagePresent) {
    std::tie(psfPixels, nColumns_psf, nRows_psf, status) = GetPsfImage(options->psfFileName);
    if (status < 0)
      exit(-1);
  }
  else
    printf("* No PSF image supplied -- no image convolution will be done!\n");

  // Read in oversampled PSF image(s), if supplied
  if ((options->psfOversampling) && (options->psfOversampledImagePresent)) {
    status = GetOversampledPsfInfo(options, 0,0, psfOversamplingInfoVect);
	if (status < 0)
	  exit(-1);
  }

  if (! options->subsamplingFlag)
    printf("* Pixel subsampling has been turned OFF.\n");


  
  // Set up the model object
  // Populate the column-and-row-numbers vector
  vector<int> nColumnsRowsVect;
  nColumnsRowsVect.push_back(nColumns);
  nColumnsRowsVect.push_back(nRows);
  nColumnsRowsVect.push_back(nColumns_psf);
  nColumnsRowsVect.push_back(nRows_psf);

  theModel = SetupModelObject(options, nColumnsRowsVect, NULL, psfPixels, NULL, NULL,
  								psfOversamplingInfoVect);

  // Add functions to the model object; also tells model object where function sets start
  status = AddFunctions(theModel, functionList, functionLabelList, functionSetIndices, 
  						options->subsamplingFlag, 0, optionalParams);
  if (status < 0) {
  	fprintf(stderr, "*** ERROR: Failure in AddFunctions!\n\n");
  	exit(-1);
  }

  theModel->PrintDescription();


  // Set up parameter vector(s), now that we know how many total parameters there are
  nParamsTot = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  if (nParamsTot != (int)parameterList.size()) {
  	fprintf(stderr, "*** ERROR: number of input parameters (%d) does not equal", 
  	       (int)parameterList.size());
  	fprintf(stderr, " required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	exit(-1);
 }
    
  // Copy parameters into C array and generate the model image
  paramsVect = (double *) calloc(nParamsTot, sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];

  if (! printFluxesOnly) {

#ifdef USE_LOGGING
    if (options->loggingOn)
      LOG_F(INFO, "Creating model image...");
#endif
    // OK, we're generating a normal model image
    theModel->CreateModelImage(paramsVect);
  
    // TESTING (remove later)
    if (options->printImages)
      theModel->PrintModelImage();

    /* Save model image: */
    if (options->saveImage) {
      PrepareImageComments(&imageCommentsList, progName, options->configFileName,
    					options->psfImagePresent, options->psfFileName);
      printf("\nSaving output model image (\"%s\") ...\n", options->outputImageName.c_str());
      status = SaveVectorAsImage(theModel->GetModelImageVector(), options->outputImageName, 
                        nColumns, nRows, imageCommentsList);
      if (status != 0) {
        fprintf(stderr,  "\n*** WARNING: Unable to save output image file \"%s\"!\n\n", 
        			options->outputImageName.c_str());
      }
      // code for checking PSF convolution fixes [May 2012]
      if (options->saveExpandedImage) {
        string  tempName = "expanded_" + options->outputImageName;
        printf("\nSaving full (expanded) output model image (\"%s\") ...\n", tempName.c_str());
        status = SaveVectorAsImage(theModel->GetExpandedModelImageVector(), tempName, 
                          nColumns + 2*nColumns_psf, nRows + 2*nRows_psf, imageCommentsList);
        if (status != 0) {
          fprintf(stderr,  "\n*** WARNING: Unable to save output image file \"%s\"!\n\n", 
          			tempName.c_str());
        }
      }
    }
  
    // Save individual-function images, if requested
    if ((options->saveImage) && (options->saveAllFunctions)) {
      vector<string> functionNames;
      vector<string> functionLabels;
      int  nFuncs = theModel->GetNFunctions();
      theModel->GetFunctionNames(functionNames);
      theModel->GetFunctionLabels(functionLabels);
      for (int i = 0; i < nFuncs; i++) {
        // Generate single-function image (exit if that failed -- e.g., due to
        // memory allocation failure)
        string currentFilename, headerString;
        singleFunctionImage = theModel->GetSingleFunctionImage(paramsVect, i);
        if (singleFunctionImage == NULL) {
          fprintf(stderr, "\n*** ERROR: Unable to generate single-function image #%d!\n\n", i);
          exit(-1);
        }
        currentFilename = PrintToString("%s%d_%s.fits", options->functionRootName.c_str(),
        								i + 1, functionNames[i].c_str());
        printf("%s\n", currentFilename.c_str());
        // Add comments for FITS header, describing this function
        headerString = PrintToString("FUNCTION %s", functionNames[i].c_str());
        if (! functionLabels[i].empty())
          headerString = PrintToString("%s # LABEL %s", headerString.c_str(), 
          								functionLabels[i].c_str());
        imageCommentsList.push_back(headerString);
        status = SaveVectorAsImage(singleFunctionImage, currentFilename, nColumns, nRows, 
        							imageCommentsList);
        if (status != 0) {
          fprintf(stderr,  "\n*** WARNING: Unable to save output single-function image file \"%s\"!\n\n", 
          			currentFilename.c_str());
        }
        imageCommentsList.pop_back();
      }
    }
  }
  
  
  // Estimate component fluxes, if requested
  if ((options->printFluxes) || (options->saveFluxes)) {
    DetermineFluxes(theModel, paramsVect, options, progName);
    printf("\n");
  }


  // Estimate image computation time
  if (options->timingIterations > 0) {
    struct timeval  timer_start, timer_end;
    double  microsecs, time_elapsed, time_per_iteration;
    gettimeofday(&timer_start, NULL);
    for (int i = 0; i < options->timingIterations; i++)
      theModel->CreateModelImage(paramsVect);
    gettimeofday(&timer_end, NULL);
    microsecs = timer_end.tv_usec - timer_start.tv_usec;
    time_elapsed = timer_end.tv_sec - timer_start.tv_sec + microsecs/1e6;
    time_per_iteration = time_elapsed / options->timingIterations;
    printf("\nElapsed time: %.6f sec\n", time_elapsed);
    printf("\tMean time per image computation (from %d iterations) = %.7f sec\n", 
    		options->timingIterations, time_per_iteration);
  }

  
  printf("Done!\n\n");


#ifdef USE_LOGGING
  if (options->loggingOn) {
    LOG_F(INFO, "*** Freeing up memory...");
    printf("\nInternal logging output saved to %s\n\n", LOG_FILENAME.c_str());
  }
#endif

  // Free up memory
  if (options->psfImagePresent)
    fftw_free(psfPixels);       // allocated in ReadImageAsVector()
  if (psfOversamplingInfoVect.size() > 0) {
    for (int nn = 0; nn < (int)psfOversamplingInfoVect.size(); nn++)
      free(psfOversamplingInfoVect[nn]);
    psfOversamplingInfoVect.clear();
  }
  free(paramsVect);
  delete theModel;
  
  return 0;
}



void ProcessInput( int argc, char *argv[], shared_ptr<MakeimageOptions> theOptions )
{

  CLineParser *optParser = new CLineParser();
  string  tempString = "";

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   makeimage [options] config-file");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine(" -v  --version                Prints version number");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --list-parameters        Prints list of parameter names for each available function");
  tempString = PrintToString("     --sample-config          Generates an example configuration file (%s)", configMakeimageFile.c_str());
  optParser->AddUsageLine(tempString);
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" -o  --output <output-image.fits>        name for output image [default = modelimage.fits]");
  optParser->AddUsageLine("     --refimage <reference-image.fits>   reference image (for image size)");
  optParser->AddUsageLine("     --psf <psf.fits>                    PSF image to use (for convolution)");
  optParser->AddUsageLine("     --no-normalize                      Do *not* normalize input PSF image");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     (Note that the following 3 options can be specified multiple times)");
  optParser->AddUsageLine("     --overpsf <psf.fits>                Oversampled PSF image to use");
  optParser->AddUsageLine("     --overpsf_scale <n>                 Oversampling scale (integer)");
  optParser->AddUsageLine("     --overpsf_region <x1:x2,y1:y2>      Section of image to convolve with oversampled PSF");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --ncols <number-of-columns>         x-size of output image");
  optParser->AddUsageLine("     --nrows <number-of-rows>            y-size of output image");
  optParser->AddUsageLine("     --no-subsampling                    Do *not* do pixel subsampling near centers");
//  optParser->AddUsageLine("     --printimage             Print out images (for debugging)");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --output-functions <root-name>      Output individual-function images");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --print-fluxes           Estimate total component fluxes (& magnitudes, if zero point is given)");
  optParser->AddUsageLine(EST_SIZE_HELP_STRING);
  optParser->AddUsageLine("     --zero-point <value>     Zero point (for estimating component & total magnitudes)");
  optParser->AddUsageLine("     --save-fluxes <filename>            Save print-fluxes output to user-specified file");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --nosave                 Do *not* save image (for testing)");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --timing <int>           Generate image specified number of times and estimate average creation time");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --max-threads <int>      Maximum number of threads to use");
  optParser->AddUsageLine("");
#ifdef USE_LOGGING
  optParser->AddUsageLine("     --logging                Save logging outputs to file");
  optParser->AddUsageLine("");
#endif
  optParser->AddUsageLine("     --debug <n>              Set the debugging level (integer)");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("EXAMPLES:");
  optParser->AddUsageLine("   makeimage model_config_a.dat");
  optParser->AddUsageLine("   makeimage model_config_b.dat --ncols 800 --nrows 800 --psf best_psf.fits -o testimage_convolved.fits");
  optParser->AddUsageLine("   makeimage bestfit_parameters.dat --print-fluxes --zero-point 26.24");
  optParser->AddUsageLine("");

  optParser->AddFlag("help", "h");
  optParser->AddFlag("version", "v");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("list-parameters");
  optParser->AddFlag("sample-config");
  optParser->AddFlag("printimage");
  optParser->AddFlag("save-expanded");
  optParser->AddFlag("no-normalize");
  optParser->AddFlag("no-subsampling");
  optParser->AddFlag("print-fluxes");
  optParser->AddFlag("nosave");
  optParser->AddOption("output", "o");
  optParser->AddOption("ncols");
  optParser->AddOption("nrows");
  optParser->AddOption("refimage");
  optParser->AddOption("psf");
  optParser->AddQueueOption("overpsf");
  optParser->AddQueueOption("overpsf_scale");
  optParser->AddQueueOption("overpsf_region");
  optParser->AddOption("zero-point");
  optParser->AddOption("estimation-size");
  optParser->AddOption("save-fluxes");
  optParser->AddOption("output-functions");
  optParser->AddOption("timing");
  optParser->AddOption("max-threads");
  optParser->AddOption("debug");
#ifdef USE_LOGGING
  optParser->AddFlag("logging");
#endif
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
    theOptions->configFileName = optParser->GetArgument(0);
    theOptions->noConfigFile = false;
  }

  /* Process the results: flags and options */
  // First two are options which print useful info and then exit the program
  if ( optParser->FlagSet("help") || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("version") ) {
    printf("makeimage version %s\n\n", VERSION_STRING);
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
    int saveStatus = SaveExampleMakeimageConfig();
    if (saveStatus == 0)
      printf("Sample configuration file \"%s\" saved.\n", configMakeimageFile.c_str());
    delete optParser;
    exit(1);
  }

  if (optParser->FlagSet("printimage")) {
    theOptions->printImages = true;
  }
  if (optParser->FlagSet("save-expanded")) {
    theOptions->saveExpandedImage = true;
  }
  if (optParser->FlagSet("no-normalize")) {
    theOptions->normalizePSF = false;
  }
  if (optParser->FlagSet("no-subsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (optParser->FlagSet("nosave")) {
    theOptions->saveImage = false;
  }
  if (optParser->FlagSet("print-fluxes")) {
    theOptions->printFluxes = true;
    theOptions->saveImage = false;
  }
  if (optParser->OptionSet("output")) {
    theOptions->outputImageName = optParser->GetTargetString("output");
    theOptions->noOutputImageName = false;
  }
  if (optParser->OptionSet("refimage")) {
    theOptions->referenceImageName = optParser->GetTargetString("refimage");
    theOptions->noRefImage = false;
  }
  if (optParser->OptionSet("psf")) {
    theOptions->psfFileName = optParser->GetTargetString("psf");
    theOptions->psfImagePresent = true;
    printf("\tPSF image = %s\n", theOptions->psfFileName.c_str());
  }
  
  // oversampled PSF(s) and region(s)
  if (optParser->OptionSet("overpsf")) {
    theOptions->psfOversampling = true;
    theOptions->psfOversampledImagePresent = true;
    for (int i = 0; i < optParser->GetNTargets("overpsf"); i++) {
      string fileName = optParser->GetTargetString("overpsf", i);
      theOptions->psfOversampledFileNames.push_back(fileName);
      printf("\tOversampled PSF image = %s\n", fileName.c_str());
    }
  }
  if (optParser->OptionSet("overpsf_scale")) {
    for (int i = 0; i < optParser->GetNTargets("overpsf_scale"); i++) {
      string scaleStr = optParser->GetTargetString("overpsf_scale", i).c_str();
      if (NotANumber(scaleStr.c_str(), 0, kPosInt)) {
        fprintf(stderr, "*** ERROR: overpsf_scale should be a positive integer!\n");
        delete optParser;
        exit(1);
      }
      int scale = atoi(scaleStr.c_str());
      theOptions->psfOversamplingScales.push_back(scale);
      printf("\tPSF oversampling scale = %d\n", scale);
    }
  }
  if (optParser->OptionSet("overpsf_region")) {
    theOptions->oversampleRegionSet = true;
    for (int i = 0; i < optParser->GetNTargets("overpsf_region"); i++) {
      string psfRegion = optParser->GetTargetString("overpsf_region", i);
      theOptions->psfOversampleRegions.push_back(psfRegion);
      theOptions->nOversampleRegions += 1;
      printf("\tPSF oversampling region = %s\n", psfRegion.c_str());
    }
  }
  
  if (optParser->OptionSet("ncols")) {
    if (NotANumber(optParser->GetTargetString("ncols").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: ncols should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->nColumns = atol(optParser->GetTargetString("ncols").c_str());
    theOptions->nColumnsSet = true;
  }
  if (optParser->OptionSet("nrows")) {
    if (NotANumber(optParser->GetTargetString("nrows").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: nrows should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->nRows = atol(optParser->GetTargetString("nrows").c_str());
    theOptions->nRowsSet = true;
  }
  if (optParser->OptionSet("zero-point")) {
    if (NotANumber(optParser->GetTargetString("zero-point").c_str(), 0, kAnyReal)) {
      fprintf(stderr, "*** ERROR: zero point should be a real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->magZeroPoint = atof(optParser->GetTargetString("zero-point").c_str());
    printf("\tmagnitude zero point = %g\n", theOptions->magZeroPoint);
  }
  if (optParser->OptionSet("estimation-size")) {
    if (NotANumber(optParser->GetTargetString("estimation-size").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: estimation size should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->estimationImageSize = atol(optParser->GetTargetString("estimation-size").c_str());
  }
  if (optParser->OptionSet("save-fluxes")) {
    theOptions->saveFluxesFileName = optParser->GetTargetString("save-fluxes");
    theOptions->saveFluxes = true;
  }
  if (optParser->OptionSet("output-functions")) {
    theOptions->functionRootName = optParser->GetTargetString("output-functions");
    theOptions->saveAllFunctions = true;
  }
  if (optParser->OptionSet("timing")) {
    if (NotANumber(optParser->GetTargetString("timing").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: timing should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->timingIterations = atol(optParser->GetTargetString("timing").c_str());
    theOptions->saveImage = false;
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
  if (optParser->OptionSet("debug")) {
    if (NotANumber(optParser->GetTargetString("debug").c_str(), 0, kAnyInt)) {
      fprintf(stderr, "*** ERROR: debug should be an integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->debugLevel = atol(optParser->GetTargetString("debug").c_str());
  }
#ifdef USE_LOGGING
  if (optParser->FlagSet("logging")) {
    theOptions->loggingOn = true;
  }
#endif

  if ((theOptions->nColumns) && (theOptions->nRows))
    theOptions->noImageDimensions = false;
  
  delete optParser;

}


// Note that we only use options from the config file if they have *not* already been set
// by the command line (i.e., command-line options override config-file values).
//void HandleConfigFileOptions( configOptions *configFileOptions, makeimageCommandOptions *mainOptions )
void HandleConfigFileOptions( configOptions *configFileOptions, 
								shared_ptr<MakeimageOptions> mainOptions )
{
	int  newIntVal;
	
  if (configFileOptions->nOptions == 0)
    return;

  for (int i = 0; i < configFileOptions->nOptions; i++) {
    
    if ((configFileOptions->optionNames[i] == kNCols1) || 
    		(configFileOptions->optionNames[i] == kNCols2)) {
      if (mainOptions->nColumnsSet) {
        printf("nColumns (x-size of image) value in config file ignored (using command-line value)\n");
      } else {
        if (NotANumber(configFileOptions->optionValues[i].c_str(), 0, kPosInt)) {
          fprintf(stderr, "*** ERROR: NCOLS should be a positive integer!\n");
          exit(1);
        }
        newIntVal = atoi(configFileOptions->optionValues[i].c_str());
        printf("Value from config file: nColumns = %d\n", newIntVal);
        mainOptions->nColumns = newIntVal;
      }
      continue;
    }
    if (configFileOptions->optionNames[i] == kNRows) {
      if (mainOptions->nRowsSet) {
        printf("nRows (y-size of image) value in config file ignored (using command-line value)\n");
      } else {
        if (NotANumber(configFileOptions->optionValues[i].c_str(), 0, kPosInt)) {
          fprintf(stderr, "*** ERROR: NROWS should be a positive integer!\n");
          exit(1);
        }
        newIntVal = atoi(configFileOptions->optionValues[i].c_str());
        printf("Value from config file: nRows = %d\n", newIntVal);
        mainOptions->nRows = newIntVal;
      }
      continue;
    }
    // we only get here if we encounter an unknown option
    printf("Unknown keyword (\"%s\") in config file ignored\n", 
    				configFileOptions->optionNames[i].c_str());
    
  }
}


void DetermineFluxes( ModelObject *theModel, double *parameters, 
					shared_ptr<MakeimageOptions> options, const string &programName )
{
  int  n;
  int  nComponents = theModel->GetNFunctions();
  double *fluxes = (double *) calloc(nComponents, sizeof(double));
  double  fraction, totalFlux, magnitude;
  vector<string> functionNames;
  vector<string> functionLabels;
  vector<string> saveFluxesLines;
  
  const string columHeaderLine = "Component                 Flux        Magnitude   Fraction   Label\n";
  const string columHeaderLine_for_file = "# Component               Flux        Magnitude   Fraction   Label\n";
  
  theModel->GetFunctionNames(functionNames);
  theModel->GetFunctionLabels(functionLabels);
  
  if (options->printFluxes) {
    printf("\nEstimating fluxes on %d x %d image", options->estimationImageSize,
  													options->estimationImageSize);
    if (options->magZeroPoint != NO_MAGNITUDES)
      printf(" (using zero point = %g)", options->magZeroPoint);
    printf("...\n\n");
  }
  
  totalFlux = theModel->FindTotalFluxes(parameters, options->estimationImageSize,
  										options->estimationImageSize, fluxes);
  for (n = 0; n < (nComponents + 1); n++) {
    double  thisFlux;
    string  thisLine = "";
    string  componentName;
    if (n == nComponents) {
      thisFlux = totalFlux;
      fraction = 1.0;
      componentName = "Total";
      thisLine = "\n";   // enforce line-skip before "Total" results
    }
    else {
      thisFlux = fluxes[n];
      fraction = thisFlux / totalFlux;
      componentName = functionNames[n];
    }
    thisLine += PrintToString("%-25s %10.4e", componentName.c_str(), thisFlux);
    if ( (options->magZeroPoint != NO_MAGNITUDES) && (thisFlux > 0.0) ) {
      magnitude = options->magZeroPoint - 2.5*log10(thisFlux);
      thisLine += PrintToString("   %6.4f", magnitude);
    }
    else
      thisLine += PrintToString("      ---");
    if ((fraction >= 0.0001) || (fraction == 0.0))
      thisLine += PrintToString("%12.5f", fraction);
    else
      thisLine += PrintToString("%12.3e", fraction);
    if (n < nComponents) {
      if (! functionLabels[n].empty())
        thisLine += PrintToString("    %s", functionLabels[n].c_str());
    }
    thisLine += "\n";
    saveFluxesLines.push_back(thisLine);
  }

  if (options->printFluxes) {
    // print results to console
    printf("%s", columHeaderLine.c_str());
    for (n = 0; n < (int)saveFluxesLines.size(); n++)
      printf("%s", saveFluxesLines[n].c_str());
    printf("\n");
  }
  if (options->saveFluxes) {
    // save results to file, if requested
    char  *timeStamp = TimeStamp();
    vector<string> headerLines;
    headerLines.push_back(PrintToString("# Fluxes estimated on %s by %s\n", timeStamp, 
    									programName.c_str()));
    headerLines.push_back(PrintToString("# Using config file %s\n", 
    									options->configFileName.c_str()));
    if (options->magZeroPoint != NO_MAGNITUDES)
    headerLines.push_back(PrintToString("# (photometric zero point = %g)\n", options->magZeroPoint));
    headerLines.push_back("#\n");
    headerLines.push_back(columHeaderLine_for_file);

    FILE * file_ptr = fopen(options->saveFluxesFileName.c_str(), "w");
    for (n = 0; n < (int)headerLines.size(); n++)
      fprintf(file_ptr, "%s", headerLines[n].c_str());
    for (n = 0; n < (int)saveFluxesLines.size(); n++)
      fprintf(file_ptr, "%s", saveFluxesLines[n].c_str());
    fclose(file_ptr);
    printf("\nEstimated fluxes saved to file %s.\n", options->saveFluxesFileName.c_str());
  }

  free(fluxes);
}



/* END OF FILE: makeimage_main.cpp --------------------------------------- */
