/* FILE: makemultiimages_main.cpp ---------------------------------------- */
/*
 * Program for generating model images, using same code and approaches as
 * imfit does (but without the image-fitting parts).
 * 
 * MODIFIED FOR MULTIMFIT -- test code for making images in multimfit mode
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

#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <sys/time.h>
#include <tuple>
#include "fftw3.h"

// logging
#ifdef USE_LOGGING
#include "loguru/loguru.hpp"
#endif

#include "definitions.h"
#include "image_io.h"
#include "getimages.h"
#include "model_object.h"
#include "model_object_multimage.h"
#include "add_functions.h"
#include "options_base.h"
#include "options_makemultimages.h"
#include "commandline_parser.h"
#include "config_file_parser.h"
#include "utilities_pub.h"
#include "sample_configs.h"
#include "psf_oversampling_info.h"
#include "store_psf_oversampling.h"
#include "setup_model_object.h"
#include "count_cpu_cores.h"

#include "imageparams_file_parser.h"


/* ---------------- Definitions ---------------------------------------- */
#define EST_SIZE_HELP_STRING "     --estimation-size <int>  Size of square image to use for estimating fluxes [default = 5000]"

const string  LOG_FILENAME = "log_makemultimages.txt";


// Option names for use in config files
static string  kNCols1 = "NCOLS";
static string  kNCols2 = "NCOLUMNS";
static string  kNRows = "NROWS";


#ifdef USE_OPENMP
#define VERSION_STRING      "0.8 (OpenMP-enabled)"
#else
#define VERSION_STRING      "0.8"
#endif





/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void ProcessInput( int argc, char *argv[], shared_ptr<MakeMultimagesOptions> theOptions );



/* ------------------------ Module Variables --------------------------- */





/* ---------------- MAIN ----------------------------------------------- */

int main( int argc, char *argv[] )
{
  int  nColumns_thisImage, nRows_thisImage;
  int  nColumns_psf, nRows_psf;
  int  nParamsTot;
  int  status;
  double  *psfPixels = nullptr;
  vector<double *> psfPixelsUsed;  // so we know how many allocated pixel-vectors to free
  vector<PsfOversamplingInfo *>  psfOversamplingInfoVect;
  vector< vector<PsfOversamplingInfo *> > allPsfOversamplingInfoVectors;
  double  *paramsVect;
  ModelObjectMultImage  *theMultImageModel;
  double *imageDescParamsVect;
  int  nImageDescParams;
  int  nGlobalParams = 0;
  vector<double>  localParamsList;
  int  nLocalParams = 0;
  vector<ImageInfo> imageInfoVect;
  int nImages = 0;
  int  nGlobalFuncs = 0;
  vector<string>  functionList;
  vector<string>  functionLabelList;
  vector<double>  parameterList;
  vector<int>  functionSetIndices;
  vector< map<string, string> > optionalParamsMap;
  vector<string>  imageCommentsList;
  shared_ptr<MakeMultimagesOptions> options;
  configOptions  userConfigOptions;
  
  /* Process command line and parse config file: */
  options = make_shared<MakeMultimagesOptions>();    
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

  
  if (! FileExists(options->configFileName.c_str())) {
    fprintf(stderr, "\n*** ERROR: Unable to find configuration file \"%s\"!\n\n", 
           options->configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options->configFileName, true, functionList, 
  							functionLabelList, parameterList, functionSetIndices, 
  							userConfigOptions);
  if (status != 0) {
    fprintf(stderr, "\n*** ERROR: Failure reading configuration file \"%s\"!\n\n", 
    			options->configFileName.c_str());
    return -1;
  }
  nGlobalFuncs = functionList.size();

  if (options->imageInfoFile_exists) {
    if (! FileExists(options->imageInfoFile.c_str())) {
      fprintf(stderr, "\n*** ERROR: Unable to find image-info file \"%s\"!\n\n", 
             options->imageInfoFile.c_str());
      return -1;
    }
    status = ReadImageParamsFile(options->imageInfoFile, nImages, imageInfoVect);
    std::tie(imageDescParamsVect, nImageDescParams) = GetImageDescriptionParams(imageInfoVect);
    printf("Image-description params from image-info file \"%s\": ", options->imageInfoFile.c_str());
    for (int i = 0; i < nImageDescParams; i++)
      printf("%.1f ", imageDescParamsVect[i]);
    printf("\n");
  }
  else {
      fprintf(stderr, "\n*** ERROR: No image-parameters file was specified! (Use --image-info)\n\n");
      return -1;
  }

  if (! options->saveImage) {
    printf("\nUser requested that no images be saved!\n\n");
  }
  

  // Creations and initial setup of ModelObjectMultImage
  theMultImageModel = new ModelObjectMultImage();
  theMultImageModel->SetupGeneralOptions(options);
  
  
  // MULTIMFIT: loop over creation of ModelObject instances, one per image
  //    Create new ModelObject instance via SetupModelObject
  //    Call AddFunctions using new ModelObject instance
  //    Add new ModelObject instance to ModelObjectMultImage
  // After loop finishes, call AddFunctions using ModelObjectMultImage
  
  vector<int> nColumnsRowsVect;
  ModelObject * newModel;
  for (int i = 0; i < nImages; i++) {
    // reset things to default values
    options->psfImagePresent = false;
    options->psfOversampling = false;
    options->psfOversampledImagePresent = false;
    options->oversampleRegionSet = false;
    nColumnsRowsVect.clear();
    psfOversamplingInfoVect.clear();
    // reset these to point to nullptr, otherwise we'll get e.g. psf pixels from
    // previous rounds passed in to SetupModelObject
    psfPixels = nullptr;

    if ((imageInfoVect[i].nColumns > 0) && (imageInfoVect[i].nRows > 0)) {
      nColumns_thisImage = imageInfoVect[i].nColumns;
      nRows_thisImage = imageInfoVect[i].nRows;
    }
    else if (imageInfoVect[i].dataImageFileName != "") {
      // get image dimensions from data image
      std::tie(nColumns_thisImage, nRows_thisImage, status) = GetImageSize(imageInfoVect[i].dataImageFileName);
      if (status != 0) {
        fprintf(stderr,  "\n*** ERROR: image #%d: Failure determining size of image file \"%s\"!\n\n", 
      				i + 1, imageInfoVect[i].dataImageFileName.c_str());
        exit(-1);
      }
    }
    else {
      fprintf(stderr, "\n*** ERROR: Unable to determine image size for image");
      fprintf(stderr, " number %d in image-info file\n", i + 1);
      exit(-1);
    }

    // Read in PSF image, if supplied
    if (imageInfoVect[i].psfImageFileName != "") {
      printf("Reading PSF image (\"%s\") ...\n", imageInfoVect[i].psfImageFileName.c_str());
      std::tie(psfPixels, nColumns_psf, nRows_psf, status) = GetPsfImage(imageInfoVect[i].psfImageFileName);
      if (status < 0)
        exit(-1);
      options->psfImagePresent = true;
      psfPixelsUsed.push_back(psfPixels);
    }

    // Read in oversampled PSF image(s), if supplied
    if (imageInfoVect[i].psfOversampledImage_present) {
      options->psfOversampling = true;
      status = ExtractAndStorePsfOversampling(imageInfoVect[i], i, psfOversamplingInfoVect);
      if (status < 0)
        exit(status);
    }
    // store copy of current psfOversamplingInfoVect so we can properly de-allocate 
    // its memory later
    allPsfOversamplingInfoVectors.push_back(psfOversamplingInfoVect);
    
    nColumnsRowsVect.push_back(nColumns_thisImage);
    nColumnsRowsVect.push_back(nRows_thisImage);
    nColumnsRowsVect.push_back(nColumns_psf);
    nColumnsRowsVect.push_back(nRows_psf);
    newModel = SetupModelObject(options, nColumnsRowsVect, nullptr, psfPixels, nullptr,
    						nullptr, psfOversamplingInfoVect);
    
    // Add functions to the model object; also tells model object where function sets start
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
    printf("function list: ");
    for (int k = 0; k < (int)currentFunctionList.size(); k++)
      printf(" %s", currentFunctionList[k].c_str());
    printf("\n");
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
    status = AddFunctions(newModel, currentFunctionList, currentFunctionLabelList, 
    					currentFunctionSetIndices, options->subsamplingFlag, 0, 
    					currentOptionalParamsMap, globalFuncFlags);
    if (status < 0) {
  	  fprintf(stderr, "*** ERROR: image #%d,  Failure in AddFunctions!\n\n", i + 1);
  	  exit(-1);
    }
    newModel->FinalModelSetup();
    // Note that ModelObjectMultImage will handle de-allocation of ModelObject instance
    theMultImageModel->AddModelObject(newModel);
  }

  printf("main: theMultImageModel has %d \"images\" (ModelObject instances)\n",
  		theMultImageModel->GetNImages());
  theMultImageModel->PrintDescription();
  
  // Finally, add (global) functions to ModelObjectMultImage object
  // NOTE: this does *not* add any of the per-image functions!
  status = AddFunctions(theMultImageModel, functionList, functionLabelList, functionSetIndices, 
  						options->subsamplingFlag, 0, optionalParamsMap);
  if (status < 0) {
  	fprintf(stderr, "*** ERROR: Failure in AddFunctions!\n\n");
  	exit(-1);
  }
  theMultImageModel->FinalModelSetup();
  


  // Set up parameter vector(s), now that we know how many total parameters there are
  nParamsTot = theMultImageModel->GetNParams();
  nGlobalParams = (int)parameterList.size();
  for (int i = 0; i < nImages; i++) {
    int  nLocalParams_thisImage = (int)imageInfoVect[i].perImageParamVals.size();
    nLocalParams += nLocalParams_thisImage;
    if (nLocalParams_thisImage > 0)
      localParamsList.insert(localParamsList.end(), imageInfoVect[i].perImageParamVals.begin(), 
      							imageInfoVect[i].perImageParamVals.end());
  }
  printf("\nModelObjectMultImage: nParamsTot = %d\n", nParamsTot);
  printf("nImageDescParams = %d\n", nImageDescParams);
  printf("nGlobalParams = %d\n", nGlobalParams);
  printf("nLocalParams = %d\n", nLocalParams);
  
  
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
    
  // Copy parameters into C array and generate the model image
  paramsVect = (double *) calloc(nInputParamsTot, sizeof(double));
  // stick image-description params at start of global parameter
  // vector, then add model params
  for (int i = 0; i < nImageDescParams; i++)
    paramsVect[i] = imageDescParamsVect[i];
  for (int i = nImageDescParams; i < nImageDescParams + nGlobalParams; i++)
    paramsVect[i] = parameterList[i - nImageDescParams];
  // finally, add per-image-function parameters, if they exist
  // nInputParams = nImageDescParams + nGlobalParams
  for (int i = nInputParams_nonLocal; i < nInputParamsTot; i++)
    paramsVect[i] = localParamsList[i - nInputParams_nonLocal];
  
  printf("\nparamsVect for CreateAllModelImages():\n");
  for (int i = 0; i < nInputParamsTot - 1; i++)
    printf(" %f,", paramsVect[i]);
  printf(" %f\n\n", paramsVect[nInputParamsTot - 1]);

#ifdef USE_LOGGING
    if (options->loggingOn)
      LOG_F(INFO, "Creating model images...");
#endif
  theMultImageModel->CreateAllModelImages(paramsVect);
  nImages = theMultImageModel->GetNImages();
  int  nColsThisImage, nRowsThisImage;
  for (int i = 0; i < nImages; i++) {
    std::tie(nColsThisImage, nRowsThisImage) = theMultImageModel->GetImageDimensions(i);
    string multImageOutputName = PrintToString("%s_%d.fits", options->outputImageName.c_str(), i);
    printf("MULTIMFIT: image %d: nCols,nRows = %d,%d\n", i, nColsThisImage, nRowsThisImage);
    printf("   Saving to: %s\n", multImageOutputName.c_str());
    status = SaveVectorAsImage(theMultImageModel->GetModelImageVector(i), 
    							multImageOutputName, nColsThisImage, nRowsThisImage,
    							imageCommentsList);
  }

  
  printf("Done!\n\n");

#ifdef USE_LOGGING
  if (options->loggingOn) {
    LOG_F(INFO, "*** Freeing up memory...");
    printf("\nInternal logging output saved to %s\n\n", LOG_FILENAME.c_str());
  }
#endif

  // Free up memory
  if (psfPixelsUsed.size() > 0) {
    for (int i = 0; i < (int)psfPixelsUsed.size(); i++)
      fftw_free(psfPixelsUsed[i]);  // allocated in ReadImageAsVector()
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
  free(paramsVect);
  free(imageDescParamsVect);
  delete theMultImageModel;
//   delete options;
  
  return 0;
}



void ProcessInput( int argc, char *argv[], shared_ptr<MakeMultimagesOptions> theOptions )
{

  CLineParser *optParser = new CLineParser();
  string  tempString = "";

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   makemultimages [options] config-file");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine(" -v  --version                Prints version number");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --list-parameters        Prints list of parameter names for each available function");
  tempString = PrintToString("     --sample-config          Generates an example configuration file (%s)", configMakeimageFile.c_str());
  optParser->AddUsageLine(tempString);
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --image-info <image-info-file> TEST-CODE: read in image info");
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" -o  --output <output-image-root>        root name for output image [default = modelimage_multi]");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --nosubsampling                     Do *not* do pixel subsampling near centers");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --output-functions <root-name>      Output individual-function images");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --nosave                 Do *not* save image (for testing, or for use with --print-fluxes)");
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
  optParser->AddUsageLine("   makemultimages model_config_a.dat");
  optParser->AddUsageLine("");

  optParser->AddFlag("help", "h");
  optParser->AddFlag("version", "v");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("list-parameters");
  optParser->AddFlag("sample-config");
  optParser->AddFlag("save-expanded");
  optParser->AddFlag("nosubsampling");
  optParser->AddFlag("nosave");
  optParser->AddOption("image-info");
  optParser->AddOption("output", "o");
  optParser->AddOption("estimation-size");
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

  /* Process the results: options */
  // First two are options which print useful info and then exit the program
  if ( optParser->FlagSet("help") || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("version") ) {
    printf("makemultimages version %s\n\n", VERSION_STRING);
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

  if (optParser->FlagSet("save-expanded")) {
    theOptions->saveExpandedImage = true;
  }
  if (optParser->FlagSet("nosubsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (optParser->FlagSet("nosave")) {
    theOptions->saveImage = false;
  }
  if (optParser->OptionSet("image-info")) {
    theOptions->imageInfoFile = optParser->GetTargetString("image-info");
    theOptions->imageInfoFile_exists = true;
  }
  if (optParser->OptionSet("output")) {
    theOptions->outputImageName = optParser->GetTargetString("output");
    theOptions->noOutputImageName = false;
  }
  if (optParser->OptionSet("estimation-size")) {
    if (NotANumber(optParser->GetTargetString("estimation-size").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** ERROR: estimation size should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->estimationImageSize = atol(optParser->GetTargetString("estimation-size").c_str());
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


/* END OF FILE: makemultiimages_main.cpp --------------------------------- */
