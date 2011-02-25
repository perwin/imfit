// Code for timing generation of images and PSF convolution.
// This is a modification of makeimage_main.cpp designed to generate images
// repeatedly and time how fast the image-generation and convolution code
// works.
//
// It is meant to test possible speedups in image generation (e.g., making use
// of multiple cores) and convolution (e.g., varying the FFTW "wisdom"
// parameters).




/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <sys/time.h>   // for timing-related functions and structs

#include "definitions.h"
#include "image_io.h"
#include "model_object.h"
#include "add_functions.h"
#include "commandline_parser.h"
#include "config_file_parser.h"
#include "utilities_pub.h"


/* ---------------- Definitions ---------------------------------------- */
#define NO_MAGNITUDES  -10000.0   /* indicated data are *not* in magnitudes */

#define CMDLINE_ERROR1 "Usage: -p must be followed by a string containing initial parameter values for the model"

#define DEFAULT_OUTPUT_FILENAME   "modelimage.fits"
#define FITS_FILENAME   "testimage_expdisk_tiny.fits"
#define FITS_ERROR_FILENAME   "tiny_uniform_image_0.1.fits"

#define VERSION_STRING      " v0.1"


typedef struct {
  std::string  outputImageName;
  bool  noImageName;
  std::string  referenceImageName;
  bool  noRefImage;
  bool  subsamplingFlag;
  bool  noImageDimensions;
  std::string  psfFileName;
  bool  psfImagePresent;
  int  nColumns;
  int  nRows;
  bool  noConfigFile;
  std::string  configFileName;
  char  modelName[MAXLINE];
  bool  noModel;
  char  paramString[MAXLINE];
  bool  newParameters;
  double  magZeroPoint;
  char  paramLimitsFileName[MAX_FILENAME_LENGTH];
  bool  noParamLimits;
  bool  printImages;
  int  nIterations;
  int  debugLevel;
} commandOptions;



/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
void ProcessInput2( int argc, char *argv[], commandOptions *theOptions );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */





int main(int argc, char *argv[])
{
  int  nPixels_tot, nColumns, nRows;
  int  nPixels_psf, nRows_psf, nColumns_psf;
  int  nParamsTot;
  int  status;
  double  *allPixels;
  double  *psfPixels;
  bool  allPixels_allocated = false;
  double  *paramsVect;
  ModelObject  *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<int>  functionSetIndices;
  commandOptions  options;
  configOptions  userConfigOptions;
  // timing-related stuff
  struct timeval  timer_start, timer_end;
  double  microsecs, time_elapsed;
  
  
  /* Process command line and parse config file: */
  options.outputImageName = DEFAULT_OUTPUT_FILENAME;
  options.noImageName = true;
  options.noRefImage = true;
  options.psfImagePresent = false;
  options.subsamplingFlag = true;
  options.noImageDimensions = true;
  options.noModel = true;
  options.noParamLimits = true;
  options.paramLimitsFileName[0] = '-';
  options.newParameters = false;
  options.magZeroPoint = NO_MAGNITUDES;
  options.printImages = false;
  options.nIterations = 1;
  options.debugLevel = 0;

  ProcessInput(argc, argv, &options);

  /* Read configuration file */
  if (! FileExists(options.configFileName.c_str())) {
    fprintf(stderr, "\n*** WARNING: Unable to find configuration file \"%s\"!\n\n", 
           options.configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options.configFileName, true, functionList, parameterList,
  							functionSetIndices, userConfigOptions);
  if (status != 0) {
    fprintf(stderr, "\n*** WARNING: Failure reading configuration file!\n\n");
    return -1;
  }

  if ((options.noRefImage) and (options.noImageDimensions)) {
    fprintf(stderr, "\n*** WARNING: Insufficient image dimensions (or no reference image) supplied!\n\n");
    return -1;
  }

  /* Get image size from reference image, if necessary */
  if (options.noImageDimensions) {
    // Note that we rely on the cfitsio library to catch errors like nonexistent files
    allPixels = ReadImageAsVector(options.referenceImageName, &nColumns, &nRows);
    allPixels_allocated = true;
    // Reminder: nColumns = n_pixels_per_row
    // Reminder: nRows = n_pixels_per_column
    printf("Reference image read: naxis1 [# rows] = %d, naxis2 [# columns] = %d\n",
           nRows, nColumns);
  }
  else {
    nColumns = options.nColumns;
    nRows = options.nRows;
  }
  nPixels_tot = nColumns * nRows;
  

  /* Read in PSF image, if supplied */
  if (options.psfImagePresent) {
    // Note that we rely on the cfitsio library to catch errors like nonexistent files
    printf("Reading PSF image (\"%s\") ...\n", options.psfFileName.c_str());
    psfPixels = ReadImageAsVector(options.psfFileName, &nColumns_psf, &nRows_psf);
    nPixels_psf = nColumns_psf * nRows_psf;
    printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %d\n", 
           nColumns_psf, nRows_psf, nPixels_psf);
  }
  else
    printf("* No PSF image supplied -- no image convolution will be done!\n");

  if (! options.subsamplingFlag)
    printf("* Pixel subsampling has been turned OFF.\n");

  /* Set up the model object */
  theModel = new ModelObject();
  theModel->SetDebugLevel(options.debugLevel);
  
  /* Add functions to the model object; also tells model object where function
     sets start */
  status = AddFunctions(theModel, functionList, functionSetIndices, options.subsamplingFlag);
  if (status < 0) {
  	fprintf(stderr, "*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
  }

  
  // Set up parameter vector(s), now that we know how many total parameters
  // there will be
  nParamsTot = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  if (nParamsTot != (int)parameterList.size()) {
  	fprintf(stderr, "*** WARNING: number of input parameters (%d) does not equal", 
  	       (int)parameterList.size());
  	fprintf(stderr, " required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	exit(-1);
 }
    
  /* Define the size of the requested model image */
  theModel->SetupModelImage(nPixels_tot, nColumns, nRows);
  theModel->PrintDescription();

  // Add PSF image vector, if present
  if (options.psfImagePresent)
    theModel->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels);


  /* Copy parameters into C array and generate the model image */
  paramsVect = (double *) calloc(nParamsTot, sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];
  
  
  // Generate the image (including convolution, if requested, repeatedly
  gettimeofday(&timer_start, NULL);
  for (int ii = 0; ii < options.nIterations; ii++)
    theModel->CreateModelImage(paramsVect);
  gettimeofday(&timer_end, NULL);
  microsecs = timer_end.tv_usec - timer_start.tv_usec;
  time_elapsed = timer_end.tv_sec - timer_start.tv_sec + microsecs/1e6;
  printf("\nELAPSED TIME FOR %d ITERATIONS: %.6f sec\n", options.nIterations, time_elapsed);
  printf("Mean time per iteration = %.7f\n", time_elapsed/options.nIterations);

  

  /* Save model image: */
  printf("\nSaving output model image (\"%s\") ...\n", options.outputImageName.c_str());
  vector<string> commentStrings;
  SaveVectorAsImage(theModel->GetModelImageVector(), options.outputImageName, 
                    nColumns, nRows, commentStrings);

  printf("Done!\n\n");


  // Free up memory
  if (allPixels_allocated)
    free(allPixels);
  if (options.psfImagePresent)
    free(psfPixels);
  free(paramsVect);
  delete theModel;
  
  return 0;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   makeimage [options] config-file");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine(" -o  --output <output-image.fits>        name for output image");
  optParser->AddUsageLine("     --refimage <reference-image.fits>   reference image (for image size)");
  optParser->AddUsageLine("     --psf <psf.fits>         PSF image (for convolution)");
  optParser->AddUsageLine("     --ncols <number-of-columns>   x-size of output image");
  optParser->AddUsageLine("     --nrows <number-of-rows>   y-size of output image");
  optParser->AddUsageLine("     --nosubsampling          Do *not* do pixel subsampling near centers");
  optParser->AddUsageLine("     --niterations <n>             number of iterations to do");
  optParser->AddUsageLine("     --debug <n>             debugging level");
  optParser->AddUsageLine("");


  /* by default all options are checked on the command line and from option/resource file */
  optParser->AddFlag("help", "h");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("nosubsampling");
  optParser->AddOption("output", "o");      /* an option (takes an argument) */
  optParser->AddOption("niterations");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("ncols");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("nrows");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("refimage");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("psf");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("debug");      /* an option (takes an argument), supporting only long form */

  /* parse the command line:  */
  optParser->ParseCommandLine( argc, argv );


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
  if (optParser->FlagSet("list-functions")) {
    PrintAvailableFunctions();
    delete optParser;
    exit(1);
  }

  if (optParser->FlagSet("nosubsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (optParser->OptionSet("output")) {
    theOptions->outputImageName = optParser->GetTargetString("output");
    theOptions->noImageName = false;
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
  if (optParser->OptionSet("niterations")) {
    if (NotANumber(optParser->GetTargetString("niterations").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: niterations should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nIterations = atol(optParser->GetTargetString("niterations").c_str());
  }
  if (optParser->OptionSet("ncols")) {
    if (NotANumber(optParser->GetTargetString("ncols").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: ncols should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nColumns = atol(optParser->GetTargetString("ncols").c_str());
  }
  if (optParser->OptionSet("nrows")) {
    if (NotANumber(optParser->GetTargetString("nrows").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: nrows should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nRows = atol(optParser->GetTargetString("nrows").c_str());
  }
  if (optParser->OptionSet("debug")) {
    if (NotANumber(optParser->GetTargetString("debug").c_str(), 0, kAnyInt)) {
      fprintf(stderr, "*** WARNING: debug should be an integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->debugLevel = atol(optParser->GetTargetString("debug").c_str());
  }

  if ((theOptions->nColumns) && (theOptions->nRows))
    theOptions->noImageDimensions = false;
  
  delete optParser;

}
