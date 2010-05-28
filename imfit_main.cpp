/* FILE: imfit_main.cpp -------------------------------------------------- */
/*
 * 
 * The proper translations are:
 * NAXIS1 = naxes[0] = nColumns = sizeX;
 * NAXIS2 = naxes[1] = nRows = sizeY.
 *
 *
 * HISTORY
 *     3--5 Dec 2009 [v0.2]: Added handling of mask images and parameter limits.
 *    10 Nov--2 Dec 2009: Early stages of developement
*/



/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>
#include <vector>

#include "definitions.h"
#include "utilities_pub.h"
#include "image_io.h"
#include "model_object.h"
#include "function_object.h"
#include "add_functions.h"
#include "param_struct.h"   // for mp_par structure
#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "diff_evoln_fit.h"
#include "anyoption.h"   // Kishan Thomas' class for command-line option parsing
#include "config_file_parser.h"
//#include "statistics.h"
#include "print_results.h"


/* ---------------- Definitions & Constants ----------------------------- */
#define MPFIT_SOLVER        1
#define DIFF_EVOLN_SOLVER   2

#define MAX_DE_GENERATIONS	600

#define NO_MAGNITUDES  -10000.0   /* indicated data are *not* in magnitudes */
#define MONTE_CARLO_ITER   100
#define BOOTSTRAP_ITER     200

#define CMDLINE_ERROR1 "Usage: -p must be followed by a string containing initial parameter values for the model"
#define CMDLINE_ERROR2 "Usage: -l must be followed by a filename for a file containing parameter limits"
#define CMDLINE_ERROR3 "For differential-evolution solver, a file containing parameter limits *must* be supplied (\"-l\" option)"
#define CMDLINE_ERROR4 "Usage: --mcoffset must be followed by a positive number"
#define CMDLINE_ERROR5 "Usage: --magzp must be followed by a positive number"
#define CMDLINE_ERROR6 "Usage: --mciter must be followed by a positive integer"
#define CMDLINE_ERROR7 "Usage: --pf must be followed by a positive real number"

#define FITS_FILENAME   "testimage_expdisk_tiny.fits"
#define FITS_ERROR_FILENAME   "tiny_uniform_image_0.1.fits"
#define DEFAULT_CONFIG_FILE   "imfit_config.dat"

#define VERSION_STRING      " v0.5"



typedef struct {
  std::string  configFileName;
  std::string  imageFileName;
  bool  noImage;
  std::string  psfFileName;
  bool  psfImagePresent;
  std::string  noiseFileName;
  bool  noiseImagePresent;
  int  errorType;
  std::string  maskFileName;
  bool  maskImagePresent;
  int  maskFormat;
  bool  subsamplingFlag;
  std::string  outputModelFileName;
  bool  outputModel;
  bool  useImageHeader;
  double  gain;
  double  readNoise;
  int  nCombined;
  double  originalSky;
  char  modelName[MAXLINE];
  bool  noModel;
  char  paramString[MAXLINE];
  bool  newParameters;
//  bool  doBootstrap;
//  int  bootstrapIterations;
//  bool  doMonteCarlo;
//  int  mcIterations;
//  double  monteCarloOffset;
  double  magZeroPoint;
  bool  noParamLimits;
  bool  printImages;
  bool printChiSquaredOnly;
  int  solver;
} commandOptions;



/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void DetermineImageOffset( const std::string &fullImageName, double *x_offset,
					double *y_offset);
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
void HandleConfigFileOptions( configOptions *configFileOptions, 
															commandOptions *mainOptions );
int myfunc( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */




/* This is the function used by mpfit() to compute the vector of deviates.
 * In our case, it's a wrapper which tells the ModelObject to compute the
 * deviates.
 */
int myfunc( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *theModel )
{

  theModel->ComputeDeviates(deviates, params);
  return 0;
}



/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  int  nPixels_tot, nColumns, nRows;
  int  nPixels_psf, nRows_psf, nColumns_psf;
  int  nErrColumns, nErrRows, nMaskColumns, nMaskRows;
  int  nDegFreedom;
  int  nParamsTot, nFreeParams;
  double  *allPixels;
  double  *psfPixels;
  double  *allErrorPixels;
  bool  errorPixels_allocated = false;
  double  *allMaskPixels;
  bool  maskAllocated = false;
  double  *paramsVect;
  double  *paramErrs;
  double  X0_offset = 0.0;
  double  Y0_offset = 0.0;
  std::string  noiseImage;
  std::string  mpfitMessage;
  std::string  baseFileName;
  ModelObject  *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  functionSetIndices;
  bool  paramLimitsExist = false;
  bool  parameterInfo_allocated = false;
  mp_result  mpfitResult;
  mp_config  mpConfig;
  mp_par  *parameterInfo;
  mp_par  *mpfitParameterConstraints;
  int  status;
  commandOptions  options;
  configOptions  userConfigOptions;
  const std::string  X0_string("X0");
  const std::string  Y0_string("Y0");
  
  
  /* Process the command line */
  /* First, set up the options structure: */
  options.configFileName = DEFAULT_CONFIG_FILE;
  options.noImage = true;
  options.psfImagePresent = false;
  options.noiseImagePresent = false;
  options.errorType = WEIGHTS_ARE_SIGMAS;
  options.maskImagePresent = false;
  options.maskFormat = MASK_ZERO_IS_GOOD;
  options.subsamplingFlag = true;
  options.outputModel = false;
  options.useImageHeader= false;
  options.gain = 1.0;
  options.readNoise = 0.0;
  options.nCombined = 1;
  options.originalSky = 0.0;
  options.noModel = true;
  options.noParamLimits = true;
  options.newParameters = false;
//  options.doBootstrap = false;
//  options.bootstrapIterations = BOOTSTRAP_ITER;
//  options.doMonteCarlo = false;
//  options.monteCarloOffset = 0.0;
//  options.mcIterations = 0;
  options.magZeroPoint = NO_MAGNITUDES;
  options.printChiSquaredOnly = false;
  options.printImages = false;
  options.solver = MPFIT_SOLVER;

  ProcessInput(argc, argv, &options);


  /* Read configuration file */
  if (! FileExists(options.configFileName.c_str())) {
    fprintf(stderr, "\n*** WARNING: Unable to find configuration file \"%s\"!\n\n", 
           options.configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options.configFileName, true, functionList, parameterList, 
  								paramLimits, functionSetIndices, paramLimitsExist, userConfigOptions);
  if (status != 0) {
    fprintf(stderr, "\n*** WARNING: Failure reading configuration file!\n\n");
    return -1;
  }

  // Parse and process user-supplied (non-function) values from config file, if any
  HandleConfigFileOptions(&userConfigOptions, &options);
  
  if (options.noImage) {
    fprintf(stderr, "*** WARNING: No image to fit!\n\n");
    return -1;
  }

  /* Get image data and sizes */
  // Note that we rely on the cfitsio library to catch errors like nonexistent files
  printf("Reading data image (\"%s\") ...\n", options.imageFileName.c_str());
  allPixels = ReadImageAsVector(options.imageFileName, &nColumns, &nRows);
  // Reminder: nColumns = n_pixels_per_row = x-size
  // Reminder: nRows = n_pixels_per_column = y-size
  nPixels_tot = nColumns * nRows;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %d\n", 
           nColumns, nRows, nPixels_tot);
  // Determine X0,Y0 pixel offset values if user specified an image section
  DetermineImageOffset(options.imageFileName, &X0_offset, &Y0_offset);

  /* Get and check mask image */
  if (options.maskImagePresent) {
    // Note that we rely on the cfitsio library to catch errors like nonexistent files
    printf("Reading mask image (\"%s\") ...\n", options.maskFileName.c_str());
    allMaskPixels = ReadImageAsVector(options.maskFileName, &nMaskColumns, &nMaskRows);
    if ((nMaskColumns != nColumns) || (nMaskRows != nRows)) {
      fprintf(stderr, "\n*** WARNING: Dimenstions of mask image (%s: %d columns, %d rows)\n",
             options.maskFileName.c_str(), nMaskColumns, nMaskRows);
      fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      return -1;
    }
    maskAllocated = true;
  }
           
  /* Get and check error image (or else tell function object to generate one) */
  if (options.noiseImagePresent) {
    // Note that we rely on the cfitsio library to catch errors like nonexistent files
    printf("Reading noise image (\"%s\") ...\n", options.noiseFileName.c_str());
    allErrorPixels = ReadImageAsVector(options.noiseFileName, &nErrColumns, &nErrRows);
    errorPixels_allocated = true;
    if ((nErrColumns != nColumns) || (nErrRows != nRows)) {
      fprintf(stderr, "\n*** WARNING: Dimenstions of error image (%s: %d columns, %d rows)\n",
             noiseImage.c_str(), nErrColumns, nErrRows);
      fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      return -1;
    }
  }
  else
    printf("* No noise image supplied ... will generate noise image from input image\n");
  
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


  /* Create the model object */
  theModel = new ModelObject();
  /* Add functions to the model object */
  status = AddFunctions(theModel, functionList, functionSetIndices, options.subsamplingFlag);
  if (status < 0) {
  	fprintf(stderr, "*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know how total # parameters
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  if (nParamsTot != (int)parameterList.size()) {
  	fprintf(stderr, "*** WARNING: number of input parameters (%d) does not equal", 
  	       (int)parameterList.size());
  	fprintf(stderr, " required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	exit(-1);
  }
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  
  
  /* Add image data, errors, and mask to the model object */
  theModel->AddImageDataVector(allPixels, nColumns, nRows, options.nCombined);
  theModel->PrintDescription();
  if (options.printImages)
    theModel->PrintInputImage();
  
  if (options.noiseImagePresent)
    theModel->AddErrorVector(nPixels_tot, nColumns, nRows, allErrorPixels,
                             options.errorType);
  else
    theModel->GenerateErrorVector(options.gain, options.readNoise, options.originalSky);

  if (maskAllocated) {
    theModel->AddMaskVector(nPixels_tot, nColumns, nRows, allMaskPixels,
                             options.maskFormat);
    theModel->ApplyMask();
  }
  nDegFreedom = theModel->GetNValidPixels() - nFreeParams;

  // Add PSF image vector, if present
  if (options.psfImagePresent)
    theModel->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels);


  
  // Parameter limits and other info:
  // OK, first we create a C-style array of mp_par structures, containing parameter constraints
  // (if any) *and* any other useful info (like X0,Y0 offset values).  This will be used
  // by DiffEvolnFit (if called) and by PrintResults.  We also decrement nFreeParams for
  // each *fixed* parameter.
  // Then we point the mp_par-array variable mpfitParameterConstraints to this array *if*
  // there are actual parameter constraints; if not, mpfitParameterConstraints is set = NULL,
  // since that's what mpfit() expects when there are no constraints.
  printf("Setting up parameter information array ...\n");
  parameterInfo = (mp_par *) calloc((size_t)nParamsTot, sizeof(mp_par));
  parameterInfo_allocated = true;
  for (int i = 0; i < nParamsTot; i++) {
    parameterInfo[i].fixed = paramLimits[i].fixed;
    if (parameterInfo[i].fixed == 1) {
    	printf("Fixed parameter detected (i = %d)\n", i);
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
  if ((options.solver == MPFIT_SOLVER) && (! paramLimitsExist)) {
    // If parameters are unconstrained, then mpfit() expects a NULL mp_par array
    printf("No parameter constraints!\n");
    mpfitParameterConstraints = NULL;
  } else {
    mpfitParameterConstraints = parameterInfo;
  }
  
  
  /* Copy parameters into C array, correcting for X0,Y0 offsets */
  paramsVect = (double *) calloc(nParamsTot, sizeof(double));
  for (int i = 0; i < nParamsTot; i++) {
    if (theModel->GetParameterName(i) == X0_string) {
      paramsVect[i] = parameterList[i] - X0_offset;
    } else if (theModel->GetParameterName(i) == Y0_string) {
      paramsVect[i] = parameterList[i] - Y0_offset;
    } else
      paramsVect[i] = parameterList[i];
  }
  
  
  if (options.printChiSquaredOnly) {
    printf("\n");
    theModel->SetupChisquaredCalcs();
    status = 1;
    PrintResults(paramsVect, 0, 0, theModel, nFreeParams, parameterInfo, status);
    printf("\n");
  }
  else {
    /* DO THE FIT! */
    if (options.solver == MPFIT_SOLVER) {
      // Some trial mpfit experiments:
      bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero results structure */
      mpfitResult.xerror = paramErrs;
      bzero(&mpConfig, sizeof(mpConfig));
      mpConfig.maxiter = 1000;
      printf("\nCalling mpfit ...\n");
      status = mpfit(myfunc, nPixels_tot, nParamsTot, paramsVect, mpfitParameterConstraints, &mpConfig, 
  		  						theModel, &mpfitResult);
      PrintResults(paramsVect, 0, &mpfitResult, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
    }
    else {
      printf("\nCalling DiffEvolnFit ..\n");
      status = DiffEvolnFit(nParamsTot, paramsVect, parameterInfo, theModel, MAX_DE_GENERATIONS);
      printf("\n");
      PrintResults(paramsVect, 0, 0, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
    }
  }


  // TESTING (remove later)
  if (options.printImages)
    theModel->PrintModelImage();

  /* TEST: save best-fit model image under new name: */
  // LATER: honor the boolean variable options.outputModel, which should
  // tell us (if it's = false) *not* to save the model image
  if (options.outputModel)
    SaveVectorAsImage(theModel->GetModelImageVector(), options.outputModelFileName, 
                      nColumns, nRows);


  // Free up memory
  free(allPixels);       // allocated in ReadImageAsVector()
  if (errorPixels_allocated)
    free(allErrorPixels);  // allocated in ReadImageAsVector()
  if (options.psfImagePresent)
    free(psfPixels);
  if (maskAllocated)
    free(allMaskPixels);
  free(paramsVect);
  free(paramErrs);
  if (parameterInfo_allocated)
    free(parameterInfo);
  delete theModel;
  
  return 0;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  AnyOption *opt = new AnyOption();
  opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* SET THE USAGE/HELP   */
  opt->addUsage("Usage: ");
  opt->addUsage("   imfit [options] imagefile.fits");
  opt->addUsage(" -h  --help                   Prints this help");
  opt->addUsage("     --list-functions         Prints list of available functions (components)");
  opt->addUsage("     --chisquare-only         Print chi^2 of input model and quit");
  opt->addUsage("     --printimage             Print out images (for debugging)");
  opt->addUsage(" -c  --config <config-file>   configuration file");
  opt->addUsage("     --de                     Use differential evolution solver (instead of L-M) [NOT IMPLEMENTED YET]");
  opt->addUsage("     --noise <noisemap.fits>  Noise image");
  opt->addUsage("     --mask <mask.fits>       Mask image");
  opt->addUsage("     --psf <psf.fits>         PSF image");
  opt->addUsage("     --nosubsampling          Do *not* do pixel subsampling near centers");
  opt->addUsage("     --save-model <outputname.fits>       Save best-fit model image");
  opt->addUsage("     --use-headers            Use image header values for gain, readnoise [NOT IMPLEMENTED YET]");
  opt->addUsage("     --sky <sky-level>        Original sky background (ADUs)");
  opt->addUsage("     --gain <value>           Image gain (e-/ADU)");
  opt->addUsage("     --readnoise <value>      Image read noise (e-)");
  opt->addUsage("     --ncombined <value>      Number of images averaged to make final image (if counts are average or median)");
  opt->addUsage("     --errors-are-variances   Indicates that values in noise image = variances");
  opt->addUsage("     --errors-are-weights     Indicates that values in noise image = weights");
  opt->addUsage("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  opt->addUsage("");


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag('q');
  opt->setFlag("help", 'h');
  opt->setFlag("list-functions");
  opt->setFlag("printimage");
  opt->setFlag("chisquare-only");
  opt->setFlag("de");
  opt->setFlag("use-headers");
  opt->setFlag("errors-are-variances");
  opt->setFlag("errors-are-weights");
  opt->setFlag("mask-zero-is-bad");
  opt->setFlag("nosubsampling");
  opt->setOption("noise");      /* an option (takes an argument), supporting only long form */
  opt->setOption("mask");      /* an option (takes an argument), supporting only long form */
  opt->setOption("psf");      /* an option (takes an argument), supporting only long form */
  opt->setOption("save-model");      /* an option (takes an argument), supporting only long form */
  opt->setOption("sky");        /* an option (takes an argument), supporting only long form */
  opt->setOption("gain");        /* an option (takes an argument), supporting only long form */
  opt->setOption("readnoise");        /* an option (takes an argument), supporting only long form */
  opt->setOption("ncombined");        /* an option (takes an argument), supporting only long form */
  opt->setOption("config", 'c');

  /* parse the command line:  */
  opt->processCommandArgs( argc, argv );


  /* Process the results: actual arguments, if any: */
  if (opt->getArgc() > 0) {
    theOptions->imageFileName = opt->getArgv(0);
    theOptions->noImage = false;
    printf("\tImage file = %s\n", theOptions->imageFileName.c_str());
  }

  /* Process the results: options */
  // First two are options which print useful info and then exit the program
  if ( opt->getFlag('q') ) {
  	printf("\n\nimfit_main: q option found!\n\n");
    delete opt;
    exit(1);
  }
  if ( opt->getFlag("help") || opt->getFlag('h') || (! opt->hasOptions()) ) {
    opt->printUsage();
    delete opt;
    exit(1);
  }
  if (opt->getFlag("list-functions")) {
    PrintAvailableFunctions();
    delete opt;
    exit(1);
  }

  if (opt->getFlag("printimage")) {
    theOptions->printImages = true;
  }
  if (opt->getFlag("chisquare-only")) {
    printf("\t No fitting will be done!\n");
    theOptions->printChiSquaredOnly = true;
  }
  if (opt->getFlag("de")) {
  	printf("\t Differential Evolution selected!\n");
  	theOptions->solver = DIFF_EVOLN_SOLVER;
  }
  if (opt->getFlag("nosubsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (opt->getFlag("use-header")) {
    theOptions->useImageHeader = true;
  }
  if (opt->getFlag("errors-are-variances")) {
    theOptions->errorType = WEIGHTS_ARE_VARIANCES;
  }
  if (opt->getFlag("errors-are-weights")) {
    theOptions->errorType = WEIGHTS_ARE_WEIGHTS;
  }
  if (opt->getFlag("mask-zero-is-bad")) {
    theOptions->maskFormat = MASK_ZERO_IS_BAD;
  }
  if (opt->getValue("config") != NULL) {
    theOptions->configFileName = opt->getValue("config");
    printf("\tconfiguration file = %s\n", theOptions->configFileName.c_str());
  }
  if (opt->getValue("noise") != NULL) {
    theOptions->noiseFileName = opt->getValue("noise");
    theOptions->noiseImagePresent = true;
    printf("\tnoise image = %s\n", theOptions->noiseFileName.c_str());
  }
  if (opt->getValue("psf") != NULL) {
    theOptions->psfFileName = opt->getValue("psf");
    theOptions->psfImagePresent = true;
    printf("\tPSF image = %s\n", theOptions->psfFileName.c_str());
  }
  if (opt->getValue("mask") != NULL) {
    theOptions->maskFileName = opt->getValue("mask");
    theOptions->maskImagePresent = true;
    printf("\tmask image = %s\n", theOptions->maskFileName.c_str());
  }
  if (opt->getValue("save-model") != NULL) {
    theOptions->outputModelFileName = opt->getValue("save-model");
    theOptions->outputModel = true;
    printf("\toutput best-fit model image = %s\n", theOptions->outputModelFileName.c_str());
  }
  if (opt->getValue("sky") != NULL) {
    if (NotANumber(opt->getValue("sky"), 0, kAnyReal)) {
      fprintf(stderr, "*** WARNING: sky should be a real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->originalSky = atof(opt->getValue("sky"));
    printf("\toriginal sky level = %g ADU\n", theOptions->originalSky);
  }
  if (opt->getValue("gain") != NULL) {
    if (NotANumber(opt->getValue("gain"), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: gain should be a positive real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->gain = atof(opt->getValue("gain"));
    printf("\tgain = %g e-/ADU\n", theOptions->gain);
  }
  if (opt->getValue("readnoise") != NULL) {
    if (NotANumber(opt->getValue("readnoise"), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: read noise should be a non-negative real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->readNoise = atof(opt->getValue("readnoise"));
    printf("\tread noise = %g e-\n", theOptions->readNoise);
  }
  if (opt->getValue("ncombined") != NULL) {
    if (NotANumber(opt->getValue("ncombined"), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: ncombined should be a positive integer!\n");
      delete opt;
      exit(1);
    }
    theOptions->nCombined = atoi(opt->getValue("ncombined"));
    printf("\tn_combined = %d\n", theOptions->nCombined);
  }

  delete opt;

}


void HandleConfigFileOptions( configOptions *configFileOptions, commandOptions *mainOptions )
{
  if (configFileOptions->nOptions == 0)
    return;
  ;
}


/* Function which takes the user-supplied image filename and determines what,
 * if any, x0 and y0 pixel offsets are implied by any section specification
 * in the filename.  Note that offsets are always >= 0.
 */
void DetermineImageOffset( const std::string &fullImageName, double *x_offset,
					double *y_offset)
{
  int  xStart, yStart;

  GetPixelStartCoords(fullImageName, &xStart, &yStart);
  *x_offset = xStart - 1;
  *y_offset = yStart - 1;
}






/* END OF FILE: imfit_main.cpp ------------------------------------------- */
