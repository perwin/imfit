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

// Solvers (optimization algorithms)
#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "diff_evoln_fit.h"
#include "nmsimplex_fit.h"

#include "commandline_parser.h"
#include "config_file_parser.h"
#include "print_results.h"


/* ---------------- Definitions & Constants ----------------------------- */
#define MPFIT_SOLVER        1
#define DIFF_EVOLN_SOLVER   2
#define NMSIMPLEX_SOLVER    3

#define DEFAULT_FTOL	1.0e-8
#define MAX_DE_GENERATIONS	600

#define NO_MAGNITUDES  -10000.0   /* indicated data are *not* in magnitudes */
//#define MONTE_CARLO_ITER   100
//#define BOOTSTRAP_ITER     200

#define DEFAULT_CONFIG_FILE   "imfit_config.dat"
#define DEFAULT_OUTPUT_PARAMETER_FILE   "bestfit_parameters_imfit.dat"


// Option names for use in config files
static string  kGainString = "GAIN";
static string  kReadNoiseString = "READNOISE";
static string  kNCombinedString = "NCOMBINED";
static string  kOriginalSkyString = "ORIGINAL_SKY";


#define VERSION_STRING      "0.9.8"



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
  bool  saveModel;
  std::string  outputModelFileName;
  bool  saveResidualImage;
  std::string  outputResidualFileName;
  bool  saveWeightImage;
  std::string  outputWeightFileName;
  bool  saveBestFitParams;
  std::string  outputParameterFileName;
  bool  useImageHeader;
  double  gain;
  bool  gainSet;
  double  readNoise;
  bool  readNoiseSet;
  int  nCombined;
  bool  nCombinedSet;
  double  originalSky;
  bool  originalSkySet;
  double  ftol;
  bool  ftolSet;
  char  modelName[MAXLINE];
  bool  noModel;
  char  paramString[MAXLINE];
  bool  newParameters;
  bool  verbose;
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
  int  maxThreads;
  bool  maxThreadsSet;
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
void PrepareImageComments( vector<string> *comments, const string &programName, 
                           commandOptions *mainOptions );


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
  vector<string>  imageCommentsList;
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
  options.saveModel = false;
  options.saveResidualImage = false;
  options.saveWeightImage = false;
  options.saveBestFitParams = true;
  options.outputParameterFileName = DEFAULT_OUTPUT_PARAMETER_FILE;
  options.useImageHeader= false;
  options.gain = 1.0;
  options.gainSet = false;
  options.readNoise = 0.0;
  options.readNoiseSet = false;
  options.nCombined = 1;
  options.nCombinedSet = false;
  options.originalSky = 0.0;
  options.originalSkySet = false;
  options.ftol = DEFAULT_FTOL;
  options.ftolSet = false;
  options.noModel = true;
  options.noParamLimits = true;
  options.newParameters = false;
  options.magZeroPoint = NO_MAGNITUDES;
  options.printChiSquaredOnly = false;
  options.printImages = false;
  options.solver = MPFIT_SOLVER;
  options.verbose = true;

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
  // (for this and other image-access, we rely on the cfitsio library to catch errors 
  // like nonexistent files)
  printf("Reading data image (\"%s\") ...\n", options.imageFileName.c_str());
  allPixels = ReadImageAsVector(options.imageFileName, &nColumns, &nRows);
  // Reminder: nColumns = n_pixels_per_row = x-size; nRows = n_pixels_per_column = y-size
  nPixels_tot = nColumns * nRows;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %d\n", 
           nColumns, nRows, nPixels_tot);
  // Determine X0,Y0 pixel offset values if user specified an image section
  DetermineImageOffset(options.imageFileName, &X0_offset, &Y0_offset);

  /* Get and check mask image */
  if (options.maskImagePresent) {
    printf("Reading mask image (\"%s\") ...\n", options.maskFileName.c_str());
    allMaskPixels = ReadImageAsVector(options.maskFileName, &nMaskColumns, &nMaskRows);
    if ((nMaskColumns != nColumns) || (nMaskRows != nRows)) {
      fprintf(stderr, "\n*** WARNING: Dimensions of mask image (%s: %d columns, %d rows)\n",
             options.maskFileName.c_str(), nMaskColumns, nMaskRows);
      fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      return -1;
    }
    maskAllocated = true;
  }
           
  /* Get and check error image, if supplied */
  if (options.noiseImagePresent) {
    printf("Reading noise image (\"%s\") ...\n", options.noiseFileName.c_str());
    allErrorPixels = ReadImageAsVector(options.noiseFileName, &nErrColumns, &nErrRows);
    errorPixels_allocated = true;
    if ((nErrColumns != nColumns) || (nErrRows != nRows)) {
      fprintf(stderr, "\n*** WARNING: Dimensions of error image (%s: %d columns, %d rows)\n",
             noiseImage.c_str(), nErrColumns, nErrRows);
      fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      return -1;
    }
  }
  else
    printf("* No noise image supplied ... will generate noise image from input image.\n");
  
  /* Read in PSF image, if supplied */
  if (options.psfImagePresent) {
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
  // Put limits on number of FFTW and OpenMP threads, if requested by user
  if (options.maxThreadsSet)
    theModel->SetMaxThreads(options.maxThreads);

  /* Add functions to the model object */
  status = AddFunctions(theModel, functionList, functionSetIndices, options.subsamplingFlag);
  if (status < 0) {
  	fprintf(stderr, "*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know total # parameters
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
  
  
  // Add PSF image vector, if present (needs to be added prior to image data, so that
  // ModelObject can figure out proper internal model-image size
  if (options.psfImagePresent)
    theModel->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels);

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


  
  /* START OF MINIMIZATION-ROUTINE-RELATED CODE */
  // Parameter limits and other info:
  // First we create a C-style array of mp_par structures, containing parameter constraints
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
    	//printf("Fixed parameter detected (i = %d)\n", i);
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
  if (options.printChiSquaredOnly) {
    printf("\n");
//    theModel->SetupChisquaredCalcs();
    status = 1;
    PrintResults(paramsVect, 0, 0, theModel, nFreeParams, parameterInfo, status);
    printf("\n");
  }
  else {
    // DO THE FIT!
    if (options.solver == MPFIT_SOLVER) {
      bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero results structure */
      mpfitResult.xerror = paramErrs;
      bzero(&mpConfig, sizeof(mpConfig));
      mpConfig.maxiter = 1000;
      if (options.ftolSet)
        mpConfig.ftol = options.ftol;
      if (options.verbose)
        mpConfig.verbose = 1;
      else
        mpConfig.verbose = 0;
      printf("\nCalling mpfit ...\n");
      status = mpfit(myfunc, nPixels_tot, nParamsTot, paramsVect, mpfitParameterConstraints,
      							&mpConfig, theModel, &mpfitResult);
      PrintResults(paramsVect, 0, &mpfitResult, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
    }
    else if (options.solver == DIFF_EVOLN_SOLVER) {
      printf("\nCalling DiffEvolnFit ..\n");
      status = DiffEvolnFit(nParamsTot, paramsVect, parameterInfo, theModel, 
      						MAX_DE_GENERATIONS, options.ftol);
      printf("\n");
      PrintResults(paramsVect, 0, 0, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
    }
    else if (options.solver == NMSIMPLEX_SOLVER) {
      printf("\nCalling Nelder-Mead Simplex solver ..\n");
      status = NMSimplexFix(nParamsTot, paramsVect, parameterInfo, theModel, options.ftol);
      printf("\n");
      PrintResults(paramsVect, 0, 0, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
    }
  }


  // TESTING (remove later)
  if (options.printImages)
    theModel->PrintModelImage();

  // Handle assorted output requests
  if (options.saveBestFitParams) {
    printf("Saving best-fit parameters in file \"%s\"\n", options.outputParameterFileName.c_str());
    SaveParameters(paramsVect, theModel, parameterInfo, options.outputParameterFileName,
    								argc, argv);
  }
  if (options.saveModel) {
    string  progName = "imfit ";
    progName += VERSION_STRING;
    PrepareImageComments(&imageCommentsList, progName, &options);
    printf("Saving model image in file \"%s\"\n", options.outputModelFileName.c_str());
    SaveVectorAsImage(theModel->GetModelImageVector(), options.outputModelFileName, 
                      nColumns, nRows, imageCommentsList);
  }
  if (options.saveResidualImage) {
    printf("Saving residual (input - model) image in file \"%s\"\n", options.outputResidualFileName.c_str());
    SaveVectorAsImage(theModel->GetResidualImageVector(), options.outputResidualFileName, 
                      nColumns, nRows, imageCommentsList);
  }
  if (options.saveWeightImage) {
    printf("Saving weight image in file \"%s\"\n", options.outputWeightFileName.c_str());
    SaveVectorAsImage(theModel->GetWeightImageVector(), options.outputWeightFileName, 
                      nColumns, nRows, imageCommentsList);
  }


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
  
  printf("Done!\n\n");
  
  return 0;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   imfit [options] imagefile.fits");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine(" -v  --version                Prints version number");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --list-parameters        Prints list of parameter names for each available function");
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" -c  --config <config-file>   configuration file");
  optParser->AddUsageLine("     --chisquare-only         Print chi^2 of input model and quit");
  optParser->AddUsageLine("     --noise <noisemap.fits>  Noise image to use");
  optParser->AddUsageLine("     --mask <mask.fits>       Mask image to use");
  optParser->AddUsageLine("     --psf <psf.fits>         PSF image to use");
  optParser->AddUsageLine("     --nosubsampling          Do *not* do pixel subsampling near centers");
  optParser->AddUsageLine("     --save-params <output-file>          Save best-fit parameters in config-file format [default = bestfit_parameters_imfit.dat]");
  optParser->AddUsageLine("     --save-model <outputname.fits>       Save best-fit model image");
  optParser->AddUsageLine("     --save-residual <outputname.fits>       Save residual (input - model) image");
  optParser->AddUsageLine("     --save-weights <outputname.fits>       Save weight image");
//  optParser->AddUsageLine("     --use-headers            Use image header values for gain, readnoise [NOT IMPLEMENTED YET]");
  optParser->AddUsageLine("     --sky <sky-level>        Original sky background (ADUs) which was subtracted from image");
  optParser->AddUsageLine("     --gain <value>           Image gain (e-/ADU)");
  optParser->AddUsageLine("     --readnoise <value>      Image read noise (e-)");
  optParser->AddUsageLine("     --ncombined <value>      Number of images averaged to make final image (if counts are average or median)");
  optParser->AddUsageLine("     --errors-are-variances   Indicates that values in noise image = variances (instead of sigmas)");
  optParser->AddUsageLine("     --errors-are-weights     Indicates that values in noise image = weights (instead of sigmas)");
  optParser->AddUsageLine("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --ftol                   Fractional tolerance in chi^2 for convergence [default = 1.0e-8]");
  optParser->AddUsageLine("     --de                     Use differential evolution solver instead of L-M");
  optParser->AddUsageLine("     --nm                     Use Nelder-Mead simplex solver instead of L-M");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --quiet                  Turn off printing of mpfit iteration updates");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --max-threads <int>      Maximum number of threads to use");
//  optParser->AddUsageLine("     --printimage             Print out images (for debugging)");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("EXAMPLES:");
  optParser->AddUsageLine("   imfit -c model_config_a.dat ngc100.fits");
  optParser->AddUsageLine("   imfit -c model_config_b.dat ngc100.fits[405:700,844:1060] --mask ngc100_mask.fits[405:700,844:1060] --gain 4.5 --readnoise 0.7");


  /* by default all options are checked on the command line and from option/resource file */
  optParser->AddFlag("help", "h");
  optParser->AddFlag("version", "v");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("list-parameters");
  optParser->AddFlag("printimage");
  optParser->AddFlag("chisquare-only");
  optParser->AddFlag("use-headers");
  optParser->AddFlag("errors-are-variances");
  optParser->AddFlag("errors-are-weights");
  optParser->AddFlag("mask-zero-is-bad");
  optParser->AddFlag("nosubsampling");
  optParser->AddFlag("de");
  optParser->AddFlag("nm");
  optParser->AddFlag("quiet");
  optParser->AddOption("noise");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("mask");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("psf");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("save-params");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("save-model");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("save-residual");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("save-weights");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("sky");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("gain");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("readnoise");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("ncombined");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("ftol");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("config", "c");
  optParser->AddOption("max-threads");      /* an option (takes an argument), supporting only long form */

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

  if (optParser->FlagSet("printimage")) {
    theOptions->printImages = true;
  }
  if (optParser->FlagSet("chisquare-only")) {
    printf("\t No fitting will be done!\n");
    theOptions->printChiSquaredOnly = true;
  }
  if (optParser->FlagSet("de")) {
  	printf("\t Differential Evolution selected!\n");
  	theOptions->solver = DIFF_EVOLN_SOLVER;
  }
  if (optParser->FlagSet("nm")) {
  	printf("\t Nelder-Mead simplex solver selected!\n");
  	theOptions->solver = NMSIMPLEX_SOLVER;
  }
  if (optParser->FlagSet("nosubsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (optParser->FlagSet("quiet")) {
    theOptions->verbose = false;
  }
  if (optParser->FlagSet("use-headers")) {
    theOptions->useImageHeader = true;
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
    printf("\toutput residual (input - model) image = %s\n", theOptions->outputResidualFileName.c_str());
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
      fprintf(stderr, "*** WARNING: sky should be a real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->originalSky = atof(optParser->GetTargetString("sky").c_str());
    theOptions->originalSkySet = true;
    printf("\toriginal sky level = %g ADU\n", theOptions->originalSky);
  }
  if (optParser->OptionSet("gain")) {
    if (NotANumber(optParser->GetTargetString("gain").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: gain should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->gain = atof(optParser->GetTargetString("gain").c_str());
    theOptions->gainSet = true;
    printf("\tgain = %g e-/ADU\n", theOptions->gain);
  }
  if (optParser->OptionSet("readnoise")) {
    if (NotANumber(optParser->GetTargetString("readnoise").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: read noise should be a non-negative real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->readNoise = atof(optParser->GetTargetString("readnoise").c_str());
    theOptions->readNoiseSet = true;
    printf("\tread noise = %g e-\n", theOptions->readNoise);
  }
  if (optParser->OptionSet("ncombined")) {
    if (NotANumber(optParser->GetTargetString("ncombined").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: ncombined should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->nCombined = atoi(optParser->GetTargetString("ncombined").c_str());
    theOptions->nCombinedSet = true;
    printf("\tn_combined = %d\n", theOptions->nCombined);
  }
  if (optParser->OptionSet("ftol")) {
    if (NotANumber(optParser->GetTargetString("ftol").c_str(), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: ftol should be a positive real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->ftol = atof(optParser->GetTargetString("ftol").c_str());
    theOptions->ftolSet = true;
    printf("\tfractional tolerance ftol for chi^2 convergence = %g\n", theOptions->ftol);
  }
  if (optParser->OptionSet("max-threads")) {
    if (NotANumber(optParser->GetTargetString("max-threads").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: max-threads should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->maxThreads = atol(optParser->GetTargetString("max-threads").c_str());
    theOptions->maxThreadsSet = true;
  }

  delete optParser;

}



// Note that we only use options from the config file if they have *not*
// already been set by the command line (i.e., command-line options override
// config-file values).
void HandleConfigFileOptions( configOptions *configFileOptions, commandOptions *mainOptions )
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
        printf("Value from config file: gain = %f\n", newDblVal);
        mainOptions->gain = newDblVal;
      }
      continue;
    }
    if (configFileOptions->optionNames[i] == kReadNoiseString) {
      if (mainOptions->readNoiseSet) {
        printf("Read-noise value in config file ignored (using command-line value)\n");
      } else {
        newDblVal = strtod(configFileOptions->optionValues[i].c_str(), NULL);
        printf("Value from config file: read noise = %f\n", newDblVal);
        mainOptions->readNoise = newDblVal;
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



void PrepareImageComments( vector<string> *comments, const string &programName, 
                           commandOptions *mainOptions )
{
  string  aString;
  char  *my_string;
  
  asprintf(&my_string, "Image generated by %s", programName.c_str());
  aString = my_string;
  comments->push_back(string(my_string));
  asprintf(&my_string, "Using parameters saved in file: %s", mainOptions->outputParameterFileName.c_str());
  comments->push_back(string(my_string));
  if (mainOptions->psfImagePresent) {
    asprintf(&my_string, "Convolved with PSF image: %s", mainOptions->psfFileName.c_str());
    comments->push_back(string(my_string));
  }
}




/* END OF FILE: imfit_main.cpp ------------------------------------------- */
