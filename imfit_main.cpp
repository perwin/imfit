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
#include <string>
#include <vector>

#include "definitions.h"
#include "utilities_pub.h"
#include "image_io.h"
#include "model_object.h"
#include "function_object.h"
#include "add_functions.h"
#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "anyoption.h"   // Kishan Thomas' class for command-line option parsing
#include "config_file_parser.h"
#include "statistics.h"


/* ---------------- Definitions & Constants ----------------------------- */
#define MPFIT_SOLVER        1
#define DIFF_EVOLN_SOLVER   2

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

#define VERSION_STRING      " v0.2"



typedef struct {
  std::string  configFileName;
  std::string  imageFileName;
  bool  noImage;
  std::string  noiseFileName;
  bool  noiseImagePresent;
  int  errorType;
  std::string  maskFileName;
  bool  maskImagePresent;
  int  maskFormat;
  std::string  outputModelFileName;
  bool  outputModel;
  bool  useImageHeader;
  double  gain;
  double  readNoise;
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
  char  paramLimitsFileName[MAX_FILENAME_LENGTH];
  bool  noParamLimits;
//  double  penaltyFunctionAmp;
  bool  printImages;
} commandOptions;



/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
void PrintParam( string& paramName, double paramValue, double paramErr );
void PrintResult( double *x, double *xact, mp_result *result,
					ModelObject *model, int nFreeParameters );
int myfunc( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */




void PrintParam( string& paramName, double paramValue, double paramErr )
{
  printf("  %10s = %f +/- %f\n", paramName.c_str(), paramValue, paramErr);
}


/* Simple routine to print the fit results [taken & modified from Craig Markwardt's
   testmpfit.c] 
*/
void PrintResult( double *x, double *xact, mp_result *result, ModelObject *model,
									int nFreeParameters )
{
  int  i;
  int  nValidPixels = model->GetNValidPixels();
  int  nDegreesFreedom = nValidPixels - nFreeParameters;
  double  aic, bic, chiSquared;
  
  if ((x == 0) || (result == 0)) return;
  printf("  CHI-SQUARE = %f    (%d DOF)\n", 
	 result->bestnorm, result->nfunc-result->nfree);
  printf("  INITIAL CHI^2 = %f\n", result->orignorm);
  printf("        NPAR = %d\n", result->npar);
  printf("       NFREE = %d\n", result->nfree);
  printf("     NPEGGED = %d\n", result->npegged);
  printf("     NITER = %d\n", result->niter);
  printf("      NFEV = %d\n", result->nfev);
  printf("\n");
  aic = AIC_corrected(result->bestnorm, nFreeParameters, nValidPixels, 1);
  bic = BIC(result->bestnorm, nFreeParameters, nValidPixels, 1);
  printf("Reduced Chi^2 = %f\n", result->bestnorm / nDegreesFreedom);
  printf("AIC = %f, BIC = %f\n\n", aic, bic);
  
  if (xact) {
    for (i=0; i < result->npar; i++) {
      printf("  P[%d] = %f +/- %f     (ACTUAL %f)\n", 
	     i, x[i], result->xerror[i], xact[i]);
    }
  } else {
    for (i = 0; i < result->npar; i++) {
      PrintParam(model->GetParameterName(i), x[i], result->xerror[i]);
//      printf("  P[%d] = %f +/- %f\n", 
//	     i, x[i], result->xerror[i]);
    }
  }    
}


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
  int  nErrColumns, nErrRows, nMaskColumns, nMaskRows;
  int  nDegFreedom;
  int  nParamsTot, nFreeParams;
  double  *allPixels;
  double  *allErrorPixels;
  double  *allMaskPixels;
  bool  maskAllocated = false;
  double  *paramsVect;
  double  *paramErrs;
  std::string  noiseImage;
  ModelObject  *theModel;
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  functionSetIndices;
  bool  paramLimitsExist = false;
  mp_result  mpfitResult;
  mp_config  mpConfig;
  mp_par  *mpParamLimits;
  int  status;
  commandOptions  options;
  
  
  /* Process the command line */
  /* First, set up the options structure: */
  options.configFileName = DEFAULT_CONFIG_FILE;
  options.noImage = true;
  options.noiseImagePresent = false;
  options.errorType = WEIGHTS_ARE_SIGMAS;
  options.maskImagePresent = false;
  options.maskFormat = MASK_ZERO_IS_GOOD;
  options.outputModel = false;
  options.useImageHeader= false;
  options.gain = 1.0;
  options.readNoise = 0.0;
  options.originalSky = 0.0;
  options.noModel = true;
  options.noParamLimits = true;
  options.paramLimitsFileName[0] = '-';
  options.newParameters = false;
//  options.penaltyFunctionAmp = -1;
//  options.doBootstrap = false;
//  options.bootstrapIterations = BOOTSTRAP_ITER;
//  options.doMonteCarlo = false;
//  options.monteCarloOffset = 0.0;
//  options.mcIterations = 0;
  options.magZeroPoint = NO_MAGNITUDES;
  options.printImages = false;

  ProcessInput(argc, argv, &options);


  /* Read configuration file */
  if (! FileExists(options.configFileName.c_str())) {
    printf("\n*** WARNING: Unable to find or open configuration file \"%s\"!\n\n", 
           options.configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options.configFileName, functionList, parameterList, 
  								paramLimits, functionSetIndices, paramLimitsExist);
  if (status != 0) {
    printf("\n*** WARNING: Failure reading configuration file!\n\n");
    return -1;
  }

  if (options.noImage) {
    printf("*** WARNING: No image to fit!\n\n");
    return -1;
  }

  /* Get image data and sizes */
  if (! FileExists(options.imageFileName.c_str())) {
    printf("\n*** WARNING: Unable to find or open input image \"%s\"!\n\n", 
           options.imageFileName.c_str());
    return -1;
  }
    printf("Reading data image (\"%s\") ...\n", options.imageFileName.c_str());
  allPixels = ReadImageAsVector(options.imageFileName, &nColumns, &nRows);
  // Reminder: nColumns = n_pixels_per_row = x-size
  // Reminder: nRows = n_pixels_per_column = y-size
  nPixels_tot = nColumns * nRows;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %d\n", 
           nColumns, nRows, nPixels_tot);
  
  /* Get and check mask image */
  if (options.maskImagePresent) {
    if (! FileExists(options.maskFileName.c_str())) {
      printf("\n*** WARNING: Unable to find or open mask image \"%s\"!\n\n", 
             options.maskFileName.c_str());
      return -1;
    }
    printf("Reading mask image (\"%s\") ...\n", options.maskFileName.c_str());
    allMaskPixels = ReadImageAsVector(options.maskFileName, &nMaskColumns, &nMaskRows);
    if ((nMaskColumns != nColumns) || (nMaskRows != nRows)) {
      printf("\n*** WARNING: Dimenstions of mask image (%s: %d columns, %d rows)\n",
             options.maskFileName.c_str(), nMaskColumns, nMaskRows);
      printf("do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      return -1;
    }
    maskAllocated = true;
  }
           
  /* Get and check error image (or else tell function object to generate one) */
  if (options.noiseImagePresent) {
    if (! FileExists(options.noiseFileName.c_str())) {
      printf("\n*** WARNING: Unable to find or open noise image \"%s\"!\n\n", 
           options.noiseFileName.c_str());
      return -1;
    }
    printf("Reading noise image (\"%s\") ...\n", options.noiseFileName.c_str());
    allErrorPixels = ReadImageAsVector(options.noiseFileName, &nErrColumns, &nErrRows);
    if ((nErrColumns != nColumns) || (nErrRows != nRows)) {
      printf("\n*** WARNING: Dimenstions of error image (%s: %d columns, %d rows)\n",
             noiseImage.c_str(), nErrColumns, nErrRows);
      printf("do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
             options.imageFileName.c_str(), nColumns, nRows);
      return -1;
    }
  }
  else
    printf("* No noise image supplied ... will generate noise image from input image\n");
           

  /* Create the model object */
  theModel = new ModelObject();
  /* Add functions to the model object */
  status = AddFunctions(theModel, functionList, functionSetIndices);
  if (status < 0) {
  	printf("*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know how total # parameters
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  if (nParamsTot != (int)parameterList.size()) {
  	printf("*** WARNING: number of input parameters (%d) does not equal", 
  	       (int)parameterList.size());
  	printf(" required number of parameters for specified functions (%d)!\n\n",
  	       nParamsTot);
  	exit(-1);
  }
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  
  
  /* Add image data, errors, and mask to the model object */
  theModel->AddImageDataVector(nPixels_tot, nColumns, nRows, allPixels);
  theModel->GetDescription();
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


  // Some trial mpfit experiments:
  bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero results structure */
  mpfitResult.xerror = paramErrs;
  bzero(&mpConfig, sizeof(mpConfig));
  mpConfig.maxiter = 1000;
  
  // Parameter limits:
  if (paramLimitsExist) {
    printf("Setting up parameter limits ...\n");
    mpParamLimits = (mp_par *) calloc((size_t)nParamsTot, sizeof(mp_par));
    for (int i = 0; i < nParamsTot; i++) {
      mpParamLimits[i].fixed = paramLimits[i].fixed;
      if (mpParamLimits[i].fixed == 1)
        nFreeParams--;
      mpParamLimits[i].limited[0] = paramLimits[i].limited[0];
      mpParamLimits[i].limited[1] = paramLimits[i].limited[1];
      mpParamLimits[i].limits[0] = paramLimits[i].limits[0];
      mpParamLimits[i].limits[1] = paramLimits[i].limits[1];
    }
  } else {
    mpParamLimits = NULL;
  }
  
  /* Copy parameters into C array and start the fitting process */
  paramsVect = (double *) calloc(nParamsTot, sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];
  
  printf("\nStarting mpfit ...\n");
  status = mpfit(myfunc, nPixels_tot, nParamsTot, paramsVect, mpParamLimits, &mpConfig, 
  								theModel, &mpfitResult);
  
  printf("*** mpfit status = %d\n", status);
  status = AddFunctions(theModel, functionList, functionSetIndices);
  PrintResult(paramsVect, 0, &mpfitResult, theModel, nFreeParams);
  
  printf("\nTrue degrees of freedom = %d\n", nDegFreedom);

  // TESTING (remove later)
  if (options.printImages)
    theModel->PrintModelImage();

  /* TEST: save best-fit model image under new name: */
  // LATER: honor the boolean variable options.outputModel, which should
  // tell us (if it's = false) *not* to save the model image
  if (options.outputModel)
    SaveVectorAsImage(theModel->GetModelImageVector(), options.outputModelFileName, 
                      nColumns, nRows);
  else   // no special name for output
    SaveVectorAsImage(theModel->GetModelImageVector(), "test_best-fit_model.fits", 
                      nColumns, nRows);


  // Free up memory
  free(allPixels);       // allocated in ReadImageAsVector()
  free(allErrorPixels);  // allocated in ReadImageAsVector()
  if (maskAllocated)
    free(allMaskPixels);
  free(paramsVect);
  free(paramErrs);
  if (paramLimitsExist)
    free(mpParamLimits);
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
  opt->addUsage("     --printimage             Print out images (for debugging)");
  opt->addUsage(" -c  --config <config-file>   configuration file");
  opt->addUsage("     --noise <noisemap.fits>  Noise image");
  opt->addUsage("     --mask <mask.fits>       Mask image");
  opt->addUsage("     --save-model <outputname.fits>       Save best-fit model image");
  opt->addUsage("     --use-headers            Use image header values for gain, readnoise");
  opt->addUsage("     --sky <sky-level>        Original sky background (ADUs)");
  opt->addUsage("     --gain <value>           Image gain (e-/ADU)");
  opt->addUsage("     --readnoise <value>      Image read noise (e-)");
  opt->addUsage("     --errors-are-variances   Indicates that values in noise image = variances");
  opt->addUsage("     --errors-are-weights     Indicates that values in noise image = weights");
  opt->addUsage("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  opt->addUsage("");


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag("help", 'h');
  opt->setFlag("printimage");
  opt->setFlag("use-headers");
  opt->setFlag("errors-are-variances");
  opt->setFlag("errors-are-weights");
  opt->setFlag("mask-zero-is-bad");
  opt->setOption("noise");      /* an option (takes an argument), supporting only long form */
  opt->setOption("mask");      /* an option (takes an argument), supporting only long form */
  opt->setOption("save-model");      /* an option (takes an argument), supporting only long form */
  opt->setOption("sky");        /* an option (takes an argument), supporting only long form */
  opt->setOption("gain");        /* an option (takes an argument), supporting only long form */
  opt->setOption("readnoise");        /* an option (takes an argument), supporting only long form */
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
  if ( opt->getFlag("help") || opt->getFlag('h') || (! opt->hasOptions()) ) {
    opt->printUsage();
    delete opt;
    exit(1);
  }
  if (opt->getFlag("printimage")) {
    theOptions->printImages = true;
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
      printf("*** WARNING: sky should be a real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->originalSky = atof(opt->getValue("sky"));
    printf("\toriginal sky level = %g ADU\n", theOptions->originalSky);
  }
  if (opt->getValue("gain") != NULL) {
    if (NotANumber(opt->getValue("gain"), 0, kPosReal)) {
      printf("*** WARNING: gain should be a positive real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->gain = atof(opt->getValue("gain"));
    printf("\tgain = %g e-/ADU\n", theOptions->gain);
  }
  if (opt->getValue("readnoise") != NULL) {
    if (NotANumber(opt->getValue("readnoise"), 0, kPosReal)) {
      printf("*** WARNING: read noise should be a non-negative real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->readNoise = atof(opt->getValue("readnoise"));
    printf("\tread noise = %g e-\n", theOptions->readNoise);
  }

  delete opt;

}






/* END OF FILE: imfit_main.cpp ------------------------------------------- */
