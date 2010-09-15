/* FILE: imfit1d_main.cpp  ----------------------------------------------- */
/*
 * This is a modified version of imfit_main.cpp, which does fitting of 1-D
 * profiles, where the profile is stored in a text file with two or three
 * columns of numbers (first column = radius or x value; second column = 
 * intensity, magnitudes per square arcsec, or some other y value; third
 * column = optional errors on y values).
*/



/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "definitions.h"
#include "utilities_pub.h"
#include "read_profile_pub.h"
#include "model_object.h"
#include "model_object_1d.h"
#include "add_functions_1d.h"
#include "function_object.h"
#include "func1d_exp.h"
#include "param_struct.h"   // for mp_par structure
#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "diff_evoln_fit.h"
#include "anyoption.h"   // Kishan Thomas' class for command-line option parsing
#include "config_file_parser.h"
//#include "statistics.h"
#include "print_results.h"


/* ---------------- Definitions & Constants ----------------------------- */
#define MAX_N_DATA_VALS   1000000   /* max # data values we'll handle (1.0e6) */

#define NO_FITTING          0 
#define MPFIT_SOLVER        1
#define DIFF_EVOLN_SOLVER   2

#define MAX_DE_GENERATIONS	800

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

#define DEFAULT_CONFIG_FILE   "sample_imfit1d_config.dat"
#define DEFAULT_MODEL_OUTPUT_FILE   "model_profile_save.dat"
#define DEFAULT_OUTPUT_PARAMETER_FILE   "bestfit_parameters_profilefit.dat"

#define VERSION_STRING      " v0.6"



typedef struct {
  std::string  configFileName;
  std::string  dataFileName;
  std::string  psfFileName;
  std::string  modelOutputFileName;
  bool  psfPresent;
  bool  noDataFile;
  bool  noConfigFile;
  bool  dataAreMagnitudes;
  int  startDataRow;
  int  endDataRow;
  bool  noErrors;
  bool  subsamplingFlag;
  bool  saveBestProfile;
  bool  saveBestFitParams;
  std::string  outputParameterFileName;
  int  solver;
} commandOptions;



/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
int myfunc( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */




int myfunc( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *theModel )
{

  theModel->ComputeDeviates(deviates, params);
  return 0;
}


/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  int  nDataVals, nStoredDataVals, nSavedRows;
  int  nPixels_psf;
  int  startDataRow, endDataRow;
  int  nParamsTot, nFreeParams;
  double  *xVals, *yVals, *yWeights;
  double  *xVals_psf, *yVals_psf;
  int  weightMode;
  FILE  *outputFile_ptr;
  ModelObject  *theModel;
//  FunctionObject  *thisFunctionObj;  
  double  *paramsVect;
  double  *paramErrs;
  vector<mp_par>  paramLimits;
  vector<int>  functionSetIndices;
  bool  paramLimitsExist = false;
  bool  parameterInfo_allocated = false;
  mp_result  mpfitResult;
  int  status;
  vector<string>  functionList;
  vector<double>  parameterList;
  mp_par  *parameterInfo;
  mp_par  *mpfitParameterConstraints;
  commandOptions  options;
  configOptions  userConfigOptions;
  
  
  /* PROCESS COMMAND-LINE: */
  /* First, set up the options structure: */
  options.configFileName = DEFAULT_CONFIG_FILE;
  options.dataFileName = "";
  options.modelOutputFileName = DEFAULT_MODEL_OUTPUT_FILE;
  options.noDataFile = true;
  options.noConfigFile = true;
  options.psfPresent = false;
  options.dataAreMagnitudes = true;   // default: assumes we usually fit mu(R) profiles
  options.startDataRow = 0;
  options.endDataRow = -1;   // default value indicating "last row in data file"
  options.noErrors = true;
  options.solver = MPFIT_SOLVER;
  options.subsamplingFlag = false;
  options.saveBestProfile = false;
  options.saveBestFitParams = true;
  options.outputParameterFileName = DEFAULT_OUTPUT_PARAMETER_FILE;

  ProcessInput(argc, argv, &options);

  if (options.noDataFile) {
    printf("*** WARNING: No data to fit!\n\n");
    return -1;
  }

  /* Read configuration file */
  if (! FileExists(options.configFileName.c_str())) {
    printf("\n*** WARNING: Unable to find or open configuration file \"%s\"!\n\n", 
           options.configFileName.c_str());
    return -1;
  }
  status = ReadConfigFile(options.configFileName, false, functionList, parameterList, paramLimits, 
  								functionSetIndices, paramLimitsExist, userConfigOptions);
  if (status < 0) {
    printf("\n*** WARNING: Problem in processing config file!\n\n");
    return -1;
  }


  /* GET THE DATA: */
  nDataVals = CountDataLines(options.dataFileName);
  if ((nDataVals < 1) || (nDataVals > MAX_N_DATA_VALS)) {
	/* file has no data *or* too much data (or an integer overflow occured 
	   in CountDataLines) */
  	printf("Something wrong: input file %s has too few or too many data points\n", 
		   options.dataFileName.c_str());
	printf("(nDataVals = %d)\n", nDataVals);
	exit(1);
  }
  printf("Data file \"%s\": %d data points\n", options.dataFileName.c_str(), nDataVals);
  /* Set default end data row (if not specified) and check for reasonable values: */
  startDataRow = options.startDataRow;
  endDataRow = options.endDataRow;
  if (endDataRow == -1)
    endDataRow = nDataVals - 1;
  if ( (startDataRow < 0) || (startDataRow >= nDataVals) ) {
    printf("Starting data row (\"--x1\") must be >= 1 and <= number of rows in data file (%d)!\n",
    				nDataVals);
    exit(-1);
  }
  if ( (endDataRow <= startDataRow) || (endDataRow >= nDataVals) ) {
    printf("Ending data row (\"--x2\") must be >= starting data row and <= number of rows in data file (%d)!\n",
    				nDataVals);
    exit(-1);
  }
  
  /* Allocate data vectors: */
  nStoredDataVals = endDataRow - startDataRow + 1;
  xVals = (double *)calloc( (size_t)nStoredDataVals, sizeof(double) );
  yVals = (double *)calloc( (size_t)nStoredDataVals, sizeof(double) );
  if ( (xVals == NULL) || (yVals == NULL) ) {
    fprintf(stderr, "\nFailure to allocate memory for input data!\n");
    exit(-1);
  }
  if (options.noErrors == 1)
    yWeights = NULL;
  else {
    yWeights = (double *)calloc( (size_t)nStoredDataVals, sizeof(double) );
  }
  
  /* Read in data */
  nSavedRows = ReadDataFile(options.dataFileName, startDataRow, endDataRow, 
                             xVals, yVals, yWeights);
  if (nSavedRows > nStoredDataVals) {
    fprintf(stderr, "\nMore data rows saved (%d) than we allocated space for (%d)!\n",
            nSavedRows, nStoredDataVals);
    exit(-1);
  }
  
  if (options.noErrors) {
    // OK, we previously had yWeights = NULL to tell ReadDataFile() to skip
    // the third column (if any); now we need to have a yWeights vector with
    // all weights = 1.0
    yWeights = (double *)calloc( (size_t)nStoredDataVals, sizeof(double) );
    for (int i = 0; i < nStoredDataVals; i++)
      yWeights[i] = 1.0;
  } else {
    if (weightMode == WEIGHTS_ARE_SIGMAS) {
      // Construct weights as 1/err
	  printf("Converting errors to weights (w = 1/err)\n");
      for (int i = 0; i < nStoredDataVals; i++)
        yWeights[i] = 1.0/yWeights[i];
    }
  }


  /* Read in PSF profile, if supplied */
  if (options.psfPresent) {
//    printf("Reading PSF profile (\"%s\") ...\n", options.psfFileName.c_str());
    nPixels_psf = CountDataLines(options.psfFileName);
    if ((nPixels_psf < 1) || (nPixels_psf > MAX_N_DATA_VALS)) {
	  /* file has no data *or* too much data (or an integer overflow occured 
	     in CountDataLines) */
  	  printf("Something wrong: input PSF file %s has too few or too many data points\n", 
		     options.psfFileName.c_str());
  	printf("(nPixels_psf = %d)\n", nPixels_psf);
	  exit(1);
    }
    printf("PSF file \"%s\": %d data points\n", options.psfFileName.c_str(), nPixels_psf);
    /* Set default end data row (if not specified) and check for reasonable values: */
    startDataRow = 0;
    endDataRow = nPixels_psf - 1;

    xVals_psf = (double *)calloc( (size_t)nPixels_psf, sizeof(double) );
    yVals_psf = (double *)calloc( (size_t)nPixels_psf, sizeof(double) );
    if ( (xVals_psf == NULL) || (yVals_psf == NULL) ) {
      fprintf(stderr, "\nFailure to allocate memory for PSF data!\n");
      exit(-1);
    }

    nSavedRows = ReadDataFile(options.psfFileName, startDataRow, endDataRow, 
                               xVals_psf, yVals_psf, NULL);
    if (nSavedRows > nStoredDataVals) {
      fprintf(stderr, "\nMore PSF rows saved (%d) than we allocated space for (%d)!\n",
              nSavedRows, nPixels_psf);
      exit(-1);
    }
  }



  /* Set up the model object */
  theModel = new ModelObject1d();
  
  /* Add functions to the model object */
  printf("Adding functions to model object...\n");
  status = AddFunctions1d(theModel, functionList, functionSetIndices);
  if (status < 0) {
  	printf("*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know how many total parameters
  // there will be
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("\t%d total parameters\n", nParamsTot);
  paramsVect = (double *) malloc(nParamsTot * sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  
  /* Add image data and errors to the model object */
  // "true" = input yVals data are magnitudes, not intensities
  theModel->AddDataVectors(nStoredDataVals, xVals, yVals, options.dataAreMagnitudes);
  theModel->AddErrorVector1D(nStoredDataVals, yWeights, WEIGHTS_ARE_SIGMAS);
  theModel->PrintDescription();
  // Add PSF vector, if present, and thereby enable convolution
  if (options.psfPresent)
    theModel->AddPSFVector1D(nPixels_psf, xVals_psf, yVals_psf);


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
  }
  if ((options.solver == MPFIT_SOLVER) && (! paramLimitsExist)) {
    // If parameters are unconstrained, then mpfit() expects a NULL mp_par array
    printf("No parameter constraints!\n");
    mpfitParameterConstraints = NULL;
  } else {
    mpfitParameterConstraints = parameterInfo;
  }

  
  switch (options.solver) {
    case MPFIT_SOLVER:
      bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero the results structure */
      mpfitResult.xerror = paramErrs;
      printf("\nCalling mpfit ...\n");
      status = mpfit(myfunc, nStoredDataVals, nParamsTot, paramsVect, mpfitParameterConstraints, 0, 
  	  								theModel, &mpfitResult);
  
      printf("\n");
      PrintResults(paramsVect, 0, &mpfitResult, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
      break;
    case DIFF_EVOLN_SOLVER:
      printf("\nCalling DiffEvolnFit ..\n");
      status = DiffEvolnFit(nParamsTot, paramsVect, parameterInfo, theModel, MAX_DE_GENERATIONS);
      printf("\n");
      PrintResults(paramsVect, 0, 0, theModel, nFreeParams, parameterInfo, status);
      printf("\n");
      break;
    case NO_FITTING:
      printf("\nNO FITTING BEING DONE!\n");
      theModel->SetupChisquaredCalcs();
      options.saveBestProfile = true;
      break;
  } // end switch

  if (options.saveBestFitParams)
    SaveParameters(paramsVect, theModel, parameterInfo, options.outputParameterFileName,
    								argc, argv);

  if (options.saveBestProfile) {
    double chisqr = theModel->ChiSquared(paramsVect);
    double *modelProfile = (double *) calloc((size_t)nStoredDataVals, sizeof(double));
    int nPts = theModel->GetModelVector(modelProfile);
    if (nPts == nStoredDataVals) {
      printf("Saving model profile to %s...\n", options.modelOutputFileName.c_str());
      outputFile_ptr = fopen(options.modelOutputFileName.c_str(), "w");
      for (int i = 0; i < nPts; i++) {
        fprintf(outputFile_ptr, "\t%f\t%f\n", xVals[i], modelProfile[i]);
      }
      fclose(outputFile_ptr);
      printf("Done.\n");
    }
    else {
      printf("WARNING -- MISMATCH BETWEEN nStoredDataVals (main) and nDataVals (ModelObject1d)!\n");
      printf("NO PROFILE SAVED!\n");
    }
    free(modelProfile);
  }

  // Free up memory
  free(xVals);
  free(yVals);
  free(yWeights);
  if (options.psfPresent) {
    free(xVals_psf);
    free(yVals_psf);
  }
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
  //opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* SET THE USAGE/HELP   */
  opt->addUsage("Usage: ");
  opt->addUsage("   profilefit [options] datafile configfile");
  opt->addUsage(" -h  --help                   Prints this help");
  opt->addUsage(" --useerrors                  Use errors from data file");
  opt->addUsage(" --intensities                Data y-values are intensities, not magnitudes");
  opt->addUsage(" --psf <psf_file>             PSF image");
  opt->addUsage(" --de                         Solve using differential evolution");
  opt->addUsage(" --no-fitting                 Don't do fitting (just save input model)");
  opt->addUsage(" --x1 <int>                   start data value");
  opt->addUsage(" --x2 <int>                   end data value");
  opt->addUsage(" --save-params <output-file>  Save best-fit parameters in config-file format");
  opt->addUsage("");


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag("help", 'h');
  opt->setFlag("useerrors");
  opt->setFlag("intensities");
  opt->setOption("psf");      /* an option (takes an argument), supporting only long form */
  opt->setFlag("de");
  opt->setFlag("no-fitting");
  opt->setOption("x1");      /* an option (takes an argument), supporting only long form */
  opt->setOption("x2");        /* an option (takes an argument), supporting only long form */
  opt->setOption("save-params");
  
  /* parse the command line:  */
  opt->processCommandArgs( argc, argv );


  /* Process the results: actual arguments, if any: */
  if (opt->getArgc() > 0) {
    theOptions->dataFileName = opt->getArgv(0);
    theOptions->noDataFile = false;
    printf("\tdata file = %s\n", theOptions->dataFileName.c_str());
  }
  if (opt->getArgc() > 1) {
    theOptions->configFileName = opt->getArgv(1);
    theOptions->noConfigFile = false;
    printf("\tconfig file = %s\n", theOptions->configFileName.c_str());
  }

  /* Process the results: options */
  if ( opt->getFlag("help") || opt->getFlag('h') || (! opt->hasOptions()) ) {
    opt->printUsage();
    delete opt;
    exit(1);
  }
  if (opt->getFlag("useerrors")) {
  	printf("\t USE ERRORS SELECTED!\n");
  	theOptions->noErrors = false;
  }
  if (opt->getFlag("intensities")) {
  	printf("\t Data values are assumed to be intensities (instead of magnitudes)!\n");
  	theOptions->dataAreMagnitudes = false;
  }
  if (opt->getValue("psf") != NULL) {
    theOptions->psfFileName = opt->getValue("psf");
    theOptions->psfPresent = true;
    printf("\tPSF profile = %s\n", theOptions->psfFileName.c_str());
  }
  if (opt->getFlag("de")) {
  	printf("\t Differential Evolution selected!\n");
  	theOptions->solver = DIFF_EVOLN_SOLVER;
  }
  if (opt->getFlag("no-fitting")) {
  	printf("\t No fitting will be done!\n");
  	theOptions->solver = NO_FITTING;
  }
  if (opt->getValue("x1") != NULL) {
    if (NotANumber(opt->getValue("x1"), 0, kPosInt)) {
      printf("*** WARNING: start data row should be a positive integer!\n");
      delete opt;
      exit(1);
    }
    theOptions->startDataRow = atol(opt->getValue("x1"));
    printf("\tstart data row = %d\n", theOptions->startDataRow);
  }
  if (opt->getValue("x2") != NULL) {
    if (NotANumber(opt->getValue("x2"), 0, kPosInt)) {
      printf("*** WARNING: end data row should be a positive integer!\n");
      delete opt;
      exit(1);
    }
    theOptions->endDataRow = atol(opt->getValue("x2"));
    printf("\tend data row = %d\n", theOptions->endDataRow);
  }
  if (opt->getValue("save-params") != NULL) {
    theOptions->outputParameterFileName = opt->getValue("save-params");
    theOptions->saveBestFitParams = true;
    printf("\toutput best-fit parameter file = %s\n", theOptions->outputParameterFileName.c_str());
  }


  delete opt;

}






/* END OF FILE: imfit1d_main.cpp  ---------------------------------------- */
