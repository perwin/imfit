/* FILE: profilefit_main.cpp  -------------------------------------------- */
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
#include "bootstrap_errors_1d.h"


// Solvers (optimization algorithms)
#include "dispatch_solver.h"
#include "solver_results.h"
#include "levmar_fit.h"
#include "diff_evoln_fit.h"
#ifndef NO_NLOPT
#include "nmsimplex_fit.h"
#include "nlopt_fit.h"
#endif

#include "commandline_parser.h"
#include "config_file_parser.h"
#include "print_results.h"


/* ---------------- Definitions & Constants ----------------------------- */
#define MAX_N_DATA_VALS   1000000   /* max # data values we'll handle (1.0e6) */

#define MONTE_CARLO_ITER   100
#define BOOTSTRAP_ITER     1000

#define CMDLINE_ERROR1 "Usage: -p must be followed by a string containing initial parameter values for the model"
#define CMDLINE_ERROR2 "Usage: -l must be followed by a filename for a file containing parameter limits"
#define CMDLINE_ERROR3 "For differential-evolution solver, a file containing parameter limits *must* be supplied (\"-l\" option)"
#define CMDLINE_ERROR4 "Usage: --mcoffset must be followed by a positive number"
#define CMDLINE_ERROR5 "Usage: --magzp must be followed by a positive number"
#define CMDLINE_ERROR6 "Usage: --mciter must be followed by a positive integer"
#define CMDLINE_ERROR7 "Usage: --pf must be followed by a positive real number"

#define DEFAULT_CONFIG_FILE   "sample_imfit1d_config.dat"
#define DEFAULT_MODEL_OUTPUT_FILE   "model_profile_save.dat"
#define DEFAULT_1D_OUTPUT_PARAMETER_FILE   "bestfit_parameters_profilefit.dat"

#define VERSION_STRING      "v1.5"

#ifndef NO_SIGNALS
volatile sig_atomic_t  stopSignal_flag = 0;
#endif



typedef struct {
  std::string  configFileName;
  std::string  dataFileName;
  std::string  psfFileName;
  std::string  modelOutputFileName;
  bool  psfPresent;
  bool  noDataFile;
  bool  noConfigFile;
  bool  dataAreMagnitudes;
  double  zeroPoint;
  int  startDataRow;
  int  endDataRow;
  bool  noErrors;
  bool  noMask;
  int  maskFormat;
  double  ftol;
  bool  ftolSet;
  bool  subsamplingFlag;
  bool  saveBestProfile;
  bool  saveBestFitParams;
  std::string  outputParameterFileName;
  bool printChiSquaredOnly;
  int  solver;
  std::string  nloptSolverName;
  bool  doBootstrap;
  int  bootstrapIterations;
  int  verbose;
  unsigned long  rngSeed;
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
  int  nDegFreedom;
  double  *xVals, *yVals, *yWeights, *maskVals;
  double  *xVals_psf, *yVals_psf;
  int  weightMode = WEIGHTS_ARE_SIGMAS;
  FILE  *outputFile_ptr;
  ModelObject  *theModel;
  double  *paramsVect;
  double  *paramErrs;
  vector<mp_par>  paramLimits;
  vector<int>  FunctionBlockIndices;
  bool  maskAllocated = false;
  bool  paramLimitsExist = false;
  vector<mp_par>  parameterInfo;
  int  status, fitStatus;
  vector<string>  functionList;
  vector<string>  functionLabelList;
  vector<double>  parameterList;
  commandOptions  options;
  configOptions  userConfigOptions;
  SolverResults  resultsFromSolver;
  string  nloptSolverName;
  vector<string> programHeader;
  string  progNameVersion = "profilefit ";
  
  
  /* PROCESS COMMAND-LINE: */
  /* First, set up the options structure: */
  options.configFileName = DEFAULT_CONFIG_FILE;
  options.dataFileName = "";
  options.modelOutputFileName = DEFAULT_MODEL_OUTPUT_FILE;
  options.noDataFile = true;
  options.noConfigFile = true;
  options.psfPresent = false;
  options.dataAreMagnitudes = true;   // default: assumes we usually fit mu(R) profiles
  options.zeroPoint = 0.0;
  options.startDataRow = 0;
  options.endDataRow = -1;   // default value indicating "last row in data file"
  options.noErrors = true;
  options.noMask = true;
  options.maskFormat = MASK_ZERO_IS_GOOD;
  options.ftol = DEFAULT_FTOL;
  options.ftolSet = false;
  options.printChiSquaredOnly = false;
  options.solver = MPFIT_SOLVER;
  options.nloptSolverName = "NM";   // default value = Nelder-Mead Simplex
  options.doBootstrap = false;
  options.bootstrapIterations = 0;
  options.subsamplingFlag = false;
  options.saveBestProfile = false;
  options.saveBestFitParams = true;
  options.outputParameterFileName = DEFAULT_1D_OUTPUT_PARAMETER_FILE;
  options.verbose = 1;
  options.nloptSolverName = "";
  options.rngSeed = 0;

  progNameVersion += VERSION_STRING;
  MakeOutputHeader(&programHeader, progNameVersion, argc, argv);

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
  status = ReadConfigFile(options.configFileName, false, functionList, functionLabelList,
  					parameterList, paramLimits, FunctionBlockIndices, paramLimitsExist, 
  					userConfigOptions);
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
  if (options.noErrors)
    yWeights = NULL;
  else
    yWeights = (double *)calloc( (size_t)nStoredDataVals, sizeof(double) );
  if (options.noMask)
    maskVals = NULL;
  else {
    maskVals = (double *)calloc( (size_t)nStoredDataVals, sizeof(double) );
    maskAllocated = true;
  }
  
  /* Read in data */
  nSavedRows = ReadDataFile(options.dataFileName, startDataRow, endDataRow, 
                             xVals, yVals, yWeights, maskVals);
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
  } 


  /* Read in PSF profile, if supplied */
  if (options.psfPresent) {
    nPixels_psf = CountDataLines(options.psfFileName);
    if ((nPixels_psf < 1) || (nPixels_psf > MAX_N_DATA_VALS)) {
      /* file has no data *or* too much data (or an integer overflow occured 
         in CountDataLines) */
      printf("Something wrong: input PSF file %s has too few or too many data points\n", 
             options.psfFileName.c_str());
      printf("(nPixels_psf (# of PSF points) = %d)\n", nPixels_psf);
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
                               xVals_psf, yVals_psf, NULL, NULL);
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
  status = AddFunctions1d(theModel, functionList, FunctionBlockIndices);
  if (status < 0) {
    printf("*** WARNING: Failure in AddFunctions!\n\n");
    exit(-1);
  }
  theModel->SetZeroPoint(options.zeroPoint);
  
  // Set up parameter vector(s), now that we know how many total parameters
  // there will be
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("\t%d total parameters\n", nParamsTot);
  paramsVect = (double *) malloc(nParamsTot * sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  
  /* Add image data, errors, and mask to the model object */
  // "true" = input yVals data are magnitudes, not intensities
  theModel->AddDataVectors(nStoredDataVals, xVals, yVals, options.dataAreMagnitudes);
  theModel->AddErrorVector1D(nStoredDataVals, yWeights, WEIGHTS_ARE_SIGMAS);
  if (maskAllocated) {
    status = theModel->AddMaskVector1D(nStoredDataVals, maskVals, options.maskFormat);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddMaskVector1D!\n\n");
  	  exit(-1);
    }
  }
  // Add PSF vector, if present, and thereby enable convolution
  if (options.psfPresent) {
    status = theModel->AddPSFVector1D(nPixels_psf, xVals_psf, yVals_psf);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddPSFVector1D!\n\n");
  	  exit(-1);
    }
  }
  
  theModel->FinalSetupForFitting();   // calls ApplyMask(), VetDataVector()
  theModel->PrintDescription();


  // Parameter limits and other info:
  printf("Setting up parameter information vector ...\n");
  mp_par newParamLimit;
  for (int i = 0; i < nParamsTot; i++) {
    memset(&newParamLimit, 0, sizeof(mp_par));
    newParamLimit.fixed = paramLimits[i].fixed;
    if (paramLimits[i].fixed == 1) {
      printf("Fixed parameter detected (i = %d)\n", i);
      nFreeParams--;
    }
    newParamLimit.limited[0] = paramLimits[i].limited[0];
    newParamLimit.limited[1] = paramLimits[i].limited[1];
    newParamLimit.limits[0] = paramLimits[i].limits[0];
    newParamLimit.limits[1] = paramLimits[i].limits[1];
    parameterInfo.push_back(newParamLimit);
  }
  nDegFreedom = theModel->GetNValidPixels() - nFreeParams;
  printf("%d free parameters (%d degrees of freedom)\n", nFreeParams, nDegFreedom);

  // tell ModelObject about parameterInfo (mainly useful for printing-related methods)
  theModel->AddParameterInfo(parameterInfo);


  // OK, now we either print chi^2 value for the input parameters and quit, or
  // else call one of the solvers!
  if (options.printChiSquaredOnly) {
    printf("\n");
    fitStatus = 1;
    PrintFitStatistic(paramsVect, theModel, nFreeParams);
    printf("\n");
    // turn off saveing of parameter file
    options.saveBestFitParams = false;
  }
  else {
    // DO THE FIT!
    fitStatus = DispatchToSolver(options.solver, nParamsTot, nFreeParams, nStoredDataVals, 
    							paramsVect, parameterInfo, theModel, options.ftol, paramLimitsExist, 
    							options.verbose, &resultsFromSolver, options.nloptSolverName,
    							options.rngSeed);
    							
    PrintResults(paramsVect, theModel, nFreeParams, fitStatus, resultsFromSolver);
  }



  // OPTIONAL HANDLING OF BOOTSTRAP RESAMPLING HERE
  if ((options.doBootstrap) && (options.bootstrapIterations > 0)) {
    printf("\nNow doing bootstrap resampling (%d iterations) to estimate errors...\n",
           options.bootstrapIterations);
    printf("[NOT YET PROPERLY IMPLEMENTED!]\n");
    BootstrapErrors(paramsVect, parameterInfo, paramLimitsExist, theModel,
                    options.ftol, options.bootstrapIterations, nFreeParams);
  }
  
  
  
  if (options.saveBestFitParams) {
    printf("Saving best-fit parameters in file \"%s\"\n", options.outputParameterFileName.c_str());
    string  progNameVer = "profilefit ";
    progNameVer += VERSION_STRING;
    SaveParameters(paramsVect, theModel, options.outputParameterFileName,
                    programHeader, nFreeParams, options.solver, fitStatus, resultsFromSolver);
  }


  if (options.saveBestProfile) {
    theModel->CreateModelImage(paramsVect);
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
  if (maskAllocated)
    free(maskVals);
  if (options.psfPresent) {
    free(xVals_psf);
    free(yVals_psf);
  }
  free(paramsVect);
  free(paramErrs);
  delete theModel;
  
  printf("All done!\n\n");
  return 0;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   profilefit [options] datafile configfile");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine(" -v  --version                Prints version number");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --list-parameters        Prints list of parameter names for each available function");
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" --useerrors                  Use errors from data file (3rd column)");
  optParser->AddUsageLine(" --usemask                    Use mask from data file (4th column)");
  optParser->AddUsageLine(" --intensities                Data y-values are intensities, not magnitudes");
  optParser->AddUsageLine(" --psf <psf_file>             PSF profile (centered on middle row, y-values = intensities)");
  optParser->AddUsageLine("");
#ifndef NO_NLOPT
  optParser->AddUsageLine(" --nm                         Use Nelder-Mead simplex solver instead of L-M");
  optParser->AddUsageLine(" --nlopt <name>               Select misc. NLopt solver");
#endif
  optParser->AddUsageLine(" --de                         Solve using differential evolution");
  optParser->AddUsageLine("");
  optParser->AddUsageLine(" --ftol                       Fractional tolerance in chi^2 for convergence [default = 1.0e-8]");
  optParser->AddUsageLine(" --chisquare-only             Print chi^2 of input model and quit");
  optParser->AddUsageLine(" --no-fitting                 Don't do fitting (just save input model)");
  optParser->AddUsageLine(" --x1 <int>                   start data value (1 = first, 2 = second, etc.)");
  optParser->AddUsageLine(" --x2 <int>                   end data value");
  optParser->AddUsageLine(" --zp <float>                 magnitude zero point of the data");
  optParser->AddUsageLine(" --bootstrap <int>            do this many iterations of bootstrap resampling to estimate errors");
  optParser->AddUsageLine(" --save-params <output-file>  Save best-fit parameters in config-file format");
  optParser->AddUsageLine(" --save-best-fit <output-file>  Save best-fit profile");
  optParser->AddUsageLine(" --seed <int>                 RNG seed (for testing purposes)");
  optParser->AddUsageLine(" --quiet                      Turn off printing of mpfit interation updates");
  optParser->AddUsageLine("");


  /* by default all options are checked on the command line and from option/resource file */
  optParser->AddFlag("help", "h");
  optParser->AddFlag("version", "v");
  optParser->AddFlag("list-functions");
  optParser->AddFlag("list-parameters");
  optParser->AddFlag("useerrors");
  optParser->AddFlag("usemask");
  optParser->AddFlag("intensities");
  optParser->AddOption("psf");      /* an option (takes an argument), supporting only long form */
#ifndef NO_NLOPT
  optParser->AddFlag("nm");
  optParser->AddOption("nlopt");
#endif
  optParser->AddFlag("de");
  optParser->AddOption("ftol");
  optParser->AddFlag("chisquare-only");
  optParser->AddFlag("no-fitting");
  optParser->AddOption("x1");      /* an option (takes an argument), supporting only long form */
  optParser->AddOption("x2");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("zp");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("bootstrap");        /* an option (takes an argument), supporting only long form */
  optParser->AddOption("save-params");
  optParser->AddOption("save-best-fit");
  optParser->AddFlag("verbose");
  optParser->AddFlag("quiet");
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
    theOptions->dataFileName = optParser->GetArgument(0);
    theOptions->noDataFile = false;
    printf("\tdata file = %s\n", theOptions->dataFileName.c_str());
  }
  if (optParser->nArguments() > 1) {
    theOptions->configFileName = optParser->GetArgument(1);
    theOptions->noConfigFile = false;
    printf("\tconfig file = %s\n", theOptions->configFileName.c_str());
  }

  /* Process the results: options */
  if ( optParser->FlagSet("help")  || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("version") ) {
    printf("profilefit version %s\n\n", VERSION_STRING);
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

  if (optParser->FlagSet("useerrors")) {
    printf("\t USE ERRORS SELECTED!\n");
    theOptions->noErrors = false;
  }
  if (optParser->FlagSet("usemask")) {
    printf("\t USE MASK SELECTED!\n");
    theOptions->noMask = false;
  }
  if (optParser->FlagSet("intensities")) {
    printf("\t Data values are assumed to be intensities (instead of magnitudes)!\n");
    theOptions->dataAreMagnitudes = false;
  }
  if (optParser->OptionSet("psf")) {
    theOptions->psfFileName = optParser->GetTargetString("psf");
    theOptions->psfPresent = true;
    printf("\tPSF profile = %s\n", theOptions->psfFileName.c_str());
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
    printf("\t Differential Evolution selected!\n");
    theOptions->solver = DIFF_EVOLN_SOLVER;
  }
  if (optParser->FlagSet("chisquare-only")) {
    printf("\t* No fitting will be done!\n");
    theOptions->printChiSquaredOnly = true;
  }
  if (optParser->FlagSet("no-fitting")) {
    printf("\t No fitting will be done!\n");
    theOptions->solver = NO_FITTING;
  }
  if (optParser->FlagSet("quiet")) {
    theOptions->verbose = 0;
  }
  if (optParser->FlagSet("verbose")) {
    theOptions->verbose = 2;
  }
  if (optParser->OptionSet("x1")) {
    if (NotANumber(optParser->GetTargetString("x1").c_str(), 0, kPosInt)) {
      printf("*** WARNING: start data row should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->startDataRow = atol(optParser->GetTargetString("x1").c_str());
    printf("\tstart data row = %d\n", theOptions->startDataRow);
  }
  if (optParser->OptionSet("x2")) {
    if (NotANumber(optParser->GetTargetString("x2").c_str(), 0, kPosInt)) {
      printf("*** WARNING: end data row should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->endDataRow = atol(optParser->GetTargetString("x2").c_str());
    printf("\tend data row = %d\n", theOptions->endDataRow);
  }
  if (optParser->OptionSet("zp")) {
    if (NotANumber(optParser->GetTargetString("zp").c_str(), 0, kAnyReal)) {
      printf("*** WARNING: zero point should be a real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->zeroPoint = atof(optParser->GetTargetString("zp").c_str());
    printf("\tmagnitude zero point = %f\n", theOptions->zeroPoint);
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
  if (optParser->OptionSet("bootstrap")) {
    if (NotANumber(optParser->GetTargetString("bootstrap").c_str(), 0, kPosInt)) {
      printf("*** WARNING: number of bootstrap iterations should be a positive integer!\n");
      delete optParser;
      exit(1);
    }
    theOptions->doBootstrap = true;
    theOptions->bootstrapIterations = atol(optParser->GetTargetString("bootstrap").c_str());
    printf("\tnumber of bootstrap iterations = %d\n", theOptions->bootstrapIterations);
  }
  if (optParser->OptionSet("save-params")) {
    theOptions->outputParameterFileName = optParser->GetTargetString("save-params");
    theOptions->saveBestFitParams = true;
    printf("\toutput best-fit parameter file = %s\n", theOptions->outputParameterFileName.c_str());
  }
  if (optParser->OptionSet("save-best-fit")) {
    theOptions->modelOutputFileName = optParser->GetTargetString("save-best-fit");
    theOptions->saveBestProfile = true;
    printf("\toutput best-fit profile to file = %s\n", theOptions->modelOutputFileName.c_str());
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




/* END OF FILE: profilefit_main.cpp  ------------------------------------- */
