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

#define VERSION_STRING      " v0.5"



typedef struct {
  std::string  configFileName;
  std::string  dataFileName;
  bool  noDataFile;
  bool  noConfigFile;
  int  startDataRow;
  int  endDataRow;
  bool  noErrors;
  bool  subsamplingFlag;
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
  int  startDataRow, endDataRow;
  int  nParamsTot, nFreeParams;
  double  *xVals, *yVals, *yWeights;
  int  weightMode;
  ModelObject  *theModel;
//  FunctionObject  *thisFunctionObj;  
  double  *paramsVect;
  double  *paramErrs;
  vector<mp_par>  paramLimits;
  vector<int>  functionSetIndices;
  bool  paramLimitsExist = false;
  mp_result  mpfitResult;
  int  status;
  vector<string>  functionList;
  vector<double>  parameterList;
  mp_par  *mpParamLimits;
  commandOptions  options;
  configOptions  userConfigOptions;
  
  
  /* PROCESS COMMAND-LINE: */
  /* First, set up the options structure: */
  options.configFileName = DEFAULT_CONFIG_FILE;
  options.dataFileName = "";
  options.noDataFile = true;
  options.noConfigFile = true;
  options.startDataRow = 0;
  options.endDataRow = -1;   // default value indicating "last row in data file"
  options.noErrors = true;
  options.solver = MPFIT_SOLVER;
  options.subsamplingFlag = false;

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
  for (int k = 0; k < parameterList.size(); k++)
  	cout << parameterList[k] << "  ";
  cout << endl;


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
  


  /* Set up the model object */
  theModel = new ModelObject1d();
  
  /* Add functions to the model object */
  status = AddFunctions1d(theModel, functionList, functionSetIndices);
  if (status < 0) {
  	printf("*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know how many total parameters
  // there will be
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  paramsVect = (double *) malloc(nParamsTot * sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  
  /* Add image data and errors to the model object */
  // "true" = input yVals data are magnitudes, not intensities
  theModel->AddDataVectors(nStoredDataVals, xVals, yVals, true);
  theModel->AddErrorVector1D(nStoredDataVals, yWeights, WEIGHTS_ARE_SIGMAS);
  theModel->PrintDescription();


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
  
  cout << "Here is the input paramsVect:" << endl;
  for (int k = 0; k < nParamsTot; k++)
  	cout << paramsVect[k] << "  ";
  cout << endl;
  printf("\nStarting the fit w/ nStoredDataVals = %d, nParamsTot = %d ...\n",
  				nStoredDataVals, nParamsTot);
  
  
  if (options.solver == MPFIT_SOLVER) {
    bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero the results structure */
    mpfitResult.xerror = paramErrs;
    printf("Calling mpfit ...\n");
    status = mpfit(myfunc, nStoredDataVals, nParamsTot, paramsVect, mpParamLimits, 0, 
  									theModel, &mpfitResult);
  
    printf("\n");
    PrintResults(paramsVect, 0, &mpfitResult, theModel, nFreeParams, mpParamLimits, status);
    printf("\n");
  }
  else {
    printf("Calling DiffEvolnFit ..\n");
    status = DiffEvolnFit(nParamsTot, paramsVect, mpParamLimits, theModel, MAX_DE_GENERATIONS);
    printf("\n");
    PrintResults(paramsVect, 0, 0, theModel, nFreeParams, mpParamLimits, status);
    printf("\n");
  }


  // Free up memory
  free(xVals);
  free(yVals);
  free(yWeights);
  free(paramsVect);
  free(paramErrs);
  free(theModel);
  
  return 0;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  AnyOption *opt = new AnyOption();
  //opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* SET THE USAGE/HELP   */
  opt->addUsage("Usage: ");
  opt->addUsage("   imfit1d [options] datafile configfile");
  opt->addUsage(" -h  --help                   Prints this help");
  opt->addUsage(" --useerrors                  Use errors from data file");
  opt->addUsage(" --de                         Solve using differential evolution");
  opt->addUsage("     --x1 <int>               start data value");
  opt->addUsage("     --x2 <int>               end data value");
  opt->addUsage("");


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag("help", 'h');
  opt->setFlag("useerrors");
  opt->setFlag("de");
  opt->setOption("x1");      /* an option (takes an argument), supporting only long form */
  opt->setOption("x2");        /* an option (takes an argument), supporting only long form */

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
  	printf("\t USEERRORS SELECTED!\n");
  	theOptions->noErrors = false;
  }
  if (opt->getFlag("de")) {
  	printf("\t Differential Evolution selected!\n");
  	theOptions->solver = DIFF_EVOLN_SOLVER;
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


  delete opt;

}






/* END OF FILE: imfit1d_main.cpp  ---------------------------------------- */
