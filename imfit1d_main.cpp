/* FILE: imfit1d_main.cpp  ----------------------------------------------- */
/*
 * *** TRAIL VERSION INCORPORATING FUNCTION OBJECTS!!! ***
 * 
 * The proper translations are:
 * NAXIS1 = naxes[0] = nColumns = sizeX;
 * NAXIS2 = naxes[1] = nRows = sizeY.
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
#include "add_functions.h"
#include "function_object.h"
#include "func1d_exp.h"
#include "mpfit_cpp.h"   // lightly modified mpfit from Craig Markwardt
#include "anyoption.h"   // Kishan Thomas' class for command-line option parsing
#include "config_file_parser.h"


/* ---------------- Definitions & Constants ----------------------------- */
#define MAX_N_DATA_VALS   1000000   /* max # data values we'll handle (1.0e6) */

#define LMDIF_SOLVER        1
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
#define DEFAULT_CONFIG_FILE   "sample_imfit1d_config.dat"

#define VERSION_STRING      " v0.05"



typedef struct {
  std::string  configFileName;
  std::string  dataFileName;
  bool  noDataFile;
  bool  noConfigFile;
  int  startDataRow;
  int  endDataRow;
  bool  noErrors;
} commandOptions;



/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
void PrintParam( string& paramName, double paramValue, double paramErr );
void PrintResult( double *x, double *xact, mp_result *result,
					ModelObject *model );
int myfunc( int nDataVals, int nParams, double *params, double *deviates,
           double **derivatives, ModelObject *aModel );


/* ------------------------ Global Variables --------------------------- */

/* ------------------------ Module Variables --------------------------- */

static char  programName[MAX_FILENAME_LENGTH];
static char  errorString[MAXLINE];



void PrintParam( string& paramName, double paramValue, double paramErr )
{
  printf("  %10s = %f +/- %f\n", paramName.c_str(), paramValue, paramErr);
}


/* Simple routine to print the fit results [taken & modified from Craig Markwardt's
   testmpfit.c] 
*/
void PrintResult( double *x, double *xact, mp_result *result, ModelObject *model )
{
  int i;

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
  if (xact) {
    for (i=0; i<result->npar; i++) {
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
  int  nDegFreedom;
  double  *xVals, *yVals, *yWeights;
  int  weightMode;
  ModelObject  *theModel;
//  FunctionObject  *thisFunctionObj;  
  double  *paramsVect;
  double  *paramErrs;
  vector<mp_par>  paramLimits;
  bool  paramLimitsExist = false;
  mp_result  mpfitResult;
//  mp_par  mpParameterBounds[10];
  int  status;
  vector<string>  functionList;
  vector<double>  parameterList;
  mp_par  *mpParamLimits;
  commandOptions  options;
  
  
  /* PROCESS COMMAND-LINE: */
  options.configFileName = DEFAULT_CONFIG_FILE;
  options.dataFileName = "";
  options.noDataFile = true;
  options.noConfigFile = true;
  options.startDataRow = 0;
  options.endDataRow = -1;   // default value indicating "last row in data file"
  options.noErrors = true;

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
  // [] READ IN CONFIGURATION FILE ...
  ReadConfigFile(options.configFileName, functionList, parameterList, paramLimits, 
  								paramLimitsExist);
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
  status = AddFunctions(theModel, functionList);
  if (status < 0) {
  	printf("*** WARNING: Failure in AddFunctions!\n\n");
  	exit(-1);
 }
  
  // Set up parameter vector(s), now that we know how many total parameters
  // there will be
  nParamsTot = nFreeParams = theModel->GetNParams();
  printf("%d total parameters\n", nParamsTot);
  // [] create & add to parameter vector here?
  paramsVect = (double *) malloc(nParamsTot * sizeof(double));
  for (int i = 0; i < nParamsTot; i++)
    paramsVect[i] = parameterList[i];
  paramErrs = (double *) malloc(nParamsTot * sizeof(double));
  
  /* Add image data and errors to the model object */
  // "true" = input yVals data are magnitudes, not intensities
  theModel->AddDataVectors(nStoredDataVals, xVals, yVals, true);
  theModel->AddErrorVector1D(nStoredDataVals, yWeights, WEIGHTS_ARE_SIGMAS);
  theModel->PrintDescription();
//  theModel->AddErrorVector(nPixels_tot, nColumns, nRows, allErrorPixels,
//                           WEIGHTS_ARE_SIGMAS);

  // TESTING (remove later)
  // Create a vector to hold deviates, then ask the model object to
  // compute the deviates, given a trial set of parameters
  //diffVector = (double *) malloc(nPixels_tot * sizeof(double));
  //theModel->ComputeDeviates(diffVector, paramsVect);
  // Ask the model object to print out the model image it just computed
  //if (options.printImages)
  //  theModel->PrintModelImage();


  // Some trial mpfit experiments:
  bzero(&mpfitResult, sizeof(mpfitResult));       /* Zero results structure */
  mpfitResult.xerror = paramErrs;
  
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
  printf("\nStarting mpfit w/ nStoredDataVals = %d, nParamsTot = %d ...\n",
  				nStoredDataVals, nParamsTot);
  status = mpfit(myfunc, nStoredDataVals, nParamsTot, paramsVect, mpParamLimits, 0, 
  								theModel, &mpfitResult);
  
  printf("*** mpfit status = %d\n", status);
  PrintResult(paramsVect, 0, &mpfitResult, theModel);

  nDegFreedom = theModel->GetNValidPixels() - nFreeParams;
  printf("\nTrue degrees of freedom = %d\n", nDegFreedom);

  // Free up memory
  free(xVals);
  free(yVals);
  free(yWeights);
  free(paramsVect);
  free(paramErrs);
  //free(mpParameterBounds);
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
  opt->addUsage(" -c  --useerrors              Use errors from data file");
  opt->addUsage("     --x1 <int>               start data value");
  opt->addUsage("     --x2 <int>               end data value");
  opt->addUsage("");


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag("help", 'h');
  opt->setFlag("useerrors");
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
