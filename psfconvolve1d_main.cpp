// Test code for image convolution

// This is a modified version of psfconvolve_main.cpp, for purposes of developing
// and testing 1-D convolution of profiles.

// NOTE: GCC support for C99 "complex" types is "broken" as of version 4.4
// (and for all earlier versions, including 4.2).

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "fftw3.h"

#include "convolver1d.h"
#include "read_profile_pub.h"
#include "commandline_parser.h"
#include "utilities_pub.h"


/* ---------------- Definitions ---------------------------------------- */
#define MAX_N_DATA_VALS   1000000   /* max # data values we'll handle (1.0e6) */

#define INPUT_PROFILE_FILENAME   "some_profile.dat"
#define PSF_FILENAME    "psf_gaussian1d.dat"
#define DEFAULT_OUTPUT_FILENAME   "convolve1d_out.dat"

#define  FILE_OPENW_ERR_STRING "\n   Couldn't open file \"%s\" for writing!\n\n"


typedef struct {
  std::string  dataFileName;
  std::string  psfProfileName;
  std::string  outputProfileName;
  int  startDataRow;
  int  endDataRow;
  bool  printProfiles;
  bool  outputPaddedProfile;
  int  debugLevel;
} commandOptions;


/* ------------------- Function Prototypes ----------------------------- */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
// void ShiftAndWrapPSF( double *psfImage, int nRows_psf, int nCols_psf,
//                       fftw_complex *destImage, int nRows_dest, int nCols_dest );
void SaveProfile( std::string& fileName, double *xVals, double *yVals, int nDataVals );


/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  int  nDataVals, nStoredDataVals, nSavedRows;
  int  startDataRow, endDataRow;
  int  nPSFVals;
  std::string  imageFilename, psfFilename, outputFilename;
  double  *xVals, *yVals;
  double  *xVals_psf, *yVals_psf;
  commandOptions  options;
  Convolver1D  psfConvolver;


  /* Process command line: */
  options.dataFileName = INPUT_PROFILE_FILENAME;
  options.psfProfileName = PSF_FILENAME;
  options.outputProfileName = DEFAULT_OUTPUT_FILENAME;
  options.startDataRow = 0;
  options.endDataRow = -1;   // default value indicating "last row in data file"
  options.printProfiles = false;
  options.outputPaddedProfile = false;
  options.debugLevel = 0;

  ProcessInput(argc, argv, &options);


  printf("\nReading input profile (\"%s\") ...\n", options.dataFileName.c_str());
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

  /* Read in data */
  nSavedRows = ReadDataFile(options.dataFileName, startDataRow, endDataRow, 
                             xVals, yVals, NULL, NULL);
  if (nSavedRows > nStoredDataVals) {
    fprintf(stderr, "\nMore data rows saved (%d) than we allocated space for (%d)!\n",
            nSavedRows, nStoredDataVals);
    exit(-1);
  }
  
  if (options.debugLevel >= 3) {
    printf("The input data:\n");
    for (int k = 0; k < nStoredDataVals; k++) {
      printf("%f\t%f\n", xVals[k], yVals[k]);
    }
  }


  printf("Reading PSF profile (\"%s\") ...\n", options.psfProfileName.c_str());
  nPSFVals = CountDataLines(options.psfProfileName);
  if ((nPSFVals < 1) || (nPSFVals > MAX_N_DATA_VALS)) {
	/* file has no data *or* too much data (or an integer overflow occured 
	   in CountDataLines) */
  	printf("Something wrong: input file %s has too few or too many data points\n", 
		   options.psfProfileName.c_str());
	printf("(nPSFVals = %d)\n", nPSFVals);
	exit(1);
  }
  printf("PSF file \"%s\": %d data points\n", options.psfProfileName.c_str(), nPSFVals);
  /* Set default end data row (if not specified) and check for reasonable values: */
  startDataRow = 0;
  endDataRow = nPSFVals - 1;

  xVals_psf = (double *)calloc( (size_t)nPSFVals, sizeof(double) );
  yVals_psf = (double *)calloc( (size_t)nPSFVals, sizeof(double) );
  if ( (xVals_psf == NULL) || (yVals_psf == NULL) ) {
    fprintf(stderr, "\nFailure to allocate memory for PSF data!\n");
    exit(-1);
  }

  nSavedRows = ReadDataFile(options.psfProfileName, startDataRow, endDataRow, 
                             xVals_psf, yVals_psf, NULL, NULL);
  if (nSavedRows > nStoredDataVals) {
    fprintf(stderr, "\nMore PSF rows saved (%d) than we allocated space for (%d)!\n",
            nSavedRows, nPSFVals);
    exit(-1);
  }

  
  // NEW: pass PSF to Convolver object
  psfConvolver.SetupPSF(yVals_psf, nPSFVals);
  

  // NEW: tell Convolver object about size of image
  psfConvolver.SetupProfile(nStoredDataVals);
  

  // NEW: tell Convolver object to finish setup work
  psfConvolver.DoFullSetup(options.debugLevel);
  

  // NEW: tell Convolver object to do the convolution
  psfConvolver.ConvolveProfile(yVals);



  printf("\nSaving output convolved profile (\"%s\") ...\n", options.outputProfileName.c_str());
  SaveProfile(options.outputProfileName, xVals, yVals, nStoredDataVals);


  printf("All done.\n");
  
  free(xVals);
  free(yVals);
  free(xVals_psf);
  free(yVals_psf);
  
  return 0;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   psfconvolve1d input-profile psf-profile [ouput-profile-name]");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine("     --printProfiles            Print out images (for debugging)");
  optParser->AddUsageLine("     --debug <int>            Set debugging output (>= 1)");
//  optParser->AddUsageLine("     --savepadded             Save zero-padded output image also");
  optParser->AddUsageLine("");


  optParser->AddFlag("help", "h");
  optParser->AddFlag("printProfiles");
  optParser->AddOption("debug");      /* an option (takes an argument), supporting only long form */
//  optParser->AddFlag("savepadded");

  /* parse the command line:  */
  optParser->ParseCommandLine( argc, argv );


  /* Process the results: actual arguments, if any: */
  if (optParser->nArguments() > 0) {
    theOptions->dataFileName = optParser->GetArgument(0);
    if (optParser->nArguments() > 1) {
      theOptions->psfProfileName = optParser->GetArgument(1);
      if (optParser->nArguments() > 2) {
        theOptions->outputProfileName = optParser->GetArgument(2);
      }
    }
  }

  /* Process the results: options */
  if ( optParser->FlagSet("help") ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if (optParser->FlagSet("printProfiles")) {
    theOptions->printProfiles = true;
    theOptions->debugLevel = 2;
  }
  if (optParser->OptionSet("debug")) {
    if (NotANumber(optParser->GetTargetString("debug").c_str(), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: debugging level should be a positive integer!\n\n");
      delete optParser;
      exit(1);
    }
    theOptions->debugLevel = atol(optParser->GetTargetString("debug").c_str());
  }
//   if (optParser->FlagSet("savepadded")) {
//     theOptions->outputPaddedProfile = true;
//   }
  
  delete optParser;

}



void SaveProfile( std::string& fileName, double *xVals, double *yVals, int nDataVals )
{
  FILE  *file_ptr;

  if ((file_ptr = fopen(fileName.c_str(), "w")) == NULL) {
    fprintf(stderr, FILE_OPENW_ERR_STRING, fileName.c_str());
    exit(-1);
  }

  for (int i = 0; i < nDataVals; i++) {
    fprintf(file_ptr, "%f\t%f\n", xVals[i], yVals[i]);
  }
  
  fclose(file_ptr);
}



