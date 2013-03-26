// Test code for image convolution

// NOTE: GCC support for C99 "complex" types is "broken" as of version 4.4
// (and for all earlier versions, including 4.2).

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>

#include "convolver.h"
#include "image_io.h"
#include "commandline_parser.h"


/* ---------------- Definitions ---------------------------------------- */
#define INPUT_IMAGE_FILENAME   "narrow_gaussian_6x6.fits"
#define PSF_FILENAME    "gaussian_5x5.fits"
#define DEFAULT_OUTPUT_FILENAME   "convolve_out.fits"


typedef struct {
  std::string  inputImageName;
  std::string  psfImageName;
  std::string  outputImageName;
  bool  copyHeader;
  bool  printImages;
  bool  outputPaddedImage;
  int  debugLevel;
} commandOptions;


/* ------------------- Function Prototypes ----------------------------- */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );



/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  int  nPixels_input, nPixels_psf;
  int  nRows, nColumns;
  int  nRows_psf, nColumns_psf;
  std::string  imageFilename, psfFilename, outputFilename;
  double  *allPixels;
  double  *psfPixels;
  commandOptions  options;
  Convolver  psfConvolver;


  /* Process command line: */
  options.inputImageName = INPUT_IMAGE_FILENAME;
  options.psfImageName = PSF_FILENAME;
  options.outputImageName = DEFAULT_OUTPUT_FILENAME;
  options.copyHeader = false;
  options.printImages = false;
  options.outputPaddedImage = false;
  options.debugLevel = 0;

  ProcessInput(argc, argv, &options);


  printf("\nReading input image (\"%s\") ...\n", options.inputImageName.c_str());
  allPixels = ReadImageAsVector(options.inputImageName, &nColumns, &nRows);
  nPixels_input = nColumns * nRows;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_input = %d\n", 
           nColumns, nRows, nPixels_input);

  if (options.printImages) {
    printf("The whole image, row by row:\n");
    // The following fetches pixels row-by-row, starting with the bottom
    // row (i.e., what we would normally like to think of as the first row)
    for (int i = 0; i < nRows; i++) {   // step by row number = y
      for (int j = 0; j < nColumns; j++)   // step by column number = x
        printf(" %f", allPixels[i*nColumns + j]);
      printf("\n");
    }
    printf("\n");
  }

  printf("Reading PSF image (\"%s\") ...\n", options.psfImageName.c_str());
  psfPixels = ReadImageAsVector(options.psfImageName, &nColumns_psf, &nRows_psf);
  nPixels_psf = nColumns_psf * nRows_psf;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %d\n", 
           nColumns_psf, nRows_psf, nPixels_psf);

  
  // NEW: pass PSF to Convolver object
  psfConvolver.SetupPSF(psfPixels, nColumns_psf, nRows_psf);
  

  // NEW: tell Convolver object about size of image
  psfConvolver.SetupImage(nColumns, nRows);
  

  // NEW: tell Convolver object to finish setup work
  psfConvolver.DoFullSetup(options.debugLevel);
  

  // NEW: tell Convolver object to do the convolution
  psfConvolver.ConvolveImage(allPixels);



  printf("\nSaving output convolved image (\"%s\") ...\n", options.outputImageName.c_str());
  vector<string>  commentStrings;
  SaveVectorAsImage(allPixels, options.outputImageName, nColumns, nRows, commentStrings);


  
 free(allPixels);
 free(psfPixels);
  
  return 0;
}




void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   psfconvolve input-image psf-image [ouput-image-name]");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine("     --copyheader             Copy input image header to output (convolved) image [NOT YET IMPLEMENTED!]");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --printimages            Print out images (for debugging)");
//  optParser->AddUsageLine("     --savepadded             Save zero-padded output image also");
  optParser->AddUsageLine("");


  optParser->AddFlag("help", "h");
  optParser->AddFlag("copyheader");
  optParser->AddFlag("printimages");
//  optParser->AddFlag("savepadded");

  /* parse the command line:  */
  optParser->ParseCommandLine( argc, argv );


  /* Process the results: actual arguments, if any: */
  int  nArgsFound = optParser->nArguments();
  if (nArgsFound > 0) {
    theOptions->inputImageName = optParser->GetArgument(0);
    if (nArgsFound > 1) {
      theOptions->psfImageName = optParser->GetArgument(1);
      if (nArgsFound > 2) {
        theOptions->outputImageName = optParser->GetArgument(2);
      }
    }
  }

  /* Process the results: options */
  if ( optParser->FlagSet("help") ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if (optParser->FlagSet("copyheader")) {
    theOptions->copyHeader = true;
  }
  if (optParser->FlagSet("printimages")) {
    theOptions->printImages = true;
    theOptions->debugLevel = 2;
  }
//   if (optParser->FlagSet("savepadded")) {
//     theOptions->outputPaddedImage = true;
//   }
  
  delete optParser;

}



