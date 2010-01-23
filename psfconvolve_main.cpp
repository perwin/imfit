// Test code for image convolution

// NOTE: GCC support for C99 "complex" types is "broken" as of version 4.4
// (and for all earlier versions, including 4.2).

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "fftw3.h"

#include "image_io.h"
#include "anyoption.h"   // Kishan Thomas' class for command-line option parsing


/* ---------------- Definitions ---------------------------------------- */
#define INPUT_IMAGE_FILENAME   "narrow_gaussian_6x6.fits"
#define PSF_FILENAME    "gaussian_5x5.fits"
#define DEFAULT_OUTPUT_FILENAME   "convolve_out.fits"


typedef struct {
  std::string  inputImageName;
  std::string  psfImageName;
  std::string  outputImageName;
  bool  printImages;
  bool  psfNormalized;
  bool  outputPaddedImage;
} commandOptions;


/* ------------------- Function Prototypes ----------------------------- */
void ProcessInput( int argc, char *argv[], commandOptions *theOptions );
void ShiftAndWrapPSF( double *psfImage, int nRows_psf, int nCols_psf,
                      fftw_complex *destImage, int nRows_dest, int nCols_dest );



/* ---------------- MAIN ----------------------------------------------- */

int main(int argc, char *argv[])
{
  int  nPixels_input, nPixels_psf, nPixels_padded;
  int  nRows, nColumns;
  int  nRows_padded, nColumns_padded;
  int  nRows_psf, nColumns_psf;
  int  ii, jj;
  std::string  imageFilename, psfFilename, outputFilename;
  std::string  outputPaddedFilename = "convolve_out_padded.fits";
  double  *allPixels;
  double  *psfPixels;
  double  *convolvedData_real, *convolvedData_padded;
  double  psfSum = 0.0;
  double  a, b, c, d, realPart, rescaleFactor;
  fftw_complex  *image_in, *image_fft, *psf_in, *psf_fft;
  fftw_complex  *multiplied, *convolvedData;
  fftw_plan  plan_inputImage, plan_psf, plan_inverse;
  commandOptions  options;

//   double  allPixels[36] =   {0, 0, 0, 0, 0, 0,
//                              0, 0, 0, 0, 0, 0,
//                              0, 0, 1, 0, 0, 0,
//                              0, 0, 0, 0, 0, 0,
//                              0, 0, 0, 0, 0, 0,
//                              0, 0, 0, 0, 0, 0};

  double  psfPixels_simple[25] =   {0,   0,    0.1, 0,    0,
                             0,   0.25, 0.5, 0.25, 0,
                             0.1, 0.5,  1.0, 0.5,  0.1,
                             0,   0.25, 0.5, 0.25, 0,
                             0,   0,    0.1, 0,    0};
  double  psfPad[100] =  {1.0,  0.5,  0.1,  0,   0,  0,  0,  0, 0.1,  0.5,
                          0.5, 0.25,    0,  0,   0,  0,  0,  0,   0, 0.25,
                          0.1,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                            0,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                            0,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                            0,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                            0,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                            0,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                          0.1,    0,    0,  0,   0,  0,  0,  0,   0,    0,
                          0.5, 0.25,    0,  0,   0,  0,  0,  0,   0, 0.25};


  /* Process command line: */
  options.inputImageName = INPUT_IMAGE_FILENAME;
  options.psfImageName = PSF_FILENAME;
  options.outputImageName = DEFAULT_OUTPUT_FILENAME;
  options.printImages = false;
  options.psfNormalized = false;
  options.outputPaddedImage = false;

  ProcessInput(argc, argv, &options);


  printf("\nReading input image (\"%s\") ...\n", options.inputImageName.c_str());
  allPixels = ReadImageAsVector(options.inputImageName, &nColumns, &nRows);
//   nColumns = 6;
//   nRows = 6;
  nPixels_input = nColumns * nRows;
  printf("naxis1 [# pixels/row] = %ld, naxis2 [# pixels/col] = %ld; nPixels_input = %ld\n", 
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

//  psfFilename = PSF_FILENAME;
  printf("Reading PSF image (\"%s\") ...\n", options.psfImageName.c_str());
  psfPixels = ReadImageAsVector(options.psfImageName, &nColumns_psf, &nRows_psf);
//  nColumns_psf = 5;
//  nRows_psf = 5;
  nPixels_psf = nColumns_psf * nRows_psf;
  printf("naxis1 [# pixels/row] = %ld, naxis2 [# pixels/col] = %ld; nPixels_tot = %ld\n", 
           nColumns_psf, nRows_psf, nPixels_psf);

  if (options.printImages) {
    printf("The whole PSF image, row by row:\n");
    for (int i = 0; i < nRows_psf; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_psf; j++)   // step by column number = x
        printf(" %f", psfPixels[i*nColumns_psf + j]);
      printf("\n");
    }
    printf("\n");
  }
  
  /* Normalize the PSF, if it's not already */
  if (! options.psfNormalized) {
    printf("Normalizing the PSF ...\n");
    for (int k = 0; k < nPixels_psf; k++)
      psfSum += psfPixels[k];
    for (int k = 0; k < nPixels_psf; k++)
      psfPixels[k] = psfPixels[k] / psfSum;
  }

  if (options.printImages) {
    printf("The whole *normalized* PSF image, row by row:\n");
    for (int i = 0; i < nRows_psf; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_psf; j++)   // step by column number = x
        printf(" %f", psfPixels[i*nColumns_psf + j]);
      printf("\n");
    }
    printf("\n");
  }
  
  
  /* Figure out how big the zero-padded images should be */
//  nRows_padded = nRows + nRows_psf - 1;
//  nColumns_padded = nColumns + nColumns_psf - 1;
  nRows_padded = nRows + nRows_psf;
  nColumns_padded = nColumns + nColumns_psf;
  nPixels_padded = nRows_padded*nColumns_padded;


  // Setup for FFT work:
  image_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  image_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  psf_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  psf_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  multiplied = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  convolvedData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  
  plan_inputImage = fftw_plan_dft_2d(nColumns_padded, nRows_padded, image_in, image_fft, FFTW_FORWARD,
                             FFTW_ESTIMATE);
  plan_psf = fftw_plan_dft_2d(nColumns_padded, nRows_padded, psf_in, psf_fft, FFTW_FORWARD,
                             FFTW_ESTIMATE);
  plan_inverse = fftw_plan_dft_2d(nColumns_padded, nRows_padded, multiplied, convolvedData, FFTW_BACKWARD, 
                             FFTW_ESTIMATE);


  // Populate (complex) input image array for FFT:
  for (ii = 0; ii < nPixels_padded; ii++) {
    image_in[ii][0] = 0.0;
    image_in[ii][1] = 0.0;
  }
  for (ii = 0; ii < nRows; ii++)
    for (jj = 0; jj < nColumns; jj++) {
      image_in[ii*nColumns_padded + jj][0] = allPixels[ii*nColumns + jj];
    }
  if (options.printImages) {
    printf("The whole padded image, row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %f", image_in[i*nColumns_padded + j][0]);
      printf("\n");
    }
    printf("\n");
  }
  
  // Populate (complex) psf array for FFT:
  for (ii = 0; ii < nPixels_padded; ii++) {
    psf_in[ii][0] = 0.0;
    psf_in[ii][1] = 0.0;
  }
  ShiftAndWrapPSF(psfPixels, nRows_psf, nColumns_psf, psf_in, nRows_padded, nColumns_padded);
  if (options.printImages) {
  printf("The whole padded PSF image, row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %f", psf_in[i*nColumns_padded + j][0]);
      printf("\n");
    }
    printf("\n");
  }
  

  /* Do the forward FFTs: */
  fftw_execute(plan_inputImage);
  fftw_execute(plan_psf);


  /* Multiply the transformed arrays together: */
  for (jj = 0; jj < nPixels_padded; jj++) {
    a = image_fft[jj][0];   // real part
    b = image_fft[jj][1];   // imaginary part
    c = psf_fft[jj][0];
    d = psf_fft[jj][1];
    multiplied[jj][0] = a*c - b*d;
    multiplied[jj][1] = b*c + a*d;
  }

  /* Do the inverse FFT on the product array */
  fftw_execute(plan_inverse);


  /* Extract & rescale the real part of the convolved image */
  convolvedData_real = (double *) calloc(nPixels_input, sizeof(double));
  rescaleFactor = 1.0 / nPixels_padded;
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++) {  // step by column number = x
      realPart = fabs(convolvedData[i*nColumns_padded + j][0]);
      convolvedData_real[i*nColumns + j] = realPart * rescaleFactor;
    }
  }
  if (options.printImages) {
    printf("The whole (padded) convolved image (rescaled), row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %10f", fabs(convolvedData[i*nColumns_padded + j][0] / nPixels_padded));
      printf("\n");
    }
    printf("\n");
  }

  printf("\nSaving output convolved image (\"%s\") ...\n", options.outputImageName.c_str());
  SaveVectorAsImage(convolvedData_real, options.outputImageName, nColumns, nRows);

  if (options.outputPaddedImage) {
    convolvedData_padded = (double *) calloc(nPixels_padded, sizeof(double));
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++) {  // step by column number = x
        realPart = fabs(convolvedData[i*nColumns_padded + j][0]);
        convolvedData_padded[i*nColumns_padded + j] = realPart / nPixels_padded;
      }
    }
    printf("\nSaving padded version of output convolved image (\"%s\") ...\n", options.outputImageName.c_str());
    SaveVectorAsImage(convolvedData_padded, outputPaddedFilename, nColumns_padded, 
                      nRows_padded);
  }


  // Free up FFTW-related stuff
  fftw_destroy_plan(plan_inputImage);
  fftw_destroy_plan(plan_psf);
  fftw_destroy_plan(plan_inverse);
  fftw_free(image_in);
  fftw_free(image_fft);
  fftw_free(psf_in);
  fftw_free(psf_fft);
  fftw_free(multiplied);
  fftw_free(convolvedData);
  
   free(allPixels);
   free(psfPixels);
   free(convolvedData_real);
   if (options.outputPaddedImage)
     free(convolvedData_padded);
  
  return 0;
}



void ShiftAndWrapPSF( double *psfImage, int nRows_psf, int nCols_psf,
                      fftw_complex *destImage, int nRows_dest, int nCols_dest )
{
  int  centerX_psf = nCols_psf / 2;
  int  centerY_psf = nRows_psf / 2;
//  int  centerX_dest = nCols_dest / 2;
//  int  centerY_dest = nRows_dest / 2;
  int  psfCol, psfRow, destCol, destRow;
  int  pos_in_psf, pos_in_dest;

  for (int i = 0; i < nRows_psf; i++) {
    for (int j = 0; j < nCols_psf; j++) {
      psfCol = j;
      psfRow = i;
      pos_in_psf = i*nCols_psf + j;
      destCol = (nCols_dest - centerX_psf + psfCol) % nCols_dest;
      destRow = (nRows_dest - centerY_psf + psfRow) % nRows_dest;
      pos_in_dest = destRow*nCols_dest + destCol;
      destImage[pos_in_dest][0] = psfImage[pos_in_psf];
    }
  }
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  AnyOption *opt = new AnyOption();
  //opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* SET THE USAGE/HELP   */
  opt->addUsage("Usage: ");
  opt->addUsage("   psfconvolve input-image psf-image [ouput-image-name]");
  opt->addUsage(" -h  --help                   Prints this help");
  opt->addUsage("     --printimages            Print out images (for debugging)");
  opt->addUsage("     --normalized             PSF image is already normalized");
  opt->addUsage("     --savepadded             Save zero-padded output image also");
  opt->addUsage("");


  opt->setFlag("help", 'h');
  opt->setFlag("printimages");
  opt->setFlag("normalized");
  opt->setFlag("savepadded");

  /* parse the command line:  */
  opt->processCommandArgs( argc, argv );


  /* Process the results: actual arguments, if any: */
  if (opt->getArgc() > 0) {
    theOptions->inputImageName = opt->getArgv(0);
    if (opt->getArgc() > 1) {
      theOptions->psfImageName = opt->getArgv(1);
      if (opt->getArgc() > 2) {
        theOptions->outputImageName = opt->getArgv(2);
      }
    }
  }

  /* Process the results: options */
  if ( opt->getFlag("help") || opt->getFlag('h') ) {
    opt->printUsage();
    delete opt;
    exit(1);
  }
  if (opt->getFlag("printimages")) {
    theOptions->printImages = true;
  }
  if (opt->getFlag("normalized")) {
    theOptions->psfNormalized = true;
  }
  if (opt->getFlag("savepadded")) {
    theOptions->outputPaddedImage = true;
  }
  
  delete opt;

}



