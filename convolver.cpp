/* FILE: convolver.cpp ------------------------------------------------- */
/* VERSION 0.01
 *
 *   Module for image convolution functions.
 *
 *   MODIFICATION HISTORY:
 *     [v0.01]: 26 Mar 2010: Created.
 */

// What we want:
// 
// SETUP:
// 	1. Read in PSF, pass to ModelObject
// 			ModelObject passes PSF vector, size/shape info to Convolver
// 
// 	2. Determine size of model image
// 			ModelObject passes size info and pointer to modelVector to Convolver
// 	
// 	3. Calculate size of padded images
// 	
// 	4. Allocate fftw_complex arrays for
// 			psf_in
// 			psf_fft
// 			image_in
// 			image_fft
// 			multiplied
// 			multiplied_fft [= convolvedData]
// 	
// 	5. Set up FFTW plans
// 			plan_psf
// 			plan_inputImage
// 			plan_inverse
// 
// 	6. Generate FFT(PSF)
// 			A. Normalize PSF
// 			B. ShiftAndWrapPSF()
// 			C. fftw_execute(plan_psf)
// 
// REPEAT FROM MODELOBJECT TILL DONE:
// 	1. Copy modelVector [double] into image_in [fftw_complex]
// 	
// 	2. fftw_execute(plan_inputImage)
// 	
// 	3. Multiply image_fft * psf_fft
// 	
// 	4. fftw_execute(plan_inverse)
// 	
// 	5. Copy & rescale convolved image (multiplied_ff) back into modelVector
// 
// CLEANUP:
// 	1. Clean up FFTW plans:
// 			A. fftw_destroy_plan(plan_inputImage)
// 			B. fftw_destroy_plan(plan_psf)
// 			C. fftw_destroy_plan(plan_inverse)
// 	2. Free fftw_complex arrays:
// 			A. fftw_free(image_in);
// 			B. fftw_free(image_fft);
// 			C. fftw_free(psf_in);
// 			D. fftw_free(psf_fft);
// 			E. fftw_free(multiplied);
// 			F. fftw_free(convolvedData);
// 

/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h"

#include "convolver.h"

//using namespace std;


			
/* ---------------- CONSTRUCTOR ---------------------------------------- */

Convolver::Convolver( )
{
  psfInfoSet = false;
  imageInfoSet = false;
  fftVectorsAllocated = false;
  fftPlansCreated = false;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

Convolver::~Convolver( )
{

  if (fftPlansCreated) {
    fftw_destroy_plan(plan_inputImage);
    fftw_destroy_plan(plan_psf);
    fftw_destroy_plan(plan_inverse);
  }
  if (fftVectorsAllocated) {
    fftw_free(image_in_cmplx);
    fftw_free(image_fft_cmplx);
    fftw_free(psf_in_cmplx);
    fftw_free(psf_fft_cmplx);
    fftw_free(multiplied_cmplx);
    fftw_free(convolvedImage_cmplx);

  }
}


/* ---------------- SetupPSF ------------------------------------------- */
// Pass in a pointer to the pixel vector for the input PSF image, as well as
// the image dimensions
void Convolver::SetupPSF( double *psfPixels_input, int nColumns, int nRows )
{

  psfPixels = psfPixels_input;
  nColumns_psf = nColumns;
  nRows_psf = nRows;
  nPixels_psf = nColumns_psf * nRows_psf;
  psfInfoSet = true;
}


/* ---------------- SetupImage ----------------------------------------- */
// Pass in the dimensions of the image we'll be convolving with the PSF
void Convolver::SetupImage( int nColumns, int nRows )
{

  nColumns_image = nColumns;
  nRows_image = nRows;
  nPixels_image = nColumns_image * nRows_image;
  imageInfoSet = true;
}


/* ---------------- DoFullSetup ---------------------------------------- */
// General setup prior to actually supplying the image data and doing the
// convolution: determine padding dimensions; allocate FFTW arrays and plans;
// normalize, shift, and Fourier transform the PSF image.
void Convolver::DoFullSetup( int debugLevel, bool doFFTWMeasure )
{
  int  k;
  unsigned  fftwFlags;
  double  psfSum;
  
  debugStatus = debugLevel;
  
  // compute padding dimensions
  if ((! psfInfoSet) || (! imageInfoSet)) {
    printf("*** WARNING: Convolver.DoFullSetup: PSF and/or image parameters not set!\n");
    exit(-1);
  }
  nRows_padded = nRows_image + nRows_psf;
  nColumns_padded = nColumns_image + nColumns_psf;
  nPixels_padded = nRows_padded*nColumns_padded;
  rescaleFactor = 1.0 / nPixels_padded;

  // allocate memory for fftw_complex arrays
  image_in_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  image_fft_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  psf_in_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  psf_fft_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  multiplied_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  convolvedImage_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  fftVectorsAllocated = true;

  // set up FFTW plans
  if (doFFTWMeasure)
    fftwFlags = FFTW_MEASURE;
  else
    fftwFlags = FFTW_ESTIMATE;
  plan_psf = fftw_plan_dft_2d(nColumns_padded, nRows_padded, psf_in_cmplx, psf_fft_cmplx, FFTW_FORWARD,
                             fftwFlags);
  plan_inputImage = fftw_plan_dft_2d(nColumns_padded, nRows_padded, image_in_cmplx, image_fft_cmplx, FFTW_FORWARD,
                             fftwFlags);
  plan_inverse = fftw_plan_dft_2d(nColumns_padded, nRows_padded, multiplied_cmplx, convolvedImage_cmplx, FFTW_BACKWARD, 
                             fftwFlags);
  fftPlansCreated = true;
  

  // Generate the Fourier transform of the PSF:
  // First, normalize the PSF
  if (debugStatus > 0) {
    printf("Normalizing the PSF ...\n");
    if (debugStatus > 1) {
      printf("The whole input PSF image, row by row:\n");
      for (int i = 0; i < nRows_psf; i++) {   // step by row number = y
        for (int j = 0; j < nColumns_psf; j++)   // step by column number = x
          printf(" %f", psfPixels[i*nColumns_psf + j]);
        printf("\n");
      }
      printf("\n");
    }
  }
  psfSum = 0.0;
  for (k = 0; k < nPixels_psf; k++)
    psfSum += psfPixels[k];
  for (k = 0; k < nPixels_psf; k++)
    psfPixels[k] = psfPixels[k] / psfSum;
  if (debugStatus > 1) {
    printf("The whole *normalized* PSF image, row by row:\n");
    for (int i = 0; i < nRows_psf; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_psf; j++)   // step by column number = x
        printf(" %f", psfPixels[i*nColumns_psf + j]);
      printf("\n");
    }
    printf("\n");
  }

  // Second, prepare (complex) psf array for FFT, and then copy input PSF into
  // it with appropriate shift/wrap:
  for (k = 0; k < nPixels_padded; k++) {
    psf_in_cmplx[k][0] = 0.0;
    psf_in_cmplx[k][1] = 0.0;
  }
  if (debugStatus > 0)
    printf("Shifting and wrapping the PSF ...\n");
  ShiftAndWrapPSF(psfPixels, nRows_psf, nColumns_psf, psf_in_cmplx, nRows_padded, nColumns_padded);
  if (debugStatus > 1) {
  printf("The whole padded, normalized PSF image, row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %f", psf_in_cmplx[i*nColumns_padded + j][0]);
      printf("\n");
    }
    printf("\n");
  }
  
  // Finally, do forward FFT on PSF image
  if (debugStatus > 0)
    printf("Performing FFT of PSF image ...\n");
  fftw_execute(plan_psf);
}


/* ---------------- ConvolveImage -------------------------------------- */
// Given an input image (pointer to its pixel vector), convolve it with the PSF
// by: 1) Copying image to fft_complex array; 2) Taking FFT of image; 3)
// Multiplying transform of image by transform of PSF; 4) Taking inverse FFT
// of product; 5) Copying (and rescaling) result back into input image.
void Convolver::ConvolveImage( double *pixelVector )
{
  int  ii, jj;
  double  a, b, c, d, realPart;
  
  // Populate (complex) input image array for FFT
  //   First, zero the complex array (especially need to do this if this isn't the
  // first time we've called this function!):
  for (ii = 0; ii < nPixels_padded; ii++) {
    image_in_cmplx[ii][0] = 0.0;
    image_in_cmplx[ii][1] = 0.0;
  }
  //   Second, copy input image array into complex array (accounting for padding):
  for (ii = 0; ii < nRows_image; ii++) {   // step by row number = y
    for (jj = 0; jj < nColumns_image; jj++) {  // step by column number = x
      image_in_cmplx[ii*nColumns_padded + jj][0] = pixelVector[ii*nColumns_image + jj];
    }
  }
  if (debugStatus > 1) {
    printf("The whole (padded) input mage [image_in_cmplx], row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %10f", image_in_cmplx[i*nColumns_padded + j][0]);
      printf("\n");
    }
    printf("\n");
  }

  // Do FFT of input image:
  if (debugStatus > 0)
    printf("Performing FFT of input image ...\n");
  fftw_execute(plan_inputImage);
  
  // Multiply transformed arrays:
  for (jj = 0; jj < nPixels_padded; jj++) {
    a = image_fft_cmplx[jj][0];   // real part
    b = image_fft_cmplx[jj][1];   // imaginary part
    c = psf_fft_cmplx[jj][0];
    d = psf_fft_cmplx[jj][1];
    multiplied_cmplx[jj][0] = a*c - b*d;
    multiplied_cmplx[jj][1] = b*c + a*d;
  }

  // Do the inverse FFT on the product array:
  if (debugStatus > 0)
    printf("Performing inverse FFT of multiplied image ...\n");
  fftw_execute(plan_inverse);

  if (debugStatus > 1) {
    printf("The whole (padded) convolved image [convolvedImage_cmplx, rescaled], row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %10f", fabs(convolvedImage_cmplx[i*nColumns_padded + j][0] / nPixels_padded));
      printf("\n");
    }
    printf("\n");
  }

  // Extract & rescale the real part of the convolved image and copy into
  // input pixel vector:
  for (ii = 0; ii < nRows_image; ii++) {   // step by row number = y
    for (jj = 0; jj < nColumns_image; jj++) {  // step by column number = x
      realPart = fabs(convolvedImage_cmplx[ii*nColumns_padded + jj][0]);
      pixelVector[ii*nColumns_image + jj] = rescaleFactor * realPart;
    }
  }
}


// ShiftAndWrapPSF: Takes an input PSF (assumed to be centered in the central pixel
// of the image) and copies it into the real part of an fftw_complex image, with the
// PSF wrapped into the corners, suitable for convolutions.
// Call with, e.g.
//   ShiftAndWrapPSF(psfPixels, nRows_psf, nColumns_psf, psf_in_cmplx, nRows_padded, nColumns_padded);
void Convolver::ShiftAndWrapPSF( double *psfImage, int nRows_psf, int nCols_psf,
                      fftw_complex *destImage, int nRows_dest, int nCols_dest )
{
  int  centerX_psf = nCols_psf / 2;
  int  centerY_psf = nRows_psf / 2;
  int  psfCol, psfRow, destCol, destRow;
  int  pos_in_psf, pos_in_dest;
  int  i, j;

  for (i = 0; i < nRows_psf; i++) {
    for (j = 0; j < nCols_psf; j++) {
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




/* END OF FILE: convolver.cpp ------------------------------------------ */

