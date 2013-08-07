/* FILE: convolver.cpp ------------------------------------------------- */
/* VERSION 0.2
 *
 *   Module for image convolution functions.
 *
 *   MODIFICATION HISTORY:
 *     [v0.2]: 31 May/1 June 2010: Fixed bug in dealing with non-square images:
 * calls to fftw_plan_* functions had nRows and nColumns in the wrong order!
 *     [v0.1]: 15 April 2010: More or less usable now (thought not thoroughly tested,
 * and I'm not sure whether the "lowered-edges" peculiarity is still present).
 *     [v0.01]: 26 Mar 2010: Created.
 */

// Copyright 2010, 2011, 2012, 2013 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.



// IMPORTANT NOTES ON ROWS, COLUMNS, ETC.
//
//    nColumns = # of possible x values = size of a row = width of image
//    nRows = # of possible y values = size of a column = height of image
//
//    FFTW more or less requires that multi-dimensional arrays be in row-major
// order, as is standard for C.  This means that an image which is W x H (W columns
// by H rows) is indexed thusly:
//    (column x, row y) = A[y][x] = A[y*Ncols + x]
// In the confusing terminology of the FFTW manual: "the last dimension has the fastest-varying
// index in the array" -- meaning that the "last dimension" corresponds to x-values,
// and thus to the number of *columns*.
//    Thus, the proper way to call fftw_plan_dft_2d(), which is described in the
// FFTW manual as:
//       fftw_plan_dft_2d(int n0, int n1, fftw_complex *in, fftw_complex *out,
//                        int sign, unsigned flags);
// is:
//       fftw_plan_dft_2d(nRows, nColumns, fftw_complex *in, fftw_complex *out,
//                        int sign, unsigned flags);



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

#ifndef USING_SCONS
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h"

#ifdef FFTW_THREADING
#include <unistd.h>
#endif  // FFTW_THREADING

#include "convolver.h"

//using namespace std;


			
/* ---------------- CONSTRUCTOR ---------------------------------------- */

Convolver::Convolver( )
{
  psfInfoSet = false;
  imageInfoSet = false;
  fftVectorsAllocated = false;
  fftPlansCreated = false;
  maxRequestedThreads = 0;   // default value --> use all available processors/cores
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


/* ---------------- SetMaxThreads -------------------------------------- */
// User specifies maximum number of FFTW threads to use (ignored if not compiled
// with multithreaded FFTW library)
void Convolver::SetMaxThreads( int maximumThreadNumber )
{
  maxRequestedThreads = maximumThreadNumber;
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
  nColumns_padded = nColumns_image + nColumns_psf - 1;
  nRows_padded = nRows_image + nRows_psf - 1;
  nPixels_padded = nColumns_padded*nRows_padded;
  rescaleFactor = 1.0 / nPixels_padded;
  if (debugStatus >= 1)
    printf("Images will be padded to %d x %d pixels in size\n", nColumns_padded, nRows_padded);

#ifdef FFTW_THREADING
  int  threadStatus;
  threadStatus = fftw_init_threads();
#endif  // FFTW_THREADING

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
  // Note that there's not much purpose in multi-threading plan_psf, since we only do
  // the FFT of the PSF once
  plan_psf = fftw_plan_dft_2d(nRows_padded, nColumns_padded, psf_in_cmplx, psf_fft_cmplx, FFTW_FORWARD,
                             fftwFlags);

#ifdef FFTW_THREADING
  int  nThreads, nCores;
  nCores = sysconf(_SC_NPROCESSORS_ONLN);
  if (maxRequestedThreads == 0) {
    // Default: 1 thread per available core
    nThreads = nCores;
  } else
    nThreads = maxRequestedThreads;
  if (nThreads < 1)
    nThreads = 1;
  fftw_plan_with_nthreads(nThreads);
#endif  // FFTW_THREADING

  plan_inputImage = fftw_plan_dft_2d(nRows_padded, nColumns_padded, image_in_cmplx, image_fft_cmplx, FFTW_FORWARD,
                             fftwFlags);
  plan_inverse = fftw_plan_dft_2d(nRows_padded, nColumns_padded, multiplied_cmplx, convolvedImage_cmplx, FFTW_BACKWARD, 
                             fftwFlags);
  fftPlansCreated = true;
  

  // Generate the Fourier transform of the PSF:
  // First, normalize the PSF
  if (debugStatus >= 1) {
    printf("Normalizing the PSF ...\n");
    if (debugStatus >= 2) {
      printf("The whole input PSF image, row by row:\n");
      PrintRealImage(psfPixels, nColumns_psf, nRows_psf);
    }
  }
  psfSum = 0.0;
  for (k = 0; k < nPixels_psf; k++)
    psfSum += psfPixels[k];
  for (k = 0; k < nPixels_psf; k++)
    psfPixels[k] = psfPixels[k] / psfSum;
  if (debugStatus >= 2) {
    printf("The whole *normalized* PSF image, row by row:\n");
    PrintRealImage(psfPixels, nColumns_psf, nRows_psf);
  }

  // Second, prepare (complex) psf array for FFT, and then copy input PSF into
  // it with appropriate shift/wrap:
  for (k = 0; k < nPixels_padded; k++) {
    psf_in_cmplx[k][0] = 0.0;
    psf_in_cmplx[k][1] = 0.0;
  }
  if (debugStatus >= 1)
    printf("Shifting and wrapping the PSF ...\n");
  ShiftAndWrapPSF();
  if (debugStatus >= 2) {
    printf("The whole padded, normalized PSF image, row by row:\n");
    PrintComplexImage_RealPart(psf_in_cmplx, nColumns_padded, nRows_padded);
  }
  
  // Finally, do forward FFT on PSF image
  if (debugStatus >= 1)
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
  if (debugStatus >= 3) {
    printf("The whole (padded) input mage [image_in_cmplx], row by row:\n");
    PrintComplexImage_RealPart(image_in_cmplx, nColumns_padded, nRows_padded);
  }

  // Do FFT of input image:
  if (debugStatus >= 2)
    printf("Performing FFT of input image ...\n");
  fftw_execute(plan_inputImage);
  if (debugStatus >= 3) {
    printf("The (modulus of the) transform of the input image [image_fft_cmplx], row by row:\n");
    PrintComplexImage_Absolute(image_fft_cmplx, nColumns_padded, nRows_padded);
  }
  
  // Multiply transformed arrays:
  for (jj = 0; jj < nPixels_padded; jj++) {
    a = image_fft_cmplx[jj][0];   // real part
    b = image_fft_cmplx[jj][1];   // imaginary part
    c = psf_fft_cmplx[jj][0];
    d = psf_fft_cmplx[jj][1];
    multiplied_cmplx[jj][0] = a*c - b*d;
    multiplied_cmplx[jj][1] = b*c + a*d;
  }
  if (debugStatus >= 3) {
    printf("The (modulus of the) product [multiplied_cmplx], row by row:\n");
    PrintComplexImage_Absolute(multiplied_cmplx, nColumns_padded, nRows_padded);
  }

  // Do the inverse FFT on the product array:
  if (debugStatus >= 2)
    printf("Performing inverse FFT of multiplied image ...\n");
  fftw_execute(plan_inverse);

  if (debugStatus >= 3) {
    printf("The whole (padded) convolved image [convolvedImage_cmplx, rescaled], row by row:\n");
    for (int i = 0; i < nRows_padded; i++) {   // step by row number = y
      for (int j = 0; j < nColumns_padded; j++)   // step by column number = x
        printf(" %9f", fabs(convolvedImage_cmplx[i*nColumns_padded + j][0] / nPixels_padded));
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


// ShiftAndWrapPSF: Takes the input PSF (assumed to be centered in the central pixel
// of the image) and copy it into the real part of the (padded) fftw_complex image, with the
// PSF wrapped into the corners, suitable for convolutions.
void Convolver::ShiftAndWrapPSF( )
{
  int  centerX_psf, centerY_psf;
  int  psfCol, psfRow, destCol, destRow;
  int  pos_in_psf, pos_in_dest;
  int  i, j;

  centerX_psf = nColumns_psf / 2;
  centerY_psf = nRows_psf / 2;
  for (i = 0; i < nRows_psf; i++) {
    for (j = 0; j < nColumns_psf; j++) {
      psfCol = j;
      psfRow = i;
      pos_in_psf = i*nColumns_psf + j;
      destCol = (nColumns_padded - centerX_psf + psfCol) % nColumns_padded;
      destRow = (nRows_padded - centerY_psf + psfRow) % nRows_padded;
      pos_in_dest = destRow*nColumns_padded + destCol;
      psf_in_cmplx[pos_in_dest][0] = psfPixels[pos_in_psf];
    }
  }
}



void PrintRealImage( double *image, int nColumns, int nRows )
{

  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++)   // step by column number = x
      printf(" %f", image[i*nColumns + j]);
    printf("\n");
  }
  printf("\n");
}


void PrintComplexImage_RealPart( fftw_complex *image_cmplx, int nColumns, int nRows )
{

  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++)   // step by column number = x
      printf(" %9f", image_cmplx[i*nColumns + j][0]);
    printf("\n");
  }
  printf("\n");
}


void PrintComplexImage_Absolute( fftw_complex *image_cmplx, int nColumns, int nRows )
{
  double  absVal;
  
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++) {   // step by column number = x
      absVal = hypot(image_cmplx[i*nColumns + j][0], image_cmplx[i*nColumns + j][1]);
      printf(" %9f", absVal);
    }
    printf("\n");
  }
  printf("\n");
}




/* END OF FILE: convolver.cpp ------------------------------------------ */

