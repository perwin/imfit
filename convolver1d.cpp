/* FILE: convolver1d.cpp ----------------------------------------------- */
/* VERSION 0.2
 *
 *   Module for profile convolution functions.
 *
 *   MODIFICATION HISTORY:
 *     [v0.01]: 13--14 Aug 2010: Created as modification of convolver.cpp.
 */




/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h"

#include "convolver1d.h"

//using namespace std;



/* ---------------- CONSTRUCTOR ---------------------------------------- */

Convolver1D::Convolver1D( )
{
  psfInfoSet = false;
  profileInfoSet = false;
  fftVectorsAllocated = false;
  fftPlansCreated = false;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

Convolver1D::~Convolver1D( )
{

  if (fftPlansCreated) {
    fftw_destroy_plan(plan_InputProfile);
    fftw_destroy_plan(plan_psf);
    fftw_destroy_plan(plan_inverse);
  }
  if (fftVectorsAllocated) {
    fftw_free(profile_in_cmplx);
    fftw_free(profile_fft_cmplx);
    fftw_free(psf_in_cmplx);
    fftw_free(psf_fft_cmplx);
    fftw_free(multiplied_cmplx);
    fftw_free(convolvedProfile_cmplx);

  }
}


/* ---------------- SetupPSF ------------------------------------------- */
// Pass in a pointer to the pixel vector for the input PSF profiles, as well as
// its size.
void Convolver1D::SetupPSF( double *psfPixels_input, int nPixels )
{

  psfPixels = psfPixels_input;
  nPixels_psf = nPixels;
  psfInfoSet = true;
}


/* ---------------- SetupProfile --------------------------------------- */
// Pass in the size of the data profiles we'll be convolving with the PSF.
void Convolver1D::SetupProfile( int nPixels )
{

  nPixels_data = nPixels;
  profileInfoSet = true;
}


/* ---------------- DoFullSetup ---------------------------------------- */
// General setup prior to actually supplying the profiles data and doing the
// convolution: determine padding size; allocate FFTW arrays and plans;
// normalize, shift, and Fourier transform the PSF profiles.
void Convolver1D::DoFullSetup( int debugLevel, bool doFFTWMeasure )
{
  int  k;
  unsigned  fftwFlags;
  double  psfSum;
  
  debugStatus = debugLevel;
  
  // compute padding dimensions
  if ((! psfInfoSet) || (! profileInfoSet)) {
    printf("*** WARNING: Convolver1D.DoFullSetup: PSF and/or data-profile parameters not set!\n");
    exit(-1);
  }
  nPixels_padded = nPixels_data + nPixels_psf - 1;
  rescaleFactor = 1.0 / nPixels_padded;
  if (debugStatus >= 1)
    printf("Profiles will be padded to %d pixels in size\n", nPixels_padded);

#ifdef FFTW_THREADING
  // TEST: multi-threaded FFTW:
  int  threadStatus;
  threadStatus = fftw_init_threads();
#endif  // FFTW_THREADING

  // allocate memory for fftw_complex arrays
  profile_in_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  profile_fft_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  psf_in_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  psf_fft_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  multiplied_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  convolvedProfile_cmplx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPixels_padded);
  fftVectorsAllocated = true;

  // set up FFTW plans
  if (doFFTWMeasure)
    fftwFlags = FFTW_MEASURE;
  else
    fftwFlags = FFTW_ESTIMATE;
  // Note that there's not much purpose in multi-threading plan_psf, since we only do
  // the FFT of the PSF once
  plan_psf = fftw_plan_dft_1d(nPixels_padded, psf_in_cmplx, psf_fft_cmplx, FFTW_FORWARD,
                             fftwFlags);

#ifdef FFTW_THREADING
  // TEST: multi-threaded FFTW:
  int nThreads = 2;
  fftw_plan_with_nthreads(nThreads);
#endif  // FFTW_THREADING

  plan_InputProfile = fftw_plan_dft_1d(nPixels_padded, profile_in_cmplx, profile_fft_cmplx, FFTW_FORWARD,
                             fftwFlags);
  plan_inverse = fftw_plan_dft_1d(nPixels_padded, multiplied_cmplx, convolvedProfile_cmplx, FFTW_BACKWARD, 
                             fftwFlags);
  fftPlansCreated = true;
  

  // Generate the Fourier transform of the PSF:
  // First, normalize the PSF
  if (debugStatus >= 1) {
    printf("Normalizing the PSF ...\n");
    if (debugStatus >= 2) {
      printf("The whole input PSF profile:\n");
      PrintRealProfile(psfPixels, nPixels_psf);
    }
  }
  psfSum = 0.0;
  for (k = 0; k < nPixels_psf; k++)
    psfSum += psfPixels[k];
  for (k = 0; k < nPixels_psf; k++)
    psfPixels[k] = psfPixels[k] / psfSum;
  if (debugStatus >= 2) {
    printf("The whole *normalized* PSF profile:\n");
    PrintRealProfile(psfPixels, nPixels_psf);
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
    printf("The whole padded, normalized PSF profile:\n");
    PrintComplexProfile_RealPart(psf_in_cmplx, nPixels_padded);
  }
  
  // Finally, do forward FFT on PSF profile
  if (debugStatus >= 1)
    printf("Performing FFT of PSF profile ...\n");
  fftw_execute(plan_psf);
}


/* ---------------- ConvolveProfile ------------------------------------ */
// Given an input profiles (pointer to its pixel vector), convolve it with the PSF
// by: 1) Copying profiles to fft_complex array; 2) Taking FFT of profiles; 3)
// Multiplying transform of profiles by transform of PSF; 4) Taking inverse FFT
// of product; 5) Copying (and rescaling) result back into input profile.
void Convolver1D::ConvolveProfile( double *pixelVector )
{
  int  ii, jj;
  double  a, b, c, d, realPart;
  
  if (debugStatus >= 3) {
    printf("nPixels_data = %d, nPixels_padded = %d\n", nPixels_data, nPixels_padded);
    printf("Original input profile [pixelVector]:\n");
    PrintRealProfile(pixelVector, nPixels_data);
  }
  // Populate (complex) input profiles array for FFT
  //   First, zero the complex array (especially need to do this if this isn't the
  // first time we've called this function!):
  for (ii = 0; ii < nPixels_padded; ii++) {
    profile_in_cmplx[ii][0] = 0.0;
    profile_in_cmplx[ii][1] = 0.0;
  }
  //   Second, copy input profile array into complex array (accounting for padding):
  for (ii = 0; ii < nPixels_data; ii++) {   // step by row number = y
    profile_in_cmplx[ii][0] = pixelVector[ii];
  }
  if (debugStatus >= 3) {
    printf("The whole (padded) input profile [profile_in_cmplx]:\n");
    PrintComplexProfile_RealPart(profile_in_cmplx, nPixels_padded);
  }

  // Do FFT of input profile:
  if (debugStatus >= 2)
    printf("Performing FFT of input profile ...\n");
  fftw_execute(plan_InputProfile);
  if (debugStatus >= 3) {
    printf("The (modulus of the) transform of the input profile [profile_fft_cmplx]:\n");
    PrintComplexProfile_Absolute(profile_fft_cmplx, nPixels_padded);
  }
  
  // Multiply transformed arrays:
  for (jj = 0; jj < nPixels_padded; jj++) {
    a = profile_fft_cmplx[jj][0];   // real part
    b = profile_fft_cmplx[jj][1];   // imaginary part
    c = psf_fft_cmplx[jj][0];
    d = psf_fft_cmplx[jj][1];
    multiplied_cmplx[jj][0] = a*c - b*d;
    multiplied_cmplx[jj][1] = b*c + a*d;
  }
  if (debugStatus >= 3) {
    printf("The (modulus of the) product [multiplied_cmplx]:\n");
    PrintComplexProfile_Absolute(multiplied_cmplx, nPixels_padded);
  }

  // Do the inverse FFT on the product array:
  if (debugStatus >= 2)
    printf("Performing inverse FFT of multiplied profile ...\n");
  fftw_execute(plan_inverse);


  // Extract & rescale the real part of the convolved profile and copy into
  // input pixel vector:
  for (ii = 0; ii < nPixels_data; ii++) {
      realPart = fabs(convolvedProfile_cmplx[ii][0]);
      pixelVector[ii] = rescaleFactor * realPart;
  }
}


// ShiftAndWrapPSF: Takes the input PSF (assumed to be centered in the central pixel
// of the profile) and copy it into the real part of the (padded) fftw_complex profile,
// with the PSF wrapped into the edges, suitable for convolutions.
void Convolver1D::ShiftAndWrapPSF( )
{
  int  centerX_psf;
  int  psfPixel, destPixel;

  centerX_psf = nPixels_psf / 2;
  for (psfPixel = 0; psfPixel < nPixels_psf; psfPixel++) {
    destPixel = (nPixels_padded - centerX_psf + psfPixel) % nPixels_padded;
    psf_in_cmplx[destPixel][0] = psfPixels[psfPixel];
  }
}




void PrintRealProfile( double *data, int nPixels )
{

  for (int i = 0; i < nPixels; i++)
    printf(" %f", data[i]);
  printf("\n");
}


void PrintComplexProfile_RealPart( fftw_complex *data_cmplx, int nPixels )
{

  for (int i = 0; i < nPixels; i++)
    printf(" %9f", data_cmplx[i][0]);
  printf("\n");
}


void PrintComplexProfile_Absolute( fftw_complex *data_cmplx, int nPixels )
{
  double  absVal;
  
  for (int i = 0; i < nPixels; i++) {
    absVal = hypot(data_cmplx[i][0], data_cmplx[i][1]);
    printf(" %9f", absVal);
  }
  printf("\n");
}




/* END OF FILE: convolver1d.cpp ---------------------------------------- */

