// Code for estimating how much memory will be needed by imfit

#include <math.h>
#include <stdio.h>


const int  FFTW_SIZE = 16;
const int  DOUBLE_SIZE = 8;


/* ------------------- Function Prototypes ----------------------------- */
long EstimateConvolverMemory( int nData_cols, int nData_rows, int nPSF_cols, int nPSF_rows );



long EstimateConvolverMemory( int nModel_cols, int nModel_rows, int nPSF_cols, int nPSF_rows )
{
  long  nPaddedPixels = 0;
  long  nPaddedPixels_cmplx = 0;
  long  nBytesNeeded = 0;
  int  nCols_padded, nRows_padded, nCols_padded_trimmed;

  nBytesNeeded += (long)(nPSF_cols * nPSF_rows);   // allocated outside
  // Convolver stuff
  nCols_padded = nModel_cols + nPSF_cols - 1;
  nRows_padded = nModel_rows + nPSF_rows - 1;
  nCols_padded_trimmed = (int)(floor(nCols_padded/2)) + 1;   // reduced size of r2c/c2r complex array
  nPaddedPixels = (long)(nCols_padded * nRows_padded);
  nPaddedPixels_cmplx = (long)(nCols_padded_trimmed * nRows_padded);
  // 3 double-precision arrays allocated in Convolver::DoFullSetup:
  nBytesNeeded += (long)(3 * nPaddedPixels * DOUBLE_SIZE);
  // 3 fftw_complex arrays allocated in Convolver::DoFullSetup:
  nBytesNeeded += (long)(3 * nPaddedPixels_cmplx * FFTW_SIZE);
  
  return nBytesNeeded;
}


long EstimateMemoryUse( int nData_cols, int nData_rows, int nPSF_cols, int nPSF_rows,
						int nFreeParams, bool levMarFit, bool cashTerms, bool outputResidual,
						bool outputModel, int nPSF_osamp_cols, int nPSF_osamp_rows,
						int deltaX_osamp, int deltaY_osamp, int oversampleScale )
{
  long  nBytesNeeded = 0.0;
  long  nPaddedPixels = 0;
  long  nPaddedPixels_cmplx = 0;
  int  nCols_padded = 0;
  int  nRows_padded = 0;
  int  nCols_padded_trimmed = 0;
  long  nDataPixels = (long)(nData_cols * nData_rows);
  long  dataSize = (long)(nDataPixels * DOUBLE_SIZE);
  int  nModel_rows, nModel_cols;
  int  nOversampModel_cols, nOversampModel_rows;
  long  nModelPixels = 0;
  long  nOversampModelPixels = 0;
  
  nBytesNeeded += dataSize;   // allocated outside
  
  if (nPSF_cols > 0) {
    // we're doing PSF convolution, so model image will be larger
    nModel_cols = nData_cols + 2*nPSF_cols;
    nModel_rows = nData_rows + 2*nPSF_rows;
    nModelPixels = (long)(nModel_cols * nModel_rows);
    // memory used by Convolver object
    nBytesNeeded += EstimateConvolverMemory(nModel_cols, nModel_rows, nPSF_cols, nPSF_rows);
    // possible extra memory used by convolution with oversampled PSF
    if (nPSF_osamp_cols > 0) {
      nOversampModel_cols = deltaX_osamp*oversampleScale + 2*nPSF_osamp_cols;
      nOversampModel_rows = deltaY_osamp*oversampleScale + 2*nPSF_osamp_rows;
      nOversampModelPixels = (long)(nOversampModel_cols * nOversampModel_rows);
      // memory for oversampled model image
      nBytesNeeded += (long)(nOversampModelPixels * DOUBLE_SIZE);
      // memory used by Convolver object for oversampled convolution
      nBytesNeeded += EstimateConvolverMemory(nOversampModel_cols, nOversampModel_rows, 
      											nPSF_osamp_cols, nPSF_osamp_rows);
    }
  }
  else
    nModelPixels = nDataPixels;
  long  modelSize = (long)(nModelPixels * DOUBLE_SIZE); 
  // the following are always allocated
  nBytesNeeded += 3*modelSize;   // modelVector, weightVector, maskVector
  // possible allocations, depending on type of fit and/or outputs requested
  int  nDataSizeAllocs = 0;
  if (levMarFit) {
    nDataSizeAllocs += 3;   // ModelObject's deviatesVector + 2 allocations [fvec, wa4] w/in mpfit.cpp
    nDataSizeAllocs += nFreeParams;   // jacobian array fjac allocated w/in mpfit.cpp
  }
  if (cashTerms)
    nDataSizeAllocs += 1;
  if (outputResidual)
    nDataSizeAllocs += 1;
  if (outputModel)
    nDataSizeAllocs += 1;
  nBytesNeeded += nDataSizeAllocs * dataSize;
  return nBytesNeeded;
}
