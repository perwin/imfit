// Code for estimating how much memory will be needed by imfit

// Reminder: the correct way to cast multiplication of int variables when
// assigning to a long variable is:
//    int nCols, nRows;
//    long nPix;
//    ...
//    nPix = (long)nCols * (long)nRows;
//
// (Otherwise, what will happen (with GCC 5, at least) is that the multiplication
// will be done using int values, will overflow, and the overflow value will be
// stored in the long variable...

// Copyright 2015--2017 by Peter Erwin.
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

#include <math.h>
#include <stdio.h>
#include <vector>
#include <tuple>

using namespace std;

#include "psf_oversampling_info.h"


const int  FFTW_SIZE = 16;
const int  DOUBLE_SIZE = 8;


/* ------------------- Function Prototypes ----------------------------- */
long EstimateConvolverMemoryUse( const int nModel_cols, const int nModel_rows, 
								const int nPSF_cols, const int nPSF_rows );



/// Returns an estimate of the number of bytes needed by a Convolver object due
/// to arrays allocated within the object.
long EstimateConvolverMemoryUse( const int nModel_cols, const int nModel_rows, 
								const int nPSF_cols, const int nPSF_rows )
{
  long  nPaddedPixels = 0;
  long  nPaddedPixels_cmplx = 0;
  long  nBytesNeeded = 0;
  int  nCols_padded, nRows_padded, nCols_padded_trimmed;

  nBytesNeeded += (long)nPSF_cols * (long)nPSF_rows;   // allocated outside
  // Convolver stuff
  nCols_padded = nModel_cols + nPSF_cols - 1;
  nRows_padded = nModel_rows + nPSF_rows - 1;
  nCols_padded_trimmed = (int)(floor(nCols_padded/2)) + 1;   // reduced size of r2c/c2r complex array
  nPaddedPixels = (long)nCols_padded * (long)nRows_padded;
  nPaddedPixels_cmplx = (long)nCols_padded_trimmed * (long)nRows_padded;
  // 3 double-precision arrays allocated in Convolver::DoFullSetup:
  nBytesNeeded += 3 * nPaddedPixels * DOUBLE_SIZE;
  // 3 fftw_complex arrays allocated in Convolver::DoFullSetup:
  nBytesNeeded += 3 * nPaddedPixels_cmplx * FFTW_SIZE;
  
  return nBytesNeeded;
}


/// Returns an estimate of the total number of bytes needed due to array allocations
/// for oversampled-PSF convolution within ModelObject (and associated Convolver objects). 
long EstimatePsfOversamplingMemoryUse( vector<PsfOversamplingInfo *> oversamplingInfoVect )
{

  long  nBytesNeeded = 0;
  int  x1, x2, y1, y2;
  
  for (int i = 0; i < (int)oversamplingInfoVect.size(); i++) {
    PsfOversamplingInfo * psfOsampInfo = oversamplingInfoVect[i];
    int  oversampleScale = psfOsampInfo->GetOversamplingScale();
    int  nPSF_osamp_cols = psfOsampInfo->GetNColumns();
    int  nPSF_osamp_rows = psfOsampInfo->GetNRows();
    std::tie(x1, x2, y1, y2) = psfOsampInfo->GetCorrectedRegionCoords();
    int  deltaX = x2 - x1 + 1;
    int  deltaY = y2 - y1 + 1;

    int  nOversampModel_cols = deltaX*oversampleScale + 2*nPSF_osamp_cols;
    int  nOversampModel_rows = deltaY*oversampleScale + 2*nPSF_osamp_rows;
    long  nOversampModelPixels = (long)nOversampModel_cols * (long)nOversampModel_rows;
    // memory for oversampled model image
    nBytesNeeded += nOversampModelPixels * DOUBLE_SIZE;
    // memory used by Convolver object for oversampled convolution
    nBytesNeeded += EstimateConvolverMemoryUse(nOversampModel_cols, nOversampModel_rows, 
   												nPSF_osamp_cols, nPSF_osamp_rows);
  }
  
  return nBytesNeeded;
}


/// Returns an estimate of the total number of bytes needed due to array allocations
/// within ModelObject (and associated Convolver objects), mpfit, and main.
long EstimateMemoryUse( int nData_cols, int nData_rows, int nPSF_cols, int nPSF_rows,
						int nFreeParams, bool levMarFit, bool cashTerms, bool outputResidual,
						bool outputModel )
{
  long  nBytesNeeded = 0.0;
//   long  nPaddedPixels = 0;
//   long  nPaddedPixels_cmplx = 0;
//   int  nCols_padded = 0;
//   int  nRows_padded = 0;
//   int  nCols_padded_trimmed = 0;
  long  nDataPixels = (long)nData_cols * (long)nData_rows;
  long  dataSize = nDataPixels * DOUBLE_SIZE;
  int  nModel_rows, nModel_cols;
//   int  nOversampModel_cols, nOversampModel_rows;
  long  nModelPixels = 0;
//   long  nOversampModelPixels = 0;
  
  nBytesNeeded += dataSize;   // allocated outside
  
  if (nPSF_cols > 0) {
    // we're doing PSF convolution, so model image will be larger
    nModel_cols = nData_cols + 2*nPSF_cols;
    nModel_rows = nData_rows + 2*nPSF_rows;
    nModelPixels = (long)nModel_cols * (long)nModel_rows;
    // memory used by Convolver object
    nBytesNeeded += EstimateConvolverMemoryUse(nModel_cols, nModel_rows, nPSF_cols, nPSF_rows);
  }
  else
    nModelPixels = nDataPixels;
  long  modelSize = nModelPixels * DOUBLE_SIZE; 
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
