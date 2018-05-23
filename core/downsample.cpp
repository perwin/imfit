/* FILE: downsample.cpp -------------------------------------------- */

// Code for downsampling an oversampled image and copying the downsampled version
// into a larger, standard-sampled image.

// Copyright 2014-2018 by Peter Erwin.
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


#include <stdio.h>
#include <string>

//using namespace std;

#include "downsample.h"

/* ------------------- Function Prototypes ----------------------------- */
/* Local Functions: */



/* ---------------- FUNCTION: DownsampleAndReplace() ------------------- */
/// This function takes an image (assumed to be an oversampled sub-region of the
/// main image), downsamples it via block-averaging, and copies the result into
/// the main image.
/// We assume the following about the images:
///    Oversampled image (oversampledImage): nOversampCols x nOversampRows in total
/// size, with assumed padding for (oversampled) PSF-convolution specified by 
/// nOversampPSFCols and nOversampPSFRows.
///    Main image (mainImage): nMainCols x nMainRows in total size, with assumed padding 
/// for (standard) PSF-convolution specified by nMainPSFCols and nMainPSFRows.
///    Target sub-region within mainImage is specified by 1-based coordinates 
/// (x,y) = startX,startY, where (x,y) = (1,1) is the lower-left pixel of the *data*
/// image; this same pixel is at (i,j) = (1 + nMainPSFCols, 1 + nMainPSFRows) within
/// the full (PSF-padded) mainImage.
///    Oversampling scale (oversampleScale) specifies the 1D oversampling, so that
/// each main-size pixel in the sub-region corresponds to oversampleScale x oversampleScale
/// subpixels (i.e., pixels in oversampledImage)

void DownsampleAndReplace( const double *oversampledImage, int nOversampCols, 
						int nOversampRows, int nOversampPSFCols, int nOversampPSFRows,	
						double *mainImage, int nMainCols, int nMainRows, int nMainPSFCols, 
						int nMainPSFRows, int startX, int startY, int oversampleScale, 
						int debugLevel )
{
  int  i, j, i_sub, j_sub, ii, jj;
  int  i1, j1, ii1, jj1;
  int  nCols_subregion, nRows_subregion;
  double  binnedFlux, oversampleArea, oversampPixFlux;
  
  // Coordinate coding:
  //    i,j = 0-based row,column within mainImage (including any PSF padding);
  //          base pixel size
  //    i_sub,j_sub = 0-based row,column within the sub-region of mainImage;
  //          base pixel size
  //    ii,jj = 0-based row,column within oversampledImage (including any PSF padding);
  //          oversampled pixel size
  
  // convert 1-based IRAF-format image-coords to 0-based and account for PSF padding
  // to get starting coords w/in full main image
  j1 = startX - 1 + nMainPSFCols;
  i1 = startY - 1 + nMainPSFRows;
  // get number of columns and rows in sub-region of main image
  nCols_subregion = (int)((nOversampCols - 2*nOversampPSFCols)/oversampleScale);
  nRows_subregion = (int)((nOversampCols - 2*nOversampPSFRows)/oversampleScale);
  oversampleArea = oversampleScale*oversampleScale;

  // iterate over sub-region within main image [i,j = indices for main image]
  if (debugLevel > 1) printf("Starting main loop (with target j1,i1 = %d,%d)...\n", j1,i1);
  for (i = i1; i < i1 + nRows_subregion; i++) {
    if (debugLevel > 1) printf("target row i = %d:\n", i);
    for (j = j1; j < j1 + nCols_subregion; j++) {
      j_sub = j - j1;  // native j (column) index for sub-region
      i_sub = i - i1;  // native i (row) index for sub-region
      if (debugLevel > 1) printf("\ttarget column j = %d: j_sub,i_sub = %d,%d\n", j,j_sub,i_sub);
      // integrate over the oversampleScale x oversampleScale subpixels within this
      // sub-region pixel
      jj1 = j_sub*oversampleScale + nOversampPSFCols;   // starting j-index (columns) within oversampled subimage
      ii1 = i_sub*oversampleScale + nOversampPSFRows;   // starting i-index (rows) within oversampled subimage
      if (debugLevel > 1) printf("\tStarting loop on oversampled image (with osampImage jj1,ii1 = %d,%d):\n", jj1,ii1);
      binnedFlux = 0.0;
      for (ii = ii1; ii < ii1 + oversampleScale; ii++) {
        for (jj = jj1; jj < jj1 + oversampleScale; jj++) {
          oversampPixFlux = oversampledImage[ii*nOversampCols + jj];
          if (debugLevel > 1) printf("\t\toversample pixel at jj,ii = %d,%d: %f\n", jj,ii,oversampPixFlux);
          binnedFlux += oversampledImage[ii*nOversampCols + jj];
        }
      }
      // normalize flux to surface-brightness value for main-image pixels & store in main image
      mainImage[i*nMainCols + j] = binnedFlux/oversampleArea;
    }
  }
  if (debugLevel > 1) printf("Done.\n");

}   



/* END OF FILE: downsample.cpp ------------------------------------- */
