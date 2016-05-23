/* FILE: oversampled_region.cpp ---------------------------------------- */
/* 
 *   Module for computing oversampled image regions, and (optionally) convolving them
 * with oversampled PSF.
 *
 *   MODIFICATION HISTORY:
 *     [v0.01]: 29 July 2014: Created.
 */

// Outline for how one of these objects should be set up (i.e., by external code using them):
//   1. theOsampRegion = new OversampledRegion();
//      [optional: theOsampRegion->SetMaxThreads(...)
//   2. theOsampRegion->AddPSFVector( ... )
//   3. theOsampRegion->SetupModelImage( ... )


// For reference: utline for how we supply info to ModelObject in makeimage_main.cpp:
//  theModel = new ModelObject();
//  // Put limits on number of FFTW and OpenMP threads, if user requested it
//  theModel->SetMaxThreads(options.maxThreads);
//  status = AddFunctions(theModel, functionList, functionSetIndices, options.subsamplingFlag);
//  theModel->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels);
//  /* Define the size of the requested model image */
//  theModel->SetupModelImage(nColumns, nRows);

// For reference: utline for how ModelObject sets ups Convolver object
//   psfConvolver->SetupPSF(psfPixels, nColumns_psf, nRows_psf);
//   psfConvolver->SetMaxThreads(maxRequestedThreads);
// 
//   if (doConvolution) {
//     nModelColumns = nDataColumns + 2*nPSFColumns;
//     nModelRows = nDataRows + 2*nPSFRows;
//     psfConvolver->SetupImage(nModelColumns, nModelRows);
//     // NOTE: for now we're ignoring the status of psfConvolver->DoFullSetup because
//     // we assume that it can't fail (we give psfConvolver the PSF info before
//     // setting doConvolution to true, and we give it the image info in the line above)
//     result = psfConvolver->DoFullSetup(debugLevel);
//     nModelVals = nModelColumns*nModelRows;
//   }


// Copyright 2014, 2015 by Peter Erwin.
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




/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>

#include "convolver.h"
#include "function_objects/function_object.h"
#include "oversampled_region.h"
#include "downsample.h"
#ifdef DEBUG
#include "image_io.h"
#endif

using namespace std;


/* ---------------- Definitions ---------------------------------------- */

// current best size for OpenMP processing
#define DEFAULT_OPENMP_CHUNK_SIZE  10


			
/* ---------------- CONSTRUCTOR ---------------------------------------- */

OversampledRegion::OversampledRegion( )
{

  doConvolution = false;
  modelVectorAllocated = false;
  setupComplete = false;
  debugLevel = 0;
  maxRequestedThreads = 0;   // default value --> use all available processors/cores
  ompChunkSize = DEFAULT_OPENMP_CHUNK_SIZE;
  
  debugImageName = "oversampled_region_testoutput";
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

OversampledRegion::~OversampledRegion( )
{
  if (modelVectorAllocated)
    free(modelVector);
  if (doConvolution)
    delete psfConvolver;
}



/* ---------------- SetDebugImageName ---------------------------------- */
// Primarily for debugging purpooses
void OversampledRegion::SetDebugImageName( const string imageName )
{
  debugImageName = imageName;
}


/* ---------------- SetDebugLevel ------------------------------------- */
void OversampledRegion::SetDebugLevel( const int debuggingLevel )
{
  debugLevel = debuggingLevel;
}


/* ---------------- SetMaxThreads -------------------------------------- */
// User specifies maximum number of FFTW threads to use (ignored if not compiled
// with multithreaded FFTW library)
void OversampledRegion::SetMaxThreads( const int maximumThreadNumber )
{
  maxRequestedThreads = maximumThreadNumber;
}


/* ---------------- SetupPSF ------------------------------------------- */
// Pass in a pointer to the pixel vector for the input PSF image, as well as
// the image dimensions.
// We assume this is an oversampled PSF, with the same oversampling scale as
// specified in the input to SetupModelImage(), below
void OversampledRegion::AddPSFVector( double *psfPixels, const int nColumns_psf, 
										const int nRows_psf )
{
  assert( (nColumns_psf >= 1) && (nRows_psf >= 1) );
  
  if (setupComplete) {
    fprintf(stderr,"OverSampledRegion::SetupPSF -- WARNING: the function must be called");
    fprintf(stderr," *before* calling OversampledRegion::SetupModelImage()!\n");
    fprintf(stderr,"Attempt add PSF ignored -- oversampled region calculations will NOT use PSF convolution!\n");
  } else {
    nPSFColumns = nColumns_psf;
    nPSFRows = nRows_psf;
    psfConvolver = new Convolver();
    psfConvolver->SetupPSF(psfPixels, nColumns_psf, nRows_psf);
    psfConvolver->SetMaxThreads(maxRequestedThreads);
    doConvolution = true;
  }
}


/* ---------------- SetupModelImage ------------------------------------ */
// Pass in the dimensions of the image region, oversample scale, etc.
//    x1,y1 = x,y location of lower-left corner of image region w/in main image (IRAF-numbering)
//    nBaseColumns,nBaseRows = x,y size of region in main ("base") image
//    nColumnsMain, nRowsMain = x,y size of full main model ("base") image
int OversampledRegion::SetupModelImage( int x1, int y1, int nBaseColumns, int nBaseRows, 
						int nColumnsMain, int nRowsMain, int nColumnsPSF_main,
						int nRowsPSF_main, int oversampScale )
{
  int  result = 0;
  assert( (nBaseColumns >= 1) && (nBaseRows >= 1) && (oversampScale >= 1) );
  assert( (nColumnsMain >= 1) && (nRowsMain >= 1) );
  assert( (nColumnsPSF_main >= 0) && (nRowsPSF_main >= 0) );
  
  // info about main image (including where LL corner of region is within main image)
  x1_region = x1;
  y1_region = y1;
  nMainImageColumns = nColumnsMain;
  nMainImageRows = nRowsMain;
  nMainPSFColumns = nColumnsPSF_main;
  nMainPSFRows = nRowsPSF_main;
  
  // oversampling info and setup
  oversamplingScale = oversampScale;
  subpixFrac = 1.0 / oversamplingScale;   // linear size of oversampled pixel relative to main-image pixel
  startX_offset = 0.5*subpixFrac - 0.5;
  startY_offset = 0.5*subpixFrac - 0.5;

  // compute dimensions of oversampled region model image
  nRegionColumns = nBaseColumns*oversamplingScale;
  nRegionRows = nBaseRows*oversamplingScale;
  nRegionVals = nRegionColumns*nRegionRows;
  
  if (doConvolution) {
    nModelColumns = nRegionColumns + 2*nPSFColumns;
    nModelRows = nRegionRows + 2*nPSFRows;
    psfConvolver->SetupImage(nModelColumns, nModelRows);
    result = psfConvolver->DoFullSetup(debugLevel);
    if (result < 0) {
      fprintf(stderr, "*** Error returned from Convolver::DoFullSetup!\n");
      return result;
    }
    nModelVals = nModelColumns*nModelRows;
  }
  else {
    nModelColumns = nRegionColumns;
    nModelRows = nRegionRows;
    nModelVals = nRegionVals;
  }
  
  // Allocate modelimage vector
  // WARNING: Possible memory leak (if this function is called more than once)!
  //    If this function *is* called again, then nModelVals could be different
  //    from the first call, in wich case we'd need to realloc modelVector
  modelVector = (double *) calloc((size_t)nModelVals, sizeof(double));
  if (modelVector == NULL) {
    fprintf(stderr, "*** ERROR: Unable to allocate memory for oversampled model image!\n");
    fprintf(stderr, "    (Requested image size was %d pixels)\n", nModelVals);
    return -1;
  }
  modelVectorAllocated = true;
  setupComplete = true;
  
  return 0;
}


/* ---------------- ComputeRegionAndDownsample ------------------------- */
// This is the main method, which computes the oversampled (sub-region) model image,
// then downsamples it to the main image pixel scale and copies it into the main
// image (mainImageVector).
// We assume that the FunctionObjects pointed to by functionObjectVect have
// already been set up with the current parameter values by calling their
// individual Setup() methods -- e.g., by the method or function that is calling
// *this* method.

void OversampledRegion::ComputeRegionAndDownsample( double *mainImageVector, 
					vector<FunctionObject *> functionObjectVect, int nFunctions  )
{
  int   i, j, n, status;
  double  x, y, newValSum, tempSum, adjVal, storedError;
  string  outputName;

// Compute oversampled-region image, using OpenMP for speed
// (possibly slower if sub-region is really small, but in that case this whole
// function will only take a small part of total runtime)
#pragma omp parallel private(i,j,n,x,y,newValSum,tempSum,adjVal,storedError)
  {
  #pragma omp for schedule (static, ompChunkSize)
  // single-loop code which is ~ same in general case as double-loop, and
  // faster for case of small image + many cores (André Luiz de Amorim suggestion)
  for (int k = 0; k < nModelVals; k++) {
    j = k % nModelColumns;
    i = k / nModelColumns;
    y = y1_region + startY_offset + (i - nPSFRows)*subpixFrac;              // Iraf counting: first row = 1
                                                 // (note that nPSFRows = 0 if not doing PSF convolution)
    x = x1_region + startX_offset + (j - nPSFColumns)*subpixFrac;           // Iraf counting: first column = 1
                                                 // (note that nPSFColumns = 0 if not doing PSF convolution)
    newValSum = 0.0;
    storedError = 0.0;
    for (n = 0; n < nFunctions; n++) {
      // Use Kahan summation algorithm
      adjVal = functionObjectVect[n]->GetValue(x, y) - storedError;
      tempSum = newValSum + adjVal;
      storedError = (tempSum - newValSum) - adjVal;
      newValSum = tempSum;
    }
    modelVector[i*nModelColumns + j] = newValSum;
  }
  
  } // end omp parallel section


#ifdef DEBUG
  if (debugLevel > 0) {
    vector<string>  imageCommentsList;
    outputName = debugImageName + ".fits";
    printf("\nOversampledRegion::ComputeRegionAndDownsample -- Saving output model image (\"%s\") ...\n", outputName.c_str());
    status = SaveVectorAsImage(modelVector, outputName, nModelColumns, nModelRows, 
    							imageCommentsList);
  }
#endif

  // Do PSF convolution, if requested
  if (doConvolution)
    psfConvolver->ConvolveImage(modelVector);

#ifdef DEBUG
  if (debugLevel > 0) {
    vector<string>  imageCommentsList;
    outputName = debugImageName + "_conv.fits";
    printf("\nOversampledRegion::ComputeRegionAndDownsample -- Saving PSF-convolved output model image (\"%s\") ...\n", outputName.c_str());
    status = SaveVectorAsImage(modelVector, outputName, 
                        nModelColumns, nModelRows, imageCommentsList);
  }
#endif

  // downsample & copy into main image
  DownsampleAndReplace(modelVector, nModelColumns,nModelRows,nPSFColumns,nPSFRows, 
  						mainImageVector, nMainImageColumns,nMainImageRows,nMainPSFColumns,
  						nMainPSFRows, x1_region,y1_region, oversamplingScale, debugLevel);
}



/* END OF FILE: oversampled_region.cpp --------------------------------- */

