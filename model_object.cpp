/* FILE: model_object.cpp ---------------------------------------------- */
/* VERSION 0.5
 *
 * This is intended to be an abstract base class for the various
 * "model" objects (e.g., image data + fitting functions).
 * 
 *
 * Places where chi^2 or components of chi^2 are calculated:
 *    [] ModelObject::ChiSquared()
 *    [] ModelObject::ComputeDeviates
 *       -- prepares deviates vector, which L-M code then squares and sums
 *
 *
 *     [v0.5]: 16 Apr 2010: Convolution with PSF now works, at least in principle.
 *     [v0.4]: 20--26 Mar 2010: Added stub functions to accomodate PSF image and
 * convolution.
 *     [v0.3.5]: 18--22 Feb 2010: Added PopulateParameterNames() method;
 * this generates the proper set of parameter names (including X0,Y0 for each function set)
 * and is now called by AddFunctions().  Added CheckWeightVector() method, to catch cases
 * where weight vector (whether use-supplied or calculated from Poisson statistics)
 * has "nan" values (or negative values).
 *     [v0.3]:  4 Dec 2009: Added handling of mask images.
 *     [v0.2]: 27 Nov 2009: Modified to include AddDataVectors function, which
 * will be used by derived class ModelObject1D
 *     [v0.1]: 13--15 Nov 2009: Created; initial development.
 *
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



/* ------------------------ Include Files (Header Files )--------------- */

#ifndef USING_SCONS
#include "config.h"
#endif

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "mersenne_twister.h"

#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object.h"
#include "mp_enorm.h"
#include "param_struct.h"


/* ---------------- Definitions ---------------------------------------- */
static string  UNDEFINED = "<undefined>";

// output formatting for printing parameters
#define X0_FORMAT_WITH_ERRS "X0\t\t%.4f # +/- %.4f\n"
#define Y0_FORMAT_WITH_ERRS "Y0\t\t%.4f # +/- %.4f\n"
#define X0_FORMAT "X0\t\t%.4f\n"
#define Y0_FORMAT "Y0\t\t%.4f\n"
#define PARAM_FORMAT_WITH_ERRS "%s\t\t%7g # +/- %7g\n"
#define PARAM_FORMAT "%s\t\t%7g\n"

// very small value for Cash statistic calculations (replaces log(x) if x <= 0)
// Based on http://cxc.harvard.edu/sherpa/ahelp/cstat.html
#define LOG_SMALL_VALUE 1.0e-25

// current best size for OpenMP processing (works well with Intel Core 2 Duo and
// Core i7 in MacBook Pro, under Mac OS X 10.6 and 10.7)
#define OPENMP_CHUNK_SIZE  10


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModelObject::ModelObject( )
{
  dataValsSet = weightValsSet = false;
  parameterBoundsSet = false;
  
  parameterBounds = NULL;
  modelVector = NULL;
  residualVector = NULL;
  outputModelVector = NULL;
  modelVectorAllocated = false;
  weightVectorAllocated = false;
  residualVectorAllocated = false;
  outputModelVectorAllocated = false;
  deviatesVectorAllocated = false;
  doBootstrap = false;
  bootstrapIndicesAllocated = false;
  setStartFlag_allocated = false;
  
  modelErrors = false;
  useCashStatistic = false;
  modelImageComputed = false;
  maskExists = false;
  doConvolution = false;
  zeroPointSet = false;
  
  nFunctions = 0;
  nFunctionSets = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
  debugLevel = 0;
  
  // default image characteristics
  gain = 1.0;
  exposureTime = 1.0;
  nCombined = 1;
  
  maxRequestedThreads = 0;   // default value --> use all available processors/cores
  
  nPSFRows = nPSFColumns = 0;
}


/* ---------------- PUBLIC METHOD: SetDebugLevel ----------------------- */

void ModelObject::SetDebugLevel( int debuggingLevel )
{
  if (debuggingLevel < 0) {
    printf("ModelObject::SetDebugLevel -- WARNING: debugging level must be > 0");
    printf(" (%d was supplied); debugging level left unchanged.\n", debuggingLevel);
  }
  else
    debugLevel = debuggingLevel;
}


/* ---------------- PUBLIC METHOD: SetMaxThreads ----------------------- */

void ModelObject::SetMaxThreads( int maxThreadNumber )
{
  assert( (maxThreadNumber >= 1) );
  maxRequestedThreads = maxThreadNumber;
#ifdef USE_OPENMP
  omp_set_num_threads(maxRequestedThreads);
#endif
}


/* ---------------- PUBLIC METHOD: AddFunction ------------------------- */

void ModelObject::AddFunction( FunctionObject *newFunctionObj_ptr )
{
  int  nNewParams;
  
  functionObjects.push_back(newFunctionObj_ptr);
  nFunctions += 1;
  nNewParams = newFunctionObj_ptr->GetNParams();
  paramSizes.push_back(nNewParams);
  nFunctionParams += nNewParams;
}



/* ---------------- PUBLIC METHOD: DefineFunctionSets ------------------ */

void ModelObject::DefineFunctionSets( vector<int>& functionStartIndices )
{
  int  nn, i;
  
  nFunctionSets = functionStartIndices.size();
    // define array of [false, false, false, ...]
  setStartFlag = (bool *)calloc(nFunctions, sizeof(bool));
  setStartFlag_allocated = true;
  for (i = 0; i < nFunctionSets; i++) {
    nn = functionStartIndices[i];
    // function number n is start of new function set; 
    // change setStartFlag[n] to true
    setStartFlag[nn] = true;
  }
  
  // total number of parameters = number of parameters for individual functions
  // plus x0/y0 pair for each function set
  nParamsTot = nFunctionParams + 2*nFunctionSets;
}



/* ---------------- PUBLIC METHOD: SetZeroPoint ----------------------- */

void ModelObject::SetZeroPoint( double zeroPointValue )
{

  zeroPoint = zeroPointValue;
  zeroPointSet = true;
}


/* ---------------- PUBLIC METHOD: AddImageDataVector ------------------ */

void ModelObject::AddImageDataVector( double *pixelVector, int nImageColumns,
                                      int nImageRows )
{
  nDataVals = nValidDataVals = nImageColumns * nImageRows;
  dataVector = pixelVector;
  dataValsSet = true;
  
  SetupModelImage(nImageColumns, nImageRows);
}


/* ---------------- PUBLIC METHOD: SetupModelImage -------------------- */
// Called by AddImageDataVector(); can also be used by itself in make-image
// mode. Tells ModelObject to allocate space for the model image.
// Note that if PSF convolution is being done, then AddPSFVector() must be
// called *before* this method.
// nImageColumns and nImageRows should refer to the size of the data image
// (in image-fitting mode) OR the requested size of the output model image
// (in make-image mode).
void ModelObject::SetupModelImage( int nImageColumns, int nImageRows )
{
  assert( (nImageColumns >= 1) && (nImageRows >= 1) );
  
  nDataColumns = nImageColumns;
  nDataRows = nImageRows;
  nDataVals = nImageColumns*nImageRows;
  
  if (doConvolution) {
    nModelColumns = nDataColumns + 2*nPSFColumns;
    nModelRows = nDataRows + 2*nPSFRows;
    psfConvolver->SetupImage(nModelColumns, nModelRows);
    psfConvolver->DoFullSetup(debugLevel);
    nModelVals = nModelColumns*nModelRows;
  }
  else {
    nModelColumns = nDataColumns;
    nModelRows = nDataRows;
    nModelVals = nDataVals;
  }
  // Allocate modelimage vector
  modelVector = (double *) calloc((size_t)nModelVals, sizeof(double));
  modelVectorAllocated = true;
}


/* ---------------- PUBLIC METHOD: AddImageCharacteristics ------------ */

void ModelObject::AddImageCharacteristics( double imageGain, double readoutNoise, double expTime, 
									int nCombinedImages, double originalSkyBackground )
{
  assert( (imageGain > 0.0) && (readoutNoise >= 0.0) );
  assert( (expTime > 0.0) && (nCombinedImages >= 1) );
  assert( (originalSkyBackground >= 0.0) );

  gain = imageGain;
  readNoise = readoutNoise;
  exposureTime = expTime;
  nCombined = nCombinedImages;
  originalSky = originalSkyBackground;

  effectiveGain = gain * exposureTime * nCombined;
  readNoise_adu_squared = readNoise*readNoise/(effectiveGain*effectiveGain);
}


/* ---------------- PUBLIC METHOD: AddErrorVector ---------------------- */

void ModelObject::AddErrorVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  assert( (nDataValues == nDataVals) && (nImageColumns == nDataColumns) && 
          (nImageRows == nDataRows) );
  weightVector = pixelVector;
  
  // Convert noise values into weights, if needed.  Our standard approach is
  // to compute weights as 1/sigma; this assumes that whatever function calls
  // ComputeDeviates() will then square the individual (weighted) deviate values
  // in order to get the proper chi^2 result.
  // Currently, we assume three possibilities for weight-map pixel values:
  //    sigma (std.dev.); variance (sigma^2); and plain weights
  //    Note that correct interpretation of chi^2 values depends on weights
  //    being based on sigmas or variances!
  switch (inputType) {
    case WEIGHTS_ARE_SIGMAS:
      for (int z = 0; z < nDataVals; z++) {
        weightVector[z] = 1.0 / weightVector[z];
      }
      break;
    case WEIGHTS_ARE_VARIANCES:
      for (int z = 0; z < nDataVals; z++) {
        weightVector[z] = 1.0 / sqrt(weightVector[z]);
      }
      break;
    default:
      // do nothing, since input values *are* weights (e.g., the input image is
      // a weight map with pixel values = 1/sigma already)
      ;
  }
  
  if (CheckWeightVector())
    weightValsSet = true;
  else {
    printf("ModelObject::AddErrorVector -- Conversion of error vector resulted in bad values!\n");
    printf("Exiting ...\n\n");
    exit(-1);
  }
}


/* ---------------- PUBLIC METHOD: GenerateErrorVector ----------------- */
// Generate an error vector based on the Gaussian approximation of Poisson statistics.
//    noise^2 = object_flux + sky + rdnoise^2
//
// Since sigma_adu = sigma_e/gain, we can go from
//    noise(e-)^2 = object_flux(e-) + sky(e-) + rdnoise^2
// to
//    noise(adu)^2 = object_flux(adu)/gain + sky(adu)/gain + rdnoise^2/gain^2
// or just
//    noise(adu)^2 = (object_flux(adu) + sky(adu))/gain + rdnoise^2/gain^2
// (assuming that read noise is in units of e-, as is usual)
//
// Exposure time and number of averaged images can be accounted for by including them in 
// the effective gain:
//    gain_eff = (gain * t_exp * N_combined)
// HOWEVER, in this case we also have to account for the multiple readouts, which means
// that the read noise term is multiplied by N_combined, so that we end up with
//    noise(adu)^2 = (object_flux(adu) + sky(adu))/gain_eff + N_combined * rdnoise^2/gain_eff^2
// (where "adu" can be adu/sec if t_exp != 1)

void ModelObject::GenerateErrorVector( )
{
  double  noise_squared, totalFlux;

  // Allocate storage for weight image:
  weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  weightVectorAllocated = true;
  
//  readNoise_adu_squared = readNoise*readNoise/(effectiveGain*effectiveGain);
  // Compute noise estimate for each pixel (see above for derivation)
  // Note that we assume a constant sky background (presumably already subtracted)
  for (int z = 0; z < nDataVals; z++) {
    // total flux = nCombined*(source + subtracted_sky_value)
    totalFlux = dataVector[z] + originalSky;
    if (totalFlux < 0.0)
      totalFlux = 0.0;
    noise_squared = totalFlux/effectiveGain + nCombined*readNoise_adu_squared;
    // Note that we store 1/sigma instead of 1/sigma^2, since the chi^2 calculation in 
    // ChiSquared() [or the equivalent in mpfit.cpp) will square the individual terms
    weightVector[z] = 1.0 / sqrt(noise_squared);
  }

  if (CheckWeightVector())
    weightValsSet = true;
  else {
    printf("ModelObject::GenerateErrorVector -- Calculation of error vector resulted in bad values!\n");
    printf("Exiting ...\n\n");
    exit(-1);
  }
}



/* ---------------- PUBLIC METHOD: AddMaskVector ----------------------- */
// Code for adding and processing a vector containing the 2D mask image.
// Note that although our default *input* format is "0 = good pixel, > 0 =
// bad pixel", internally we convert all bad pixels to 0 and all good pixels
// to 1, so that we can multiply the weight vector by the (internal) mask values.
// The mask is applied to the weight vector by calling the ApplyMask() method
// for a given ModelObject instance.
void ModelObject::AddMaskVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  assert( (nDataValues == nDataVals) && (nImageColumns == nDataColumns) && 
          (nImageRows == nDataRows) );

  maskVector = pixelVector;
  nValidDataVals = 0;   // Since there's a mask, not all pixels from the original
                        // image will be valid
    
  // We need to convert the mask values so that good pixels = 1 and bad
  // pixels = 0.
  switch (inputType) {
    case MASK_ZERO_IS_GOOD:
      // This is our "standard" input mask: good pixels are zero, bad pixels
      // are positive integers
      printf("ModelObject::AddMaskVector -- treating zero-valued pixels as good ...\n");
      for (int z = 0; z < nDataVals; z++) {
        if (maskVector[z] > 0.0) {
          maskVector[z] = 0.0;
        } else {
          maskVector[z] = 1.0;
          nValidDataVals++;
        }
      }
      break;
    case MASK_ZERO_IS_BAD:
      // Alternate form for input masks: good pixels are 1, bad pixels are 0
      printf("ModelObject::AddMaskVector -- treating zero-valued pixels as bad ...\n");
      for (int z = 0; z < nDataVals; z++) {
        if (maskVector[z] < 1.0)
          maskVector[z] = 0.0;
        else {
          maskVector[z] = 1.0;
          nValidDataVals++;
        }
      }
      break;
    default:
      printf("ModelObject::AddMaskVector -- WARNING: unknown inputType detected!\n\n");
      exit(-1);
  }
      
  maskExists = true;
}


/* ---------------- PUBLIC METHOD: ApplyMask --------------------------- */

void ModelObject::ApplyMask( )
{
  if ( (weightValsSet) && (maskExists) ) {
    for (int z = 0; z < nDataVals; z++) {
      weightVector[z] = maskVector[z] * weightVector[z];
    }
    printf("ModelObject: mask vector applied to weight vector. ");
    printf("(%d valid pixels remain)\n", nValidDataVals);
  }
  else {
    printf(" ** ALERT: ModelObject::ApplyMask() called, but we are missing either\n");
    printf("    error image or mask image, or both!  ApplyMask() ignored ...\n");
  }
}


/* ---------------- PUBLIC METHOD: AddPSFVector ------------------------ */
// This function is called to pass in the PSF image and dimensions; doing so
// automatically triggers setup of a Convolver object to do convolutions (including
// prep work such as computing the FFT of the PSF).
// This function must be called *before* SetupModelImage() is called (to ensure
// that we know the proper model-image dimensions
void ModelObject::AddPSFVector(int nPixels_psf, int nColumns_psf, int nRows_psf,
                         double *psfPixels)
{
  assert( (nPixels_psf >= 1) && (nColumns_psf >= 1) && (nRows_psf >= 1) );
  
  nPSFColumns = nColumns_psf;
  nPSFRows = nRows_psf;
  psfConvolver = new Convolver();
  psfConvolver->SetupPSF(psfPixels, nColumns_psf, nRows_psf);
  psfConvolver->SetMaxThreads(maxRequestedThreads);
  doConvolution = true;
}



/* ---------------- PUBLIC METHOD: FinalSetup -------------------------- */
void ModelObject::FinalSetup( )
{
  if (maskExists)
    ApplyMask();
  bool dataOK = VetDataVector();
  if (! dataOK) {
    fprintf(stderr, "ERROR: bad (non-masked) data values!\n\n");
    exit(-1);
  }
}



/* ---------------- PUBLIC METHOD: CreateModelImage -------------------- */

void ModelObject::CreateModelImage( double params[] )
{
  double  x0, y0, x, y, newValSum;
  int  i, j, n;
  int  offset = 0;
  
  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    printf("** ModelObject::CreateModelImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g, ", parameterLabels[0].c_str(), params[0]);
    for (int z = 1; z < nParamsTot; z++)
      printf(", %s = %g", parameterLabels[z].c_str(), params[z]);
    printf("\n");
#endif
    printf("Exiting ...\n\n");
    exit(-1);
  }

  // Separate out the individual-component parameters and tell the
  // associated function objects to do setup work.
  // The first component's parameters start at params[0]; the second's
  // start at params[paramSizes[0]], the third at 
  // params[paramSizes[0] + paramSizes[1]], and so forth...
  for (n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0, y0);
    offset += paramSizes[n];
  }
  
  double  tempSum, adjVal, storedError;
  
  // OK, populate modelVector with the model image
  // OpenMP Parallel Section
  int  chunk = OPENMP_CHUNK_SIZE;
// Note that we cannot specify modelVector as shared [or private] bcs it is part
// of a class (not an independent variable); happily, by default all references in
// an omp-parallel section are shared unless specified otherwise
#pragma omp parallel private(i,j,n,x,y,newValSum,tempSum,adjVal,storedError)
  {
  #pragma omp for schedule (static, chunk)
  for (i = 0; i < nModelRows; i++) {   // step by row number = y
    y = (double)(i - nPSFRows + 1);              // Iraf counting: first row = 1 
                                                 // (note that nPSFRows = 0 if not doing PSF convolution)
    for (j = 0; j < nModelColumns; j++) {   // step by column number = x
      x = (double)(j - nPSFColumns + 1);                 // Iraf counting: first column = 1
                                                         // (note that nPSFColumns = 0 if not doing PSF convolution)
      newValSum = 0.0;
      storedError = 0.0;
      for (n = 0; n < nFunctions; n++) {
        // Use Kahan summation algorithm
        adjVal = functionObjects[n]->GetValue(x, y) - storedError;
        tempSum = newValSum + adjVal;
        storedError = (tempSum - newValSum) - adjVal;
        newValSum = tempSum;
      }
      modelVector[i*nModelColumns + j] = newValSum;
    }
  }
  
  } // end omp parallel section
  
  
  // Do PSF convolution, if requested
  if (doConvolution)
    psfConvolver->ConvolveImage(modelVector);
  
  modelImageComputed = true;
}


/* ---------------- PUBLIC METHOD: SingleFunctionImage ----------------- */
// Generate a model image using *one* of the FunctionObjects (the one indicated by
// functionIndex) and the input parameter vector; returns pointer to modelVector.
// If PSF convolution is requested, then a new output modelVector is created and
// returned, excluding the expanded margin used for PSF convolution (so that the
// returned image will be the same size as the data image).
//
// Meant to be called *externally* (i.e., do NOT call this from within ModelObject,
// unless you are aware that it will NOT return the full (expanded) model image.
double * ModelObject::GetSingleFunctionImage( double params[], int functionIndex )
{
  double  x0, y0, x, y, newVal;
  int  i, j, n;
  int  offset = 0;
  int  iDataRow, iDataCol, z, zModel;
  
  assert( (functionIndex >= 0) );
  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    printf("** ModelObject::SingleFunctionImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g, ", parameterLabels[0].c_str(), params[0]);
    for (z = 1; z < nParamsTot; z++)
      printf(", %s = %g", parameterLabels[z].c_str(), params[z]);
    printf("\n");
#endif
    printf("Exiting ...\n\n");
    exit(-1);
  }

  // Separate out the individual-component parameters and tell the
  // associated function objects to do setup work.
  // The first component's parameters start at params[0]; the second's
  // start at params[paramSizes[0]], the third at 
  // params[paramSizes[0] + paramSizes[1]], and so forth...
  for (n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0, y0);
    offset += paramSizes[n];
  }
  
  // OK, populate modelVector with the model image
  // OpenMP Parallel Section
  int  chunk = OPENMP_CHUNK_SIZE;
// Note that we cannot specify modelVector as shared [or private] bcs it is part
// of a class (not an independent variable); happily, by default all references in
// an omp-parallel section are shared unless specified otherwise
#pragma omp parallel private(i,j,n,x,y,newVal)
  {
  #pragma omp for schedule (static, chunk)
  for (i = 0; i < nModelRows; i++) {   // step by row number = y
    y = (double)(i - nPSFRows + 1);              // Iraf counting: first row = 1
    for (j = 0; j < nModelColumns; j++) {   // step by column number = x
      x = (double)(j - nPSFColumns + 1);                 // Iraf counting: first column = 1
      newVal = functionObjects[functionIndex]->GetValue(x, y);
      modelVector[i*nModelColumns + j] = newVal;
    }
  }
  
  } // end omp parallel section
  
  
  // Do PSF convolution, if requested
  if (doConvolution) {
    if (! outputModelVectorAllocated) {
      outputModelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
      outputModelVectorAllocated = true;
    }
    psfConvolver->ConvolveImage(modelVector);
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels output image
    for (z = 0; z < nDataVals; z++) {
      iDataRow = z / nDataColumns;
      iDataCol = z - iDataRow*nDataColumns;
      zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
      outputModelVector[z] = modelVector[zModel];
    }
    return outputModelVector;
  }
  else
    return modelVector;
}


/* ---------------- PUBLIC METHOD: UpdateWeightVector ------------------ */
/* This function computes new error-based weights using the current model
 * image and the Gaussian approximation to Poisson statistics. Used if we
 * are doing chi^2 minimization with sigma estimated from *model* values
 * (Pearson's chi^2).
 */
void ModelObject::UpdateWeightVector(  )
{
  int  iDataRow, iDataCol, z, zModel;
  double  totalFlux, noise_squared;
	
  if (doConvolution) {
    for (z = 0; z < nDataVals; z++) {
      if ( (! maskExists) || (maskExists && (maskVector[z] > 0)) ) {
        // only update values that aren't masked out
        // (don't rely on previous weightVector[z] value, since sometimes model flux
        // might be zero for an unmasked pixel)
        iDataRow = z / nDataColumns;
        iDataCol = z - iDataRow*nDataColumns;
        zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        totalFlux = modelVector[zModel] + originalSky;
        noise_squared = totalFlux/effectiveGain + nCombined*readNoise_adu_squared;
        // POSSIBLE PROBLEM: if originalSky = model flux = read noise = 0, we'll have /0 error!
        weightVector[z] = 1.0 / sqrt(noise_squared);
      }
    }
  }
  else {
    // No convolution, so model image is same size & shape as data and weight images
    for (z = 0; z < nDataVals; z++) {
      if ( (! maskExists) || (maskExists && (maskVector[z] > 0)) ) {
        // only update values that aren't masked out
        // (don't rely on previous weightVector[z] value, since sometimes model flux
        // might be zero for an unmasked pixel)
        totalFlux = modelVector[z] + originalSky;
        noise_squared = totalFlux/effectiveGain + nCombined*readNoise_adu_squared;
        // POSSIBLE PROBLEM: if originalSky = model flux = read noise = 0, we'll have /0 error!
        weightVector[z] = 1.0 / sqrt(noise_squared);
      }
    }
  }  
}


/* ---------------- PUBLIC METHOD: ComputeDeviates --------------------- */
/* This function computes the vector of weighted deviates (differences between
 * model and data pixel values).  Note that a proper chi^2 sum requires *squaring*
 * each deviate value before summing them; we assume this is done elsewhere, by 
 * whatever function calls ComputeDeviates().
 *
 * Primarily for use by Levenberg-Marquardt solver (mpfit.cpp); for standard
 * chi^2 calculations, use ChiSquared().
 */
void ModelObject::ComputeDeviates( double yResults[], double params[] )
{
  int  iDataRow, iDataCol, z, zModel, b, bModel;
  
#ifdef DEBUG
  printf("ComputeDeviates: Input parameters: ");
  for (int n = 0; n < nParamsTot; n++)
    printf("p[%d] = %g, ", n, params[n]);
  printf("\n");
#endif

  CreateModelImage(params);
  if (modelErrors)
    UpdateWeightVector();

  // In standard case, z = index into dataVector, weightVector, and yResults; it comes 
  // from linearly stepping through (0, ..., nDataVals).
  // In the bootstrap case, z = index into yResults and bootstrapIndices vector;
  // b = bootstrapIndices[z] = index into dataVector and weightVector
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in data and weight images (excluding the outer borders which are only
    // for ensuring proper PSF convolution)
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        iDataRow = b / nDataColumns;
        iDataCol = b - iDataRow*nDataColumns;
        bModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        yResults[z] = weightVector[b] * (dataVector[b] - modelVector[bModel]);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        iDataRow = z / nDataColumns;
        iDataCol = z - iDataRow*nDataColumns;
        zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        yResults[z] = weightVector[z] * (dataVector[z] - modelVector[zModel]);
      }
    }
  }
  else {
    // No convolution, so model image is same size & shape as data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        yResults[z] = weightVector[b] * (dataVector[b] - modelVector[b]);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        yResults[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
      }
    }
  }

}


/* ---------------- PUBLIC METHOD: UseModelErrors --------==----------- */

void ModelObject::UseModelErrors( )
{
  modelErrors = true;
  // Allocate storage for weight image (do this here because we assume that
  // AddErrorVector() or GenerateErrorVector() will NOT be called if we're using 
  // Cash statistic), and set all values = 1
  if (weightVectorAllocated) {
    printf("ERROR: ModelImage::UseModelErrors -- weight vector already allocated!\n");
    printf("Exiting ...\n\n");
    exit(-1);
  }
  weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  for (int z = 0; z < nDataVals; z++) {
    weightVector[z] = 1.0;
  }
  weightVectorAllocated = true;
  weightValsSet = true;
}


/* ---------------- PUBLIC METHOD: UseCashStatistic ------------------- */

void ModelObject::UseCashStatistic( )
{
  useCashStatistic = true;

  // Allocate storage for weight image (do this here because we assume that
  // AddErrorVector() or GenerateErrorVector() will NOT be called if we're using 
  // Cash statistic), and set all values = 1
  if (weightVectorAllocated) {
    printf("ERROR: ModelImage::UseCashStatistic -- weight vector already allocated!\n");
    printf("Exiting ...\n\n");
    exit(-1);
  }
  weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  for (int z = 0; z < nDataVals; z++) {
    weightVector[z] = 1.0;
  }
  weightVectorAllocated = true;
  weightValsSet = true;
}


/* ---------------- PUBLIC METHOD: UsingCashStatistic ------------------ */

bool ModelObject::UsingCashStatistic( )
{
  return useCashStatistic;
}


/* ---------------- PUBLIC METHOD: GetFitStatistic --------------------- */
/* Function for calculating chi^2 value for a model.
 *
 */
double ModelObject::GetFitStatistic( double params[] )
{
  if (useCashStatistic)
    return CashStatistic(params);
  else
    return ChiSquared(params);
}


/* ---------------- PUBLIC METHOD: ChiSquared -------------------------- */
/* Function for calculating chi^2 value for a model.
 *
 */
double ModelObject::ChiSquared( double params[] )
{
  int  iDataRow, iDataCol, z, zModel, b, bModel;
  double  chi;
  
  if (! deviatesVectorAllocated) {
//    deviatesVector = (double *) malloc(nDataVals * sizeof(double));
    deviatesVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    deviatesVectorAllocated = true;
  }
  
  CreateModelImage(params);
  if (modelErrors)
    UpdateWeightVector();
  
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        iDataRow = b / nDataColumns;
        iDataCol = b - iDataRow*nDataColumns;
        bModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        deviatesVector[z] = weightVector[b] * (dataVector[b] - modelVector[bModel]);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        iDataRow = z / nDataColumns;
        iDataCol = z - iDataRow*nDataColumns;
        zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        deviatesVector[z] = weightVector[z] * (dataVector[z] - modelVector[zModel]);
      }
    }
  }
  else {
    // Model image is same size & shape as data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        deviatesVector[z] = weightVector[b] * (dataVector[b] - modelVector[b]);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        deviatesVector[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
      }
    }
  }
  
  // Recall that mp_enorm returns sqrt( Sum_i(chi_i^2) ) = sqrt( Sum_i(deviatesVector[i]^2) )
  if (doBootstrap)
    chi = mp_enorm(nValidDataVals, deviatesVector);
  else
    chi = mp_enorm(nDataVals, deviatesVector);
  
  return (chi*chi);
}


/* ---------------- PUBLIC METHOD: CashStatistic ----------------------- */
/* Function for calculating Cash statistic for a model
 *
 */
double ModelObject::CashStatistic( double params[] )
{
  int  iDataRow, iDataCol, z, zModel, b, bModel;
  double  modVal, dataVal, logModel;
  double  cashStat = 0.0;
  
  CreateModelImage(params);
  
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        iDataRow = b / nDataColumns;
        iDataCol = b - iDataRow*nDataColumns;
        bModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        modVal = effectiveGain*(modelVector[bModel] + originalSky);
        dataVal = effectiveGain*(dataVector[b] + originalSky);
        if (modVal <= 0)
          logModel = LOG_SMALL_VALUE;
        else
          logModel = log(modVal);
        cashStat += weightVector[b] * (modVal - dataVal*logModel);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        iDataRow = z / nDataColumns;
        iDataCol = z - iDataRow*nDataColumns;
        zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        // Mi − Di + DilogDi − DilogMi
        modVal = effectiveGain*(modelVector[zModel] + originalSky);
        dataVal = effectiveGain*(dataVector[z] + originalSky);
        if (modVal <= 0)
          logModel = LOG_SMALL_VALUE;
        else
          logModel = log(modVal);
        cashStat += weightVector[z] * (modVal - dataVal*logModel);
      }
    }
  }
  else {
    // Model image is same size & shape as data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        modVal = effectiveGain*(modelVector[b] + originalSky);
        dataVal = effectiveGain*(dataVector[b] + originalSky);
        if (modVal <= 0)
          logModel = LOG_SMALL_VALUE;
        else
          logModel = log(modVal);
        cashStat += weightVector[z] * (modVal - dataVal*logModel);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        modVal = effectiveGain*(modelVector[z] + originalSky);
        dataVal = effectiveGain*(dataVector[z] + originalSky);
        if (modVal <= 0)
          logModel = LOG_SMALL_VALUE;
        else
          logModel = log(modVal);
        cashStat += weightVector[z] * (modVal - dataVal*logModel);
      }
    }
  }
  
  return (2.0*cashStat);
}


/* ---------------- PUBLIC METHOD: PrintDescription ------------------- */

void ModelObject::PrintDescription( )
{
  printf("Model Object: %d data values (pixels)\n", nDataVals);
}


/* ---------------- PUBLIC METHOD: GetFunctionNames ------------------- */

void ModelObject::GetFunctionNames( vector<string>& functionNames )
{
  functionNames.clear();
  for (int n = 0; n < nFunctions; n++) {
    functionNames.push_back(functionObjects[n]->GetShortName());
  }
}


/* ---------------- PUBLIC METHOD: PrintModelParams --------=---------- */
// Basic function which prints to a file (or, e.g. stdout) a summary of the
// best-fitting model, in form suitable for future use as an input config file.
// If parameterInfo != NULL, then x0,y0 are corrected for any offsets.
// If errs != NULL, then +/- errors are printed as well

void ModelObject::PrintModelParams( FILE *output_ptr, double params[], mp_par *parameterInfo,
																		double errs[] )
{
  double  x0, y0, paramVal;
  int nParamsThisFunc, k;
  int  indexOffset = 0;
  string  funcName, paramName;

  for (int n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      k = indexOffset;
      x0 = params[k];
      y0 = params[k + 1];
      if (parameterInfo != NULL) {
        x0 += parameterInfo[k].offset;
        y0 += parameterInfo[k + 1].offset;
      }
      if (errs != NULL) {
        fprintf(output_ptr, X0_FORMAT_WITH_ERRS, x0, errs[k]);
        fprintf(output_ptr, Y0_FORMAT_WITH_ERRS, y0, errs[k + 1]);
      } else {
        fprintf(output_ptr, X0_FORMAT, x0);
        fprintf(output_ptr, Y0_FORMAT, y0);
      }
      indexOffset += 2;
    }
    
    // Now print the function and its parameters
    nParamsThisFunc = paramSizes[n];
    funcName = functionObjects[n]->GetShortName();
    fprintf(output_ptr, "FUNCTION %s\n", funcName.c_str());
    for (int i = 0; i < nParamsThisFunc; i++) {
      paramName = GetParameterName(indexOffset + i);
      paramVal = params[indexOffset + i];
      if (errs != NULL)
        fprintf(output_ptr, PARAM_FORMAT_WITH_ERRS, paramName.c_str(), paramVal, errs[indexOffset + i]);
      else
        fprintf(output_ptr, PARAM_FORMAT, paramName.c_str(), paramVal);
    }
    indexOffset += paramSizes[n];
  }
}


/* ---------------- PUBLIC METHOD: UseBootstrap ------------------------ */
// Tells ModelObject1d object that from now on we'll operate in bootstrap
// resampling mode, so that bootstrapIndices vector is used to access the
// data and model values (and weight values, if any).

void ModelObject::UseBootstrap( )
{
  doBootstrap = true;
  MakeBootstrapSample();
}


/* ---------------- PUBLIC METHOD: MakeBootstrapSample ----------------- */
// Generate a new bootstrap resampling of the data (actually, we generate a
// bootstrap resampling of the data *indices*)

void ModelObject::MakeBootstrapSample( )
{
  int  n;
  bool  badIndex;
  
  if (! bootstrapIndicesAllocated) {
    bootstrapIndices = (int *) calloc((size_t)nValidDataVals, sizeof(int));
    bootstrapIndicesAllocated = true;
  }
  for (int i = 0; i < nValidDataVals; i++) {
    // pick random data point between 0 and nDataVals - 1, inclusive;
    // reject masked pixels
    //n = round( (random()/MAX_RANDF)*(nDataVals - 1) );
    badIndex = true;
    do {
      n = (int)floor( genrand_real2()*nDataVals );
      if (weightVector[n] > 0.0)
        badIndex = false;
    } while (badIndex);
    bootstrapIndices[i] = n;
  }
}




/* ---------------- PUBLIC METHOD: PrintImage ------------------------- */
// Basic function which prints an image to stdout.  Mainly meant to be
// called by PrintInputImage, PrintModelImage, and PrintWeights.

void ModelObject::PrintImage( double *pixelVector, int nColumns, int nRows )
{

  // The following fetches pixels row-by-row, starting with the bottom
  // row (i.e., what we would normally like to think of as the first row)
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++)   // step by column number = x
      printf(" %f", pixelVector[i*nColumns + j]);
    printf("\n");
  }
  printf("\n");
}


/* ---------------- PUBLIC METHOD: PrintInputImage -------------------- */
void ModelObject::PrintInputImage( )
{

  if (! dataValsSet) {
    printf("* ModelObject::PrintInputImage -- No image data supplied!\n\n");
    return;
  }
  printf("The whole input image, row by row:\n");
  PrintImage(dataVector, nDataColumns, nDataRows);
}



/* ---------------- PUBLIC METHOD: PrintModelImage -------------------- */

void ModelObject::PrintModelImage( )
{

  if (! modelImageComputed) {
    printf("* ModelObject::PrintMoelImage -- Model image has not yet been computed!\n\n");
    return;
  }
  printf("The model image, row by row:\n");
  PrintImage(modelVector, nModelColumns, nModelRows);
}


/* ---------------- PUBLIC METHOD: PrintWeights ----------------------- */

void ModelObject::PrintWeights( )
{

  if (! weightValsSet) {
    printf("* ModelObject::PrintWeights -- Weight vector has not yet been computed!\n\n");
    return;
  }
  printf("The weight image, row by row:\n");
  PrintImage(weightVector, nDataColumns, nDataRows);
}


/* ---------------- PUBLIC METHOD: PopulateParameterNames -------------- */

void ModelObject::PopulateParameterNames( )
{
  int  n;

  for (n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0,y0
      parameterLabels.push_back("X0");
      parameterLabels.push_back("Y0");
    }
    functionObjects[n]->GetParameterNames(parameterLabels);
  }
}


/* ---------------- PUBLIC METHOD: GetParameterName -------------------- */

string& ModelObject::GetParameterName( int i )
{
  if (i < nParamsTot) {
    return parameterLabels[i];
  }
  else
	return UNDEFINED;
}


/* ---------------- PUBLIC METHOD: GetNFunctions ----------------------- */

int ModelObject::GetNFunctions( )
{
  return nFunctions;
}


/* ---------------- PUBLIC METHOD: GetNParams -------------------------- */

int ModelObject::GetNParams( )
{
  return nParamsTot;
}


/* ---------------- PUBLIC METHOD: GetNDataValues ---------------------- */

int ModelObject::GetNDataValues( )
{
  return nDataVals;
}


/* ---------------- PUBLIC METHOD: GetNValidPixels --------------------- */

int ModelObject::GetNValidPixels( )
{
  return nValidDataVals;
}


/* ---------------- PUBLIC METHOD: GetModelImageVector ----------------- */

double * ModelObject::GetModelImageVector( )
{
  int  iDataRow, iDataCol, z, zModel;

  if (! modelImageComputed) {
    printf("* ModelObject::GetModelImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  outputModelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  outputModelVectorAllocated = true;
  
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels output image
    for (z = 0; z < nDataVals; z++) {
      iDataRow = z / nDataColumns;
      iDataCol = z - iDataRow*nDataColumns;
      zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
      outputModelVector[z] = modelVector[zModel];
    }
    return outputModelVector;
  }
  else
    return modelVector;
}


/* ---------------- PUBLIC METHOD: GetExpandedModelImageVector --------- */

// This differs from GetModelImageVector() in that it always returns the full
// model image, even in the case of PSF convolution (where the full model
// image will be larger than the data image!)
double * ModelObject::GetExpandedModelImageVector( )
{

  if (! modelImageComputed) {
    printf("* ModelObject::GetExpandedModelImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  return modelVector;
}


/* ---------------- PUBLIC METHOD: GetResidualImageVector -------------- */

double * ModelObject::GetResidualImageVector( )
{
  int  iDataRow, iDataCol, z, zModel;

  if (! modelImageComputed) {
    printf("* ModelObject::GetResidualImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  residualVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  residualVectorAllocated = true;
  
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in data and weight images
    for (z = 0; z < nDataVals; z++) {
      iDataRow = z / nDataColumns;
      iDataCol = z - iDataRow*nDataColumns;
      zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
      residualVector[z] = (dataVector[z] - modelVector[zModel]);
    }
  }
  else {
    // Model image is same size & shape as data and weight images
    for (z = 0; z < nDataVals; z++) {
      residualVector[z] = (dataVector[z] - modelVector[z]);
    }
  }

  return residualVector;
}


/* ---------------- PUBLIC METHOD: GetWeightImageVector ---------------- */

double * ModelObject::GetWeightImageVector( )
{
  if (! weightValsSet) {
    printf("* ModelObject::GetWeightImageVector -- Weight image has not yet been computed!\n\n");
    return NULL;
  }
  
  return weightVector;
}


/* ---------------- PUBLIC METHOD: FindTotalFluxes --------------------- */
// Estimate total fluxes for individual components (and entire model) by integrating
// over a very large image, with each component/function centered in the image.
double ModelObject::FindTotalFluxes( double params[], int xSize, int ySize,
                       double individualFluxes[] )
{
  double  x0_all, y0_all, x, y;
  double  totalModelFlux, totalComponentFlux;
  int  i, j, n;
  int  offset = 0;

  assert( (xSize >= 1) && (ySize >= 1) );
  
  // Define x0_all, y0_all as center of nominal giant image
  x0_all = 0.5*xSize;
  y0_all = 0.5*ySize;
  
  for (n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: skip over existing x0,y0 values
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0_all, y0_all);
    offset += paramSizes[n];
  }

  int  chunk = OPENMP_CHUNK_SIZE;

  totalModelFlux = 0.0;
  // Integrate over the image, once per function
  for (n = 0; n < nFunctions; n++) {
    if (functionObjects[n]->CanCalculateTotalFlux()) {
      totalComponentFlux = functionObjects[n]->TotalFlux();
    } else {
      totalComponentFlux = 0.0;
// OpenMP code currently produces wrong answers!
#pragma omp parallel private(i,j,x,y) reduction(+:totalComponentFlux)
      {
      #pragma omp for schedule (static, chunk)
      for (i = 0; i < ySize; i++) {   // step by row number = y
        y = (double)(i + 1);              // Iraf counting: first row = 1
        for (j = 0; j < xSize; j++) {   // step by column number = x
          x = (double)(j + 1);                 // Iraf counting: first column = 1
          totalComponentFlux += functionObjects[n]->GetValue(x, y);
        }
      }
    } // end omp parallel section
    } // end else [integrate total flux for component]
    individualFluxes[n] = totalComponentFlux;
    totalModelFlux += totalComponentFlux;
  }  // end for loop over functions

  return totalModelFlux;
}


/* ---------------- PROTECTED METHOD: CheckParamVector ----------------- */
// The purpose of this method is to check the parameter vector to ensure
// that all values are finite.
bool ModelObject::CheckParamVector( int nParams, double paramVector[] )
{
  bool  vectorOK = true;
  
  for (int z = 0; z < nParams; z++) {
    if (! finite(paramVector[z]))
      vectorOK = false;
  }
  
  return vectorOK;
}


/* ---------------- PROTECTED METHOD: VetDataVector -------------------- */
// The purpose of this method is to check the data vector (profile or image)
// to ensure that all non-masked pixels are finite; any non-finite pixels 
// which *are* masked will be set = 0.
bool ModelObject::VetDataVector( )
{
  bool  nonFinitePixels = false;
  bool  vectorOK = true;
  
  for (int z = 0; z < nDataVals; z++) {
    if (! finite(dataVector[z])) {
      if (weightVector[z] == 0.0)
        dataVector[z] = 0.0;
      else
        nonFinitePixels = true;
    }
  }
  
  if (nonFinitePixels) {
    printf("\n** WARNING: one or more (non-masked) pixel values in dataVector[] are non-finite!\n");
    vectorOK = false;
  }
  return vectorOK;
}


/* ---------------- PROTECTED METHOD: CheckWeightVector ---------------- */
// The purpose of this method is to check the weight vector (image) to ensure
// that all pixels are finite *and* positive.
bool ModelObject::CheckWeightVector( )
{
  bool  nonFinitePixels = false;
  bool  negativePixels = false;
  bool  weightVectorOK = true;
  
  for (int z = 0; z < nDataVals; z++) {
    if (! finite(weightVector[z]))
      nonFinitePixels = true;
    else if (weightVector[z] < 0.0)
      negativePixels = true;
  }
  
  if (nonFinitePixels) {
    printf("\n** WARNING: one or more pixel values in weightVector[] are non-finite!\n");
    weightVectorOK = false;
  }
  if (negativePixels) {
    printf("\n** WARNING: one or more pixel values in weightVector[] are < 0\n");
    weightVectorOK = false;
  }
  return weightVectorOK;
}




/* ---------------- DESTRUCTOR ----------------------------------------- */

ModelObject::~ModelObject()
{
  if (modelVectorAllocated)
    free(modelVector);
  if (weightVectorAllocated)
    free(weightVector);
  if (deviatesVectorAllocated)
    free(deviatesVector);
  if (residualVectorAllocated)
    free(residualVector);
  if (outputModelVectorAllocated)
    free(outputModelVector);
  
  if (nFunctions > 0)
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
  
  if (setStartFlag_allocated)
    free(setStartFlag);
  
  if (doConvolution)
    delete psfConvolver;

  if (bootstrapIndicesAllocated) {
    free(bootstrapIndices);
    bootstrapIndicesAllocated = false;
  }
}



/* END OF FILE: model_object.cpp --------------------------------------- */
