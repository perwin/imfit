/* FILE: model_object.cpp ---------------------------------------------- */
/* 
 * This is intended to be the main class for the "model" object (i.e., image data + 
 * fitting functions); it can also serve as a base class for derived versions thereof 
 *(e.g., fitting 1D models to profiles).
 * 
 *
 *
 * Older history:
 *     [v0.5]: 16 Apr 2010: Convolution with PSF now works, at least in principle.
 *     [v0.4]: 20--26 Mar 2010: Added stub functions to accomodate PSF image and
 * convolution.
 *     [v0.3.5]: 18--22 Feb 2010: Added PopulateParameterNames() method;
 * this generates the proper set of parameter names (including X0,Y0 for each function block)
 * and is now called by AddFunctions().  Added CheckWeightVector() method, to catch cases
 * where weight vector (whether use-supplied or calculated from Poisson statistics)
 * has "nan" values (or negative values).
 *     [v0.3]:  4 Dec 2009: Added handling of mask images.
 *     [v0.2]: 27 Nov 2009: Modified to include AddDataVectors function, which
 * will be used by derived class ModelObject1D
 *     [v0.1]: 13--15 Nov 2009: Created; initial development.
 *
 * May 2014: Now includes checks for accidental re-allocation of memory in certain
 * cases, as suggested by André Luiz de Amorim.
 */

// Copyright 2010--2015 by Peter Erwin.
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

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "mersenne_twister.h"

#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object.h"
#include "oversampled_region.h"
#include "mp_enorm.h"
#include "param_struct.h"
#include "utilities_pub.h"


/* ---------------- Definitions ---------------------------------------- */
static string  UNDEFINED = "<undefined>";

// we define these here so we don't have to worry about casting the return
// values of sizeof(), or needing to include the FFTW header. (For some reason,
// using sizeof(double) in a simple math expression *occasionally* produces
// ridiculously large values.)
const int  FFTW_SIZE = 16;
const int  DOUBLE_SIZE = 8;

// output formatting for printing parameters
#define X0_FORMAT_WITH_ERRS "%sX0\t\t%.4f # +/- %.4f\n"
#define Y0_FORMAT_WITH_ERRS "%sY0\t\t%.4f # +/- %.4f\n"
#define X0_FORMAT "%sX0\t\t%.4f\n"
#define Y0_FORMAT "%sY0\t\t%.4f\n"
#define PARAM_FORMAT_WITH_ERRS "%s%s\t\t%7g # +/- %7g\n"
#define PARAM_FORMAT "%s%s\t\t%7g\n"

// very small value for Cash statistic calculations (replaces log(m) if m <= 0)
// Based on http://cxc.harvard.edu/sherpa/ahelp/cstat.html
#define LOG_SMALL_VALUE 1.0e-25

// current best size for OpenMP processing (works well with Intel Core 2 Duo and
// Core i7 in MacBook Pro, under Mac OS X 10.6 and 10.7)
#define DEFAULT_OPENMP_CHUNK_SIZE  10


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
  maskVectorAllocated = false;
  weightVectorAllocated = false;
  standardWeightVectorAllocated = false;
  residualVectorAllocated = false;
  outputModelVectorAllocated = false;
  deviatesVectorAllocated = false;
  extraCashTermsVectorAllocated = false;
  
  doBootstrap = false;
  bootstrapIndicesAllocated = false;
  fblockStartFlags_allocated = false;
  
  // default setup = use data-based Gaussian errors + chi^2 minimization
  dataErrors = true;
  externalErrorVectorSupplied = false;
  modelErrors = false;
  useCashStatistic = false;
  poissonMLR = false;
  
  modelImageSetupDone = false;
  
  modelImageComputed = false;
  maskExists = false;
  doConvolution = false;
  oversampledRegionsExist = false;
  oversampledRegionAllocated = false;
  zeroPointSet = false;
  
  nFunctions = 0;
  nFunctionBlocks = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
  debugLevel = 0;
  verboseLevel = 0;
  
  // default image characteristics
  gain = 1.0;
  exposureTime = 1.0;
  nCombined = 1;
  effectiveGain = 1.0;
  readNoise_adu_squared = 0.0;
  
  maxRequestedThreads = 0;   // default value --> use all available processors/cores
  ompChunkSize = DEFAULT_OPENMP_CHUNK_SIZE;
  
  nDataColumns = nDataRows = 0;
  nModelColumns = nModelRows = 0;
  nPSFColumns = nPSFRows = 0;
}


/* ---------------- PUBLIC METHOD: SetDebugLevel ----------------------- */

void ModelObject::SetDebugLevel( int debuggingLevel )
{
  if (debuggingLevel < 0) {
    fprintf(stderr, "ModelObject::SetDebugLevel -- WARNING: debugging level must be > 0");
    fprintf(stderr, " (%d was supplied); debugging level left unchanged.\n", debuggingLevel);
  }
  else
    debugLevel = debuggingLevel;

  if (oversampledRegionAllocated)
    oversampledRegion->SetDebugLevel(debugLevel);
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


/* ---------------- PUBLIC METHOD: SetOMPChunkSize --------------------- */
void ModelObject::SetOMPChunkSize( int chunkSize )
{
  assert( (chunkSize >= 1) );
  ompChunkSize = chunkSize;
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



/* ---------------- PUBLIC METHOD: DefineFunctionBlocks --------------- */

void ModelObject::DefineFunctionBlocks( vector<int>& functionStartIndices )
{
  int  nn, i;
  
  nFunctionBlocks = functionStartIndices.size();
  
  // define array of [false, false, false, ...]
  // WARNING: Possible memory leak (if this function is called more than once)!
  //    If this function *is* called again, nFunctions and/or nFunctionBlocks could
  //    be different than the first call, in which we'd need to realloc fblockStartFlags
  fblockStartFlags = (bool *)calloc(nFunctions, sizeof(bool));
  fblockStartFlags_allocated = true;

  // just to be extra safe, ensure that the very first parameter is indeed start 
  // of a function-block (this will be taken care of during the loop as well)
  fblockStartFlags[0] = true;
  for (i = 0; i < nFunctionBlocks; i++) {
    nn = functionStartIndices[i];
    // function number nn is start of new function block; change fblockStartFlags[n] to true
    fblockStartFlags[nn] = true;
  }
  
  // total number of parameters = number of parameters for individual functions
  // plus x0/y0 pair for each function block
  nParamsTot = nFunctionParams + 2*nFunctionBlocks;
}



/* ---------------- PUBLIC METHOD: SetZeroPoint ----------------------- */

void ModelObject::SetZeroPoint( double zeroPointValue )
{

  zeroPoint = zeroPointValue;
  zeroPointSet = true;
}


/* ---------------- PUBLIC METHOD: AddImageDataVector ------------------ */

int ModelObject::AddImageDataVector( double *pixelVector, int nImageColumns,
                                      int nImageRows )
{
  int  status = 0;
  
  nDataVals = nValidDataVals = nImageColumns * nImageRows;
  dataVector = pixelVector;
  dataValsSet = true;
  
  status = SetupModelImage(nImageColumns, nImageRows);
  if (status < 0) {
    fprintf(stderr, "*** ERROR: AddImageDataVector: Call to SetupModelImage failed!n");
    return -1;
  }
  return 0;
}


/* ---------------- PUBLIC METHOD: SetupModelImage -------------------- */
// Called by AddImageDataVector(); can also be used by itself in make-image
// mode. Tells ModelObject to allocate space for the model image.
// Note that if PSF convolution is being done, then AddPSFVector() must be
// called *before* this method.
// nImageColumns and nImageRows should refer to the size of the data image
// (in image-fitting mode) OR the requested size of the output model image
// (in make-image mode).
int ModelObject::SetupModelImage( int nImageColumns, int nImageRows )
{
  int  result;
  assert( (nImageColumns >= 1) && (nImageRows >= 1) );
  
  nDataColumns = nImageColumns;
  nDataRows = nImageRows;
  nDataVals = nImageColumns*nImageRows;
  
  if (doConvolution) {
    nModelColumns = nDataColumns + 2*nPSFColumns;
    nModelRows = nDataRows + 2*nPSFRows;
    psfConvolver->SetupImage(nModelColumns, nModelRows);
    // NOTE: for now we're ignoring the status of psfConvolver->DoFullSetup because
    // we assume that it can't fail (we give psfConvolver the PSF info before
    // setting doConvolution to true, and we give it the image info in the line above)
    result = psfConvolver->DoFullSetup(debugLevel);
    nModelVals = nModelColumns*nModelRows;
  }
  else {
    nModelColumns = nDataColumns;
    nModelRows = nDataRows;
    nModelVals = nDataVals;
  }
  // Allocate modelimage vector
  // WARNING: Possible memory leak (if this function is called more than once)!
  //    If this function *is* called again, then nModelVals could be different
  //    from the first call, in wich case we'd need to realloc modelVector
  modelVector = (double *) calloc((size_t)nModelVals, sizeof(double));
  if (modelVector == NULL) {
    fprintf(stderr, "*** ERROR: Unable to allocate memory for model image!\n");
    fprintf(stderr, "    (Requested image size was %d x %d = %d pixels)\n", nModelRows,
    		nModelColumns, nModelVals);
    return -1;
  }
  modelVectorAllocated = true;
  modelImageSetupDone = true;
  return 0;
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

  // Avoid memory leak if pre-existing weight vector was internally allocated
  if (weightVectorAllocated) {
    free(weightVector);
    weightVectorAllocated = false;
  }
  weightVector = pixelVector;
  
  weightValsSet = true;
  externalErrorVectorSupplied = true;

  // Convert noise values into weights, if needed.  Our normal ("internal") approach is
  // to compute & store weights as 1/sigma; this assumes that whatever function calls
  // ComputeDeviates() will then square the individual (weighted) deviate values
  // in order to get the proper chi^2 result.
  // Currently, we assume three possibilities for the input weight-map pixel values:
  //    sigma (std.dev.); variance (sigma^2); and "standard" weights (1/sigma^2)
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
    case WEIGHTS_ARE_WEIGHTS:   // convert external "normal" weights to internal weights
      for (int z = 0; z < nDataVals; z++) {
        weightVector[z] = sqrt(weightVector[z]);
      }
      break;
    default:
      fprintf(stderr, "ERROR: incorrect input-type specification in ModelObject::AddErrorVector!\n");
      weightValsSet = false;
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

int ModelObject::GenerateErrorVector( )
{
  double  noise_squared, totalFlux;

  // Allocate storage for weight image:
  // WARNING: If we are calling this function for a second or subsequent time,
  // nDataVals *might* have changed; we are currently assuming it hasn't!
  if (! weightVectorAllocated) {
    weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    if (weightVector == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for weight image!\n");
      fprintf(stderr, "    (Requested image size was %d pixels)\n", nDataVals);
      return -1;
    }
    weightVectorAllocated = true;
  }
  
//  readNoise_adu_squared = readNoise*readNoise/(effectiveGain*effectiveGain);
  // Compute noise estimate for each pixel (see above for derivation)
  // Note that we assume a constant sky background (presumably already subtracted)
  for (int z = 0; z < nDataVals; z++) {
    totalFlux = dataVector[z] + originalSky;
    if (totalFlux < 0.0)
      totalFlux = 0.0;
    noise_squared = totalFlux/effectiveGain + nCombined*readNoise_adu_squared;
    // Note that we store 1/sigma instead of 1/sigma^2, since the chi^2 calculation in 
    // ChiSquared() [or the equivalent in mpfit.cpp) will square the individual terms
    weightVector[z] = 1.0 / sqrt(noise_squared);
//    printf("z = %d, noise_squared = %f\n", z, noise_squared);
  }

  weightValsSet = true;
  return 0;
}



/* ---------------- PUBLIC METHOD: GenerateExtraCashTerms -------------- */
// Generate a vector of extra terms to be added to the Cash statistic calculation
// for the case of "modified Cash statistic".
//
// This is based on treating the Cash statistic as a *likelihood ratio*, formed by
// dividing the normal Poisson likelihood by the same term evaluated for the case
// of model = data exactly. 
// The result is a likelihood function which includes an extra d_i * log d_i term for 
// each pixel i. Although we could remove them from the sum (thus getting the original
// Cash statistic), leaving it in has two advantages:
//    1. The sum (-2 log LR) has chi^2-distribution-like properties
//    2. The per-pixel values are always >= 0, and so can be used with the L-M minimizer.
//
// For the standard/original Cash statistic, all elements of this vector should be = 0
// (which is the default when we allocate the vector in UseCashStatistic).
//
// (See, e.g., Dolphin 2002, MNRAS 332: 91)
// (Suggested by David Streich, Aug 2014)
void ModelObject::GenerateExtraCashTerms( )
{
  double dataVal, extraTerm;
  
  for (int z = 0; z < nDataVals; z++) {
    dataVal = effectiveGain*(dataVector[z] + originalSky);
    // the following is strictly OK only for dataVal == 0 (lim_{x -> 0} x ln(x) = 0); 
    // the case of dataVal < 0 is undefined
    if (dataVal <= 0.0)
      extraTerm = 0.0;
    else
      extraTerm = dataVal*log(dataVal) - dataVal;
    extraCashTermsVector[z] = extraTerm;
  }
}


/* ---------------- PUBLIC METHOD: AddMaskVector ----------------------- */
// Code for adding and processing a vector containing the 2D mask image.
// Note that although our default *input* format is "0 = good pixel, > 0 =
// bad pixel", internally we convert all bad pixels to 0 and all good pixels
// to 1, so that we can multiply the weight vector by the (internal) mask values.
// The mask is applied to the weight vector by calling the ApplyMask() method
// for a given ModelObject instance.
//
// Pixels with non-finite values (e.g. NaN) are converted to "bad" (0-valued).
int ModelObject::AddMaskVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  int  returnStatus = 0;
  
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
      if (verboseLevel >= 0)
        printf("ModelObject::AddMaskVector -- treating zero-valued pixels as good ...\n");
      for (int z = 0; z < nDataVals; z++) {
        // Values of NaN or -infinity will fail > 0 test, but we want them masked, too
        if ( (! isfinite(maskVector[z])) || (maskVector[z] > 0.0) )
          maskVector[z] = 0.0;
        else {
          maskVector[z] = 1.0;
          nValidDataVals++;
        }
      }
      maskExists = true;
      break;
    case MASK_ZERO_IS_BAD:
      // Alternate form for input masks: good pixels are 1, bad pixels are 0
      if (verboseLevel >= 0)
        printf("ModelObject::AddMaskVector -- treating zero-valued pixels as bad ...\n");
      for (int z = 0; z < nDataVals; z++) {
        // Values of NaN or +infinity will fail < 1 test, but we want them masked, too
        if ( (! isfinite(maskVector[z])) || (maskVector[z] < 1.0) )
          maskVector[z] = 0.0;
        else {
          maskVector[z] = 1.0;
          nValidDataVals++;
        }
      }
      maskExists = true;
      break;
    default:
      fprintf(stderr, "ModelObject::AddMaskVector -- WARNING: unknown inputType detected!\n\n");
      returnStatus = -1;
      maskExists = false;
  }
      
  return returnStatus;
}


/* ---------------- PUBLIC METHOD: ApplyMask --------------------------- */
// Assuming both mask and weight vectors exist, this function applies the mask
// to the weight vector by multiplying the weight vector by the mask.
// (E.g., pixels with mask = 0 have weight set = 0.)
void ModelObject::ApplyMask( )
{
  double  newVal;
  
  if ( (weightValsSet) && (maskExists) ) {
    for (int z = 0; z < nDataVals; z++) {
      newVal = maskVector[z] * weightVector[z];
      // check to make sure that masked non-finite values (e.g. NaN) get zeroed
      // (because if weightVector[z] = NaN, then product will automatically be NaN)
      if ( (! isfinite(newVal)) && (maskVector[z] == 0.0) )
        newVal = 0.0;
      weightVector[z] = newVal;
    }
    if (verboseLevel >= 0) {
      printf("ModelObject: mask vector applied to weight vector. ");
      printf("(%d valid pixels remain)\n", nValidDataVals);
    }
  }
  else {
    fprintf(stderr, " ** ALERT: ModelObject::ApplyMask() called, but we are missing either\n");
    fprintf(stderr, "    error image or mask image, or both!  ApplyMask() ignored ...\n");
  }
}



/* ---------------- PUBLIC METHOD: AddPSFVector ------------------------ */
// This function is called to pass in the PSF image and dimensions; doing so
// automatically triggers setup of a Convolver object to do convolutions (including
// prep work such as computing the FFT of the PSF).
// This function must be called *before* SetupModelImage() is called (to ensure
// that we know the proper model-image dimensions), so we return an error if 
// SetupModelImage() hasn't been called yet.
int ModelObject::AddPSFVector(int nPixels_psf, int nColumns_psf, int nRows_psf,
                         double *psfPixels)
{
  int  returnStatus = 0;
  
  assert( (nPixels_psf >= 1) && (nColumns_psf >= 1) && (nRows_psf >= 1) );
  
  nPSFColumns = nColumns_psf;
  nPSFRows = nRows_psf;
  psfConvolver = new Convolver();
  psfConvolver->SetupPSF(psfPixels, nColumns_psf, nRows_psf);
  psfConvolver->SetMaxThreads(maxRequestedThreads);
  doConvolution = true;
  
  if (modelImageSetupDone) {
    fprintf(stderr, "** ERROR: PSF was added to ModelObject after SetupModelImage() was already called!\n");
    returnStatus = -1;
  }
  
  return returnStatus;
}



/* ---------------- PUBLIC METHOD: AddOversampledPSFVector ------------- */
// This function is called to pass in an oversampled PSF image, its dimensions, and
// information about the oversampling scale and the region of the image for which
// convolution of the oversampled model image with the oversampled PSF will be done.
// Calling this function automatically triggers setup of a Convolver object to handle
// oversampled convolutions (including prep work such as computing the FFT of the PSF).
// It *also* allocates space for the oversampled model sub-image vector.
//    This function *must* be called *after* SetupModelImage() [or after AddImageDataVector(),
// which itself calls SetupModelImage()], otherwise the necessary information about the 
// size of the main model image (nModelColumns, nModelRows) will not be known.
void ModelObject::AddOversampledPSFVector( int nPixels, int nColumns_psf, 
						int nRows_psf, double *psfPixels_osamp, int oversampleScale,
						int x1, int x2, int y1, int y2 )
{
  int  deltaX, deltaY, nCols_osamp, nRows_osamp;
  
  assert( (nPixels >= 1) && (nColumns_psf >= 1) && (nRows_psf >= 1) );
  assert( (oversampleScale >= 1) );
  // assertion to check that nModelColumns and nModelRows *have* been set to good values
  assert( (nModelColumns > 0) && (nModelRows > 0));

  // restrict region to be oversampled (in data or output image) to lie within image bounds
  if (x1 < 1)
    x1 = 1;
  if (y1 < 1)
    y1 = 1;
  if (x2 > nDataColumns)
    x2 = nDataColumns;
  if (y2 > nDataRows)
    y2 = nDataRows;
  // size of oversampling region
  oversamplingScale = oversampleScale;
  deltaX = x2 - x1 + 1;
  deltaY = y2 - y1 + 1;
  nCols_osamp = oversamplingScale * deltaX;
  nRows_osamp = oversamplingScale * deltaY;
  
  // oversampled PSF and corresponding Convolver object
  nPSFColumns_osamp = nColumns_psf;
  nPSFRows_osamp = nRows_psf;
  oversampledRegionsExist = true;

  // Size of actual oversampled model sub-image (including padding for PSF conv.)
  nOversampledModelColumns = nCols_osamp + 2*nPSFColumns_osamp;
  nOversampledModelRows = nRows_osamp + 2*nPSFRows_osamp;
  nOversampledModelVals = nOversampledModelColumns*nOversampledModelRows;

  // Allocate OversampledRegion object and give it necessary info
  oversampledRegion = new OversampledRegion();
  oversampledRegion->AddPSFVector(psfPixels_osamp, nPSFColumns_osamp, nPSFRows_osamp);
  oversampledRegion->SetupModelImage(x1, y1, deltaX, deltaY, nModelColumns, nModelRows, 
  									nPSFColumns, nPSFRows, oversamplingScale);
  oversampledRegionAllocated = true;
}



/* ---------------- PUBLIC METHOD: FinalSetupForFitting ---------------- */
// Call this when using ModelObject for fitting. Not necessary when just using
// ModelObject for generating model image or vector.
//    Generates blank mask vector if none already exists
//    Masks non-finite data pixels (if not already masked)
//    Generates error-based weight vector from data (if using data-based chi^2 
//       and no such  vector was supplied by user)
//    If modified Cash statistic is being used, generates extra terms from
//       data vector
//    Finally, applies mask vector to weight vector and does final vetting of
//       unmasked data values.
int ModelObject::FinalSetupForFitting( )
{
  int  nNonFinitePixels = 0;
  int  nNonFiniteErrorPixels = 0;
  int  returnStatus = 0;
  int  status = 0;
  
  // Create a default all-pixels-valid mask if no mask already exists
  if (! maskExists) {
    maskVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    if (maskVector == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for mask image!\n");
      fprintf(stderr, "    (Requested vector size was %d pixels)\n", nDataVals);
      // go ahead and return now, otherwise we'll be trying to access a null
      // pointer in the very next step
      return -1;
    }
    for (int z = 0; z < nDataVals; z++) {
      maskVector[z] = 1.0;
    }
    maskVectorAllocated = true;
    maskExists = true;
  }

  // Identify currently unmasked data pixels which have non-finite values and 
  // add those pixels to the mask
  for (int z = 0; z < nDataVals; z++) {
    if ( (maskVector[z] > 0.0) && (! isfinite(dataVector[z])) ) {
      maskVector[z] = 0.0;
      nNonFinitePixels++;
      nValidDataVals--;
    }
  }
  if ((nNonFinitePixels > 0) && (verboseLevel >= 0)) {
    if (nNonFinitePixels == 1)
      printf("ModelObject: One pixel with non-finite value found (and masked) in data image\n");
    else
      printf("ModelObject: %d pixels with non-finite values found (and masked) in data image\n", nNonFinitePixels);
  }
  
  // Generate weight vector from data-based Gaussian errors, if using chi^2 + data errors
  // and no external error map was supplied
  if ((! useCashStatistic) && (dataErrors) && (! externalErrorVectorSupplied)) {
    status = GenerateErrorVector();
    if (status < 0)  // go ahead and return now (standard behavior for memory allocation failure
      return -1;
  }
  
  // Generate extra terms vector from data for modified Cash statistic, if using latter
  if ((useCashStatistic) && (poissonMLR))
    GenerateExtraCashTerms();

  // If an external error map was supplied, identify currently unmasked *error* pixels 
  // which have non-finite values and add those pixels to the mask
  //   Possible sources of bad pixel values:
  //      1. NaN or +/-infinity in input image
  //      2. 0-valued pixels in WEIGHTS_ARE_SIGMAS case
  //      3. 0-valued or negative pixels in WEIGHTS_ARE_VARIANCES case
  //      4. Negative pixels in WEIGHTS_ARE_WEIGHTS case
  // Check only pixels which are still unmasked
  if (externalErrorVectorSupplied) {
    for (int z = 0; z < nDataVals; z++) {
      if ( (maskVector[z] > 0.0) && (! isfinite(weightVector[z])) ) {
        maskVector[z] = 0.0;
        weightVector[z] = 0.0;
        nNonFiniteErrorPixels++;
        nValidDataVals--;
      }
    }
    if ((nNonFiniteErrorPixels > 0) && (verboseLevel >= 0)) {
      if (nNonFiniteErrorPixels == 1)
        printf("ModelObject: One pixel with non-finite value found (and masked) in noise/weight image\n");
      else
        printf("ModelObject: %d pixels with non-finite values found (and masked) in noise/weight image\n", nNonFiniteErrorPixels);
    }
  }

#ifdef DEBUG
  PrintWeights();
#endif

  // Apply mask to weight vector (i.e., weight -> 0 for masked pixels)
  if (CheckWeightVector())
    ApplyMask();
  else {
    fprintf(stderr, "** ModelObject::FinalSetup -- bad values detected in weight vector!\n");
    returnStatus = -1;
  }
#ifdef DEBUG
  PrintWeights();
#endif

  if (dataValsSet) {
    bool dataOK = VetDataVector();
    if (! dataOK) {
      fprintf(stderr, "** ModelObject::FinalSetup -- bad (non-masked) data values!\n\n");
      returnStatus = -2;
    }
  }
  
#ifdef DEBUG
  PrintInputImage();
  PrintMask();
  PrintWeights();
#endif

  if (nValidDataVals < 1) {
    fprintf(stderr, "** ModelObject::FinalSetup -- not enough valid data values available for fitting!\n\n");
    returnStatus = -3;
  }

  return returnStatus;
}



/* ---------------- PUBLIC METHOD: CreateModelImage -------------------- */

void ModelObject::CreateModelImage( double params[] )
{
  double  x0, y0, x, y, newValSum;
  int  i, j, n;
  int  offset = 0;
  
  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    fprintf(stderr, "** ModelObject::CreateModelImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g, ", parameterLabels[0].c_str(), params[0]);
    for (int z = 1; z < nParamsTot; z++)
      printf(", %s = %g", parameterLabels[z].c_str(), params[z]);
    printf("\n");
#endif
  }


  // 0. Separate out the individual-component parameters and tell the associated
  // function objects to do setup work.
  // The first component's parameters start at params[0]; the second's start at
  // params[paramSizes[0]], the third at params[paramSizes[0] + paramSizes[1]], and so forth...
  for (n = 0; n < nFunctions; n++) {
    if (fblockStartFlags[n] == true) {
      // start of new function block: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0, y0);
    offset += paramSizes[n];
  }
  
  double  tempSum, adjVal, storedError;
  
  
  // 1. OK, populate modelVector with the model image -- standard pixel scaling

// Note that we cannot specify modelVector as shared [or private] bcs it is part
// of a class (not an independent variable); happily, by default all references in
// an omp-parallel section are shared unless specified otherwise
#pragma omp parallel private(i,j,n,x,y,newValSum,tempSum,adjVal,storedError)
  {
  #pragma omp for schedule (static, ompChunkSize)
//   for (i = 0; i < nModelRows; i++) {   // step by row number = y
//     y = (double)(i - nPSFRows + 1);              // Iraf counting: first row = 1 
//                                                  // (note that nPSFRows = 0 if not doing PSF convolution)
//     for (j = 0; j < nModelColumns; j++) {   // step by column number = x
//       x = (double)(j - nPSFColumns + 1);                 // Iraf counting: first column = 1
//                                                          // (note that nPSFColumns = 0 if not doing PSF convolution)
  // single-loop code which is ~ same in general case as double-loop, and
  // faster for case of small image + many cores (André Luiz de Amorim suggestion)
  for (int k = 0; k < nModelVals; k++) {
    j = k % nModelColumns;
    i = k / nModelColumns;
    y = (double)(i - nPSFRows + 1);              // Iraf counting: first row = 1
                                                 // (note that nPSFRows = 0 if not doing PSF convolution)
    x = (double)(j - nPSFColumns + 1);           // Iraf counting: first column = 1
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
  
  } // end omp parallel section
  
  
  // 2. Do PSF convolution (using standard pixel scale), if requested
  if (doConvolution)
    psfConvolver->ConvolveImage(modelVector);
  
  
  // 3. Optional generation of oversampled sub-image and convolution with oversampled PSF
  if (oversampledRegionsExist)
    oversampledRegion->ComputeRegionAndDownsample(modelVector, functionObjects, nFunctions);
  
  // [4. Possible location for charge-diffusion and other post-pixelization processing]
  
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
    fprintf(stderr, "** ModelObject::SingleFunctionImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g, ", parameterLabels[0].c_str(), params[0]);
    for (z = 1; z < nParamsTot; z++)
      printf(", %s = %g", parameterLabels[z].c_str(), params[z]);
    printf("\n");
#endif
    fprintf(stderr, "Exiting ...\n\n");
    exit(-1);
  }

  // Separate out the individual-component parameters and tell the
  // associated function objects to do setup work.
  // The first component's parameters start at params[0]; the second's
  // start at params[paramSizes[0]], the third at 
  // params[paramSizes[0] + paramSizes[1]], and so forth...
  for (n = 0; n < nFunctions; n++) {
    if (fblockStartFlags[n] == true) {
      // start of new function block: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0, y0);
    offset += paramSizes[n];
  }
  
  // OK, populate modelVector with the model image
  // OpenMP Parallel section; see CreateModelImage() for general notes on this
  // Note that since we expect this code to be called only occasionally, we have
  // not converted it to the fast-for-small-images, single-loop version used in
  // CreateModelImages()
#pragma omp parallel private(i,j,n,x,y,newVal)
  {
  #pragma omp for schedule (static, ompChunkSize)
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
      if (outputModelVector == NULL) {
        fprintf(stderr, "*** ERROR: Unable to allocate memory for output model image!\n");
        fprintf(stderr, "    (Requested image size was %d pixels)\n", nDataVals);
        return NULL;
      }
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



/* ---------------- PRIVATE METHOD: ComputePoissonMLRDeviate ----------- */
double ModelObject::ComputePoissonMLRDeviate( int i, int i_model )
{
  double   modVal, dataVal, logModel, extraTerms, deviateVal;
  
  modVal = effectiveGain*(modelVector[i_model] + originalSky);
  dataVal = effectiveGain*(dataVector[i] + originalSky);
  if (modVal <= 0)
    logModel = LOG_SMALL_VALUE;
  else
    logModel = log(modVal);
  extraTerms = extraCashTermsVector[i];
  // Note use of fabs(), to ensure that possible tiny negative values (due to
  // rounding errors when modVal =~ dataVal) don't turn into NaN
  deviateVal = sqrt(2.0 * weightVector[i] * fabs(modVal - dataVal*logModel + extraTerms));
  return deviateVal;
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
    // pixels in data and weight images (excluding the outer borders of the model image,
    // which are only for ensuring proper PSF convolution)
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        iDataRow = b / nDataColumns;
        iDataCol = b - iDataRow*nDataColumns;
        bModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        if (poissonMLR)
          yResults[z] = ComputePoissonMLRDeviate(b, bModel);
        else   // standard chi^2 term
          yResults[z] = weightVector[b] * (dataVector[b] - modelVector[bModel]);
      }
    }
    else {
      for (z = 0; z < nDataVals; z++) {
        iDataRow = z / nDataColumns;
        iDataCol = z - iDataRow*nDataColumns;
        zModel = nModelColumns*(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        if (poissonMLR)
          yResults[z] = ComputePoissonMLRDeviate(z, zModel);
        else   // standard chi^2 term
          yResults[z] = weightVector[z] * (dataVector[z] - modelVector[zModel]);
      }
    }
  }   // end if convolution case
  
  else {
    // No convolution, so model image is same size & shape as data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        if (poissonMLR)
          yResults[z] = ComputePoissonMLRDeviate(b, b);
        else   // standard chi^2 term
          yResults[z] = weightVector[b] * (dataVector[b] - modelVector[b]);
       }
    }
    else {
      for (z = 0; z < nDataVals; z++) {
        if (poissonMLR)
          yResults[z] = ComputePoissonMLRDeviate(z, z);
        else   // standard chi^2 term
          yResults[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
      }
    }
    
  }  // end else (non-convolution case)

}


/* ---------------- PUBLIC METHOD: UseModelErrors --------==----------- */

int ModelObject::UseModelErrors( )
{
  modelErrors = true;
  dataErrors = false;

  // Allocate storage for weight image (do this here because we assume that
  // AddErrorVector() will NOT be called if we're using model-based errors).
  // Set all values = 1 to start with, since we'll update this later with
  // model-based error values using UpdateWeightVector.
  
  // On the off-hand chance someone might deliberately call this after previously
  // supplying an error vector or requesting data errors (e.g., re-doing the fit
  // with only the errors changed), we allow the user to proceed even if the weight
  // vector already exists (it will be reset to all pixels = 1).
  
  // WARNING: If we are calling this function for a second or subsequent time,
  // nDataVals *might* have changed; we are currently assuming it hasn't!
  if (! weightVectorAllocated) {
    weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    if (weightVector == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for weight vector!\n");
      fprintf(stderr, "    (Requested image size was %d pixels)\n", nModelVals);
      return -1;
    }
    weightVectorAllocated = true;
  }
  else {
    fprintf(stderr, "WARNING: ModelImage::UseModelErrors -- weight vector already allocated!\n");
  }

  for (int z = 0; z < nDataVals; z++) {
    weightVector[z] = 1.0;
  }
  weightValsSet = true;
  return 0;
}


/* ---------------- PUBLIC METHOD: UseCashStatistic ------------------- */

int ModelObject::UseCashStatistic( )
{
  useCashStatistic = true;

  // On the off-hand chance someone might deliberately call this after previously
  // supplying an error vector or requesting data errors (e.g., re-doing the fit
  // with only the errors changed), we allow the user to proceed even if the weight
  // vector already exists; similarly for the extra-terms vector (the former will be 
  // reset to all pixels = 1, the latter to all pixels = 0).
  
  // WARNING: If we are calling this function for a second or subsequent time,
  // nDataVals *might* have changed; we are currently assuming it hasn't!
  if (! weightVectorAllocated) {
    weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    weightVectorAllocated = true;
  }
  else {
    fprintf(stderr, "WARNING: ModelImage::UseCashStatistic -- weight vector already allocated!\n");
  }
  if (! extraCashTermsVectorAllocated) {
    extraCashTermsVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    if (extraCashTermsVector == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for extra Cash terms vector!\n");
      fprintf(stderr, "    (Requested vector size was %d pixels)\n", nDataVals);
      return -1;
    }
    extraCashTermsVectorAllocated = true;
  }
  else {
    fprintf(stderr, "WARNING: ModelImage::UseCashStatistic -- extra-terms vector already allocated!\n");
  }

  for (int z = 0; z < nDataVals; z++) {
    weightVector[z] = 1.0;
  }
  weightValsSet = true;
  return 0;
}


/* ---------------- PUBLIC METHOD: UsePoissonMLR ----------------------- */

void ModelObject::UsePoissonMLR( )
{
  poissonMLR = true;
  UseCashStatistic();
}


/* ---------------- PUBLIC METHOD: UsingCashStatistic ------------------ */
// DEPRECATED! (Use WhichStatistic instead)

bool ModelObject::UsingCashStatistic( )
{
  return useCashStatistic;
}


/* ---------------- PUBLIC METHOD: WhichStatistic ---------------------- */

int ModelObject::WhichFitStatistic( bool verbose )
{
  if (useCashStatistic) {
    if (poissonMLR)
      return FITSTAT_POISSON_MLR;
    else
      return FITSTAT_CASH;
  }
  else
  {
    if (verbose) {
      if (modelErrors)
        return FITSTAT_CHISQUARE_MODEL;
      else if (externalErrorVectorSupplied)
        return FITSTAT_CHISQUARE_USER;
      else
        return FITSTAT_CHISQUARE_DATA;
    }
    else
      return FITSTAT_CHISQUARE;
  }
}


/* ---------------- PUBLIC METHOD: GetFitStatistic --------------------- */
/* Function for calculating chi^2 value for a model.
 *
 */
double ModelObject::GetFitStatistic( double params[] )
{
  if (useCashStatistic)
    return CashStatistic(params);  // works for both standard & modified Cash stat
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
  else {   // Model image is same size & shape as data and weight images
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
  
  // mp_enorm returns sqrt( Sum_i(chi_i^2) ) = sqrt( Sum_i(deviatesVector[i]^2) )
  if (doBootstrap)
    chi = mp_enorm(nValidDataVals, deviatesVector);
  else
    chi = mp_enorm(nDataVals, deviatesVector);
  
  return (chi*chi);
}


/* ---------------- PUBLIC METHOD: CashStatistic ----------------------- */
// Function for calculating Cash statistic for a model
//
// Note that weightVector is used here *only* for its masking purposes
//
// In the case of using Poisson MLR statistic, the extraCashTermsVector
// will be pre-populated with the appropriate terms (and will be = 0 for the
// classical Cash statistic).
//
double ModelObject::CashStatistic( double params[] )
{
  int  iDataRow, iDataCol, z, zModel, b, bModel;
  double  modVal, dataVal, logModel, extraTerms;
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
        extraTerms = extraCashTermsVector[b];   // = 0 for Cash stat
        cashStat += weightVector[b] * (modVal - dataVal*logModel + extraTerms);
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
        extraTerms = extraCashTermsVector[z];   // = 0 for Cash stat
        cashStat += weightVector[z] * (modVal - dataVal*logModel + extraTerms);
      }
    }
  }
  else {   // Model image is same size & shape as data and weight images
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        modVal = effectiveGain*(modelVector[b] + originalSky);
        dataVal = effectiveGain*(dataVector[b] + originalSky);
        if (modVal <= 0)
          logModel = LOG_SMALL_VALUE;
        else
          logModel = log(modVal);
        extraTerms = extraCashTermsVector[b];   // = 0 for Cash stat
        cashStat += weightVector[b] * (modVal - dataVal*logModel + extraTerms);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        modVal = effectiveGain*(modelVector[z] + originalSky);
        dataVal = effectiveGain*(dataVector[z] + originalSky);
        if (modVal <= 0)
          logModel = LOG_SMALL_VALUE;
        else
          logModel = log(modVal);
        extraTerms = extraCashTermsVector[z];   // = 0 for Cash stat
        cashStat += weightVector[z] * (modVal - dataVal*logModel + extraTerms);
      }
    }
  }
  
  return (2.0*cashStat);
}


/* ---------------- PUBLIC METHOD: PrintDescription ------------------- */

void ModelObject::PrintDescription( )
{
  // Don't test for verbose level, since we assume user only calls this method
  // if they *want* printed output
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
/// Basic function which prints to a file (or, e.g. stdout) a summary of the
/// best-fitting model, in form suitable for future use as an input config file.
/// If parameterInfo != NULL, then x0,y0 are corrected for any positional offsets.
///
/// If errs != NULL, then +/- errors are printed as well
///
/// If prefix != NULL, then the specified character (e.g., '#') is prepended to
/// each output line.
void ModelObject::PrintModelParams( FILE *output_ptr, double params[], 
									mp_par *parameterInfo, double errs[], 
									const char *prefix )
{
  double  x0, y0, paramVal;
  int nParamsThisFunc, k;
  int  indexOffset = 0;
  string  funcName, paramName;

  for (int n = 0; n < nFunctions; n++) {
    if (fblockStartFlags[n] == true) {
      // start of new function block: extract x0,y0 and then skip over them
      k = indexOffset;
      x0 = params[k];
      y0 = params[k + 1];
      if (parameterInfo != NULL) {
        x0 += parameterInfo[k].offset;
        y0 += parameterInfo[k + 1].offset;
      }
      if (errs != NULL) {
        fprintf(output_ptr, "\n");
        fprintf(output_ptr, X0_FORMAT_WITH_ERRS, prefix, x0, errs[k]);
        fprintf(output_ptr, Y0_FORMAT_WITH_ERRS, prefix, y0, errs[k + 1]);
      } else {
        fprintf(output_ptr, X0_FORMAT, prefix, x0);
        fprintf(output_ptr, Y0_FORMAT, prefix, y0);
      }
      indexOffset += 2;
    }
    
    // Now print the function and its parameters
    nParamsThisFunc = paramSizes[n];
    funcName = functionObjects[n]->GetShortName();
    fprintf(output_ptr, "%sFUNCTION %s\n", prefix, funcName.c_str());
    for (int i = 0; i < nParamsThisFunc; i++) {
      paramName = GetParameterName(indexOffset + i);
      paramVal = params[indexOffset + i];
      if (errs != NULL)
        fprintf(output_ptr, PARAM_FORMAT_WITH_ERRS, prefix, paramName.c_str(), paramVal, 
        		errs[indexOffset + i]);
      else
        fprintf(output_ptr, PARAM_FORMAT, prefix, paramName.c_str(), paramVal);
    }
    indexOffset += paramSizes[n];
  }
}


/* ---------------- PUBLIC METHOD: GetParamHeader ---------------------- */
/// Prints all function and parameter names in order all on one line; e.g., for use as 
/// header in bootstrap-parameters output file.
string ModelObject::GetParamHeader( )
{
  int nParamsThisFunc, nBlock;
  int  indexOffset = 0;
  string  paramName, headerLine, newString;

  headerLine = "# ";
  nBlock = 0;
  for (int n = 0; n < nFunctions; n++) {
    if (fblockStartFlags[n] == true) {
      // start of new function block: extract x0,y0 and then skip over them
      nBlock += 1;
      newString = PrintToString("X0_%d\t\tY0_%d\t\t", nBlock, nBlock);
      headerLine += newString;
      indexOffset += 2;
    }
    
    // Now print the names of the function and its parameters
    nParamsThisFunc = paramSizes[n];
    for (int i = 0; i < nParamsThisFunc; i++) {
      paramName = GetParameterName(indexOffset + i);
      newString = PrintToString("%s_%d\t", paramName.c_str(), n + 1);
      headerLine += newString;
    }
    indexOffset += paramSizes[n];
  }
  return headerLine;
}


/* ---------------- PUBLIC METHOD: UseBootstrap ------------------------ */
/// Tells ModelObject1d object that from now on we'll operate in bootstrap
/// resampling mode, so that bootstrapIndices vector is used to access the
/// data and model values (and weight values, if any).
void ModelObject::UseBootstrap( )
{
  doBootstrap = true;
  MakeBootstrapSample();
}


/* ---------------- PUBLIC METHOD: MakeBootstrapSample ----------------- */
/// Generate a new bootstrap resampling of the data (more precisely, this generate a
/// bootstrap resampling of the data *indices*)
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
/// Basic function which prints an image to stdout.  Mainly meant to be
/// called by PrintInputImage, PrintModelImage, and PrintWeights.
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
    fprintf(stderr, "* ModelObject::PrintInputImage -- No image data supplied!\n\n");
    return;
  }
  printf("The whole input image, row by row:\n");
  PrintImage(dataVector, nDataColumns, nDataRows);
}



/* ---------------- PUBLIC METHOD: PrintModelImage -------------------- */

void ModelObject::PrintModelImage( )
{

  if (! modelImageComputed) {
    fprintf(stderr, "* ModelObject::PrintMoelImage -- Model image has not yet been computed!\n\n");
    return;
  }
  printf("The model image, row by row:\n");
  PrintImage(modelVector, nModelColumns, nModelRows);
}


/* ---------------- PUBLIC METHOD: PrintMask ------------------------- */

void ModelObject::PrintMask( )
{

  if (! maskExists) {
    fprintf(stderr, "* ModelObject::PrintMask -- Mask vector does not exist!\n\n");
    return;
  }
  printf("The mask image, row by row:\n");
  PrintImage(maskVector, nDataColumns, nDataRows);
}


/* ---------------- PUBLIC METHOD: PrintWeights ----------------------- */

void ModelObject::PrintWeights( )
{

  if (! weightValsSet) {
    fprintf(stderr, "* ModelObject::PrintWeights -- Weight vector has not yet been computed!\n\n");
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
    if (fblockStartFlags[n] == true) {
      // start of new function block: extract x0,y0
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
    fprintf(stderr, "* ModelObject::GetModelImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  if (! outputModelVectorAllocated) {
    outputModelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    outputModelVectorAllocated = true;
  }
  
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

/// This differs from GetModelImageVector() in that it always returns the full
/// model image, even in the case of PSF convolution (where the full model
/// image will be larger than the data image!)
double * ModelObject::GetExpandedModelImageVector( )
{

  if (! modelImageComputed) {
    fprintf(stderr, "* ModelObject::GetExpandedModelImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  return modelVector;
}


/* ---------------- PUBLIC METHOD: GetResidualImageVector -------------- */

double * ModelObject::GetResidualImageVector( )
{
  int  iDataRow, iDataCol, z, zModel;

  if (! modelImageComputed) {
    fprintf(stderr, "* ModelObject::GetResidualImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  // WARNING: If we are calling this function for a second or subsequent time,
  // nDataVals *might* have changed; we are currently assuming it hasn't!
  if (! residualVectorAllocated) {
    residualVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    residualVectorAllocated = true;
  }
    
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
/// Returns the weightVector converted to 1/sigma^2 (i.e., 1/variance) form
double * ModelObject::GetWeightImageVector( )
{
  if (! weightValsSet) {
    fprintf(stderr, "* ModelObject::GetWeightImageVector -- Weight image has not yet been computed!\n\n");
    return NULL;
  }
  
  standardWeightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  for (int z = 0; z < nDataVals; z++) {
    double  w_sqrt = weightVector[z];   // internal weight value (sqrt of formal weight)
  	standardWeightVector[z] = w_sqrt*w_sqrt;
  }
  standardWeightVectorAllocated = true;
  return standardWeightVector;
}


/* ---------------- PUBLIC METHOD: GetDataVector ----------------------- */

double * ModelObject::GetDataVector( )
{
  if (! dataValsSet) {
    fprintf(stderr, "* ModelObject::GetDataVector -- Image data values have not yet been supplied!\n\n");
    return NULL;
  }
  
  return dataVector;
}


/* ---------------- PUBLIC METHOD: FindTotalFluxes --------------------- */
/// Estimate total fluxes for individual components (and entire model) by integrating
/// over a very large image, with each component/function centered in the image.
/// Total flux is returned by the function; fluxes for individual components are
/// returned in individualFluxes.
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
    if (fblockStartFlags[n] == true) {
      // start of new function block: skip over existing x0,y0 values
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0_all, y0_all);
    offset += paramSizes[n];
  }

//  int  chunk = OPENMP_CHUNK_SIZE;

  totalModelFlux = 0.0;
  // Integrate over the image, once per function
  for (n = 0; n < nFunctions; n++) {
    if (functionObjects[n]->CanCalculateTotalFlux()) {
      totalComponentFlux = functionObjects[n]->TotalFlux();
    } else {
      totalComponentFlux = 0.0;
// Note: since this bit of OpenMP code explicitly involves integrating over a very large
// image, we don't bother using the fast-for-small-images, single-loop version that's 
// used in CreateModelImage()
#pragma omp parallel private(i,j,x,y) reduction(+:totalComponentFlux)
      {
      #pragma omp for schedule (static, ompChunkSize)
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
/// Returns true if all values in the parameter vector are finite.
bool ModelObject::CheckParamVector( int nParams, double paramVector[] )
{
  bool  vectorOK = true;
  
  for (int z = 0; z < nParams; z++) {
    if (! isfinite(paramVector[z]))
      vectorOK = false;
  }
  
  return vectorOK;
}


/* ---------------- PROTECTED METHOD: VetDataVector -------------------- */
/// Returns true if all non-masked pixels in the image data vector are finite;
/// returns false if one or more are not, and prints an error message to stderr.
/// ALSO sets any masked pixels which are non-finite to 0.
///
/// More simply: the purpose of this method is to check the data vector (profile or image)
/// to ensure that all non-masked pixels are finite.
/// Any non-finite pixels which *are* masked will be set = 0.
bool ModelObject::VetDataVector( )
{
  bool  nonFinitePixels = false;
  bool  vectorOK = true;
  
  for (int z = 0; z < nDataVals; z++) {
    if (! isfinite(dataVector[z])) {
      if (maskVector[z] > 0.0)
        nonFinitePixels = true;
      else
        dataVector[z] = 0.0;
    }
  }
  
  if (nonFinitePixels) {
    fprintf(stderr, "\n** WARNING: one or more (non-masked) pixel values in data image are non-finite!\n");
    vectorOK = false;
  }
  return vectorOK;
}


/* ---------------- PROTECTED METHOD: CheckWeightVector ---------------- */
/// Returns true if all pixels in the weight vector are finite *and* nonnegative.
bool ModelObject::CheckWeightVector( )
{
  bool  nonFinitePixels = false;
  bool  negativePixels = false;
  bool  weightVectorOK = true;
  
  // check individual pixels in weightVector, but only if they aren't masked by maskVector
  if (maskExists) {
    for (int z = 0; z < nDataVals; z++) {
      if (maskVector[z] > 0.0) {
        if (! isfinite(weightVector[z]))
          nonFinitePixels = true;
        else if (weightVector[z] < 0.0)
          negativePixels = true;
      }
    }  
  }
  else {
    for (int z = 0; z < nDataVals; z++) {
      if (! isfinite(weightVector[z]))
        nonFinitePixels = true;
    }
  }
  
  if (nonFinitePixels) {
    fprintf(stderr, "\n** WARNING: one or more pixel values in weightVector[] are non-finite!\n");
    if (externalErrorVectorSupplied)
      fprintf(stderr, "     (Bad values in external noise or weight image)\n");
    else
      fprintf(stderr, "     (Negative pixel values in data image -- missing sky background?)\n");
    weightVectorOK = false;
  }
  if (negativePixels) {
    fprintf(stderr, "\n** WARNING: one or more pixel values in weightVector[] are < 0\n");
    fprintf(stderr, "     (Negative pixel values in noise or weight image?)\n");
    if (originalSky <= 0)
    	fprintf(stderr, "     (original-sky-background = %f -- missing or wrong value?\n", originalSky);
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
  if (standardWeightVectorAllocated)
    free(standardWeightVector);
  if (maskVectorAllocated)   // only true if we construct mask vector internally
    free(maskVector);
  if (deviatesVectorAllocated)
    free(deviatesVector);
  if (residualVectorAllocated)
    free(residualVector);
  if (outputModelVectorAllocated)
    free(outputModelVector);
  if (extraCashTermsVectorAllocated)
    free(extraCashTermsVector);
  
  if (nFunctions > 0)
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
  
  if (fblockStartFlags_allocated)
    free(fblockStartFlags);
  
  if (doConvolution)
    delete psfConvolver;
  if (oversampledRegionsExist)
    delete oversampledRegion;

  if (bootstrapIndicesAllocated) {
    free(bootstrapIndices);
    bootstrapIndicesAllocated = false;
  }
}



/* END OF FILE: model_object.cpp --------------------------------------- */
