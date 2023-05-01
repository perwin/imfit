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

// Copyright 2010--2020 by Peter Erwin.
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
//#include <math.h>
#include <cmath>
#include <iostream>
#include <tuple>

using namespace std;

#include "mersenne_twister.h"
#include "definitions.h"
#include "function_object.h"
#include "model_object.h"
#include "oversampled_region.h"
#include "psf_oversampling_info.h"
#include "psf_interpolators.h"
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

// very small value for Cash statistic calculations (replaces log(m) if m <= 0)
// Based on http://cxc.harvard.edu/sherpa/ahelp/cstat.html
#define LOG_SMALL_VALUE 1.0e-25

// current best size for OpenMP processing (works well with Intel Core 2 Duo and
// Core i7 in MacBook Pro, under Mac OS X 10.6 and 10.7)
#define DEFAULT_OPENMP_CHUNK_SIZE  10


// for use in ModelObject::AddFunction()
map<string, int> interpolationMap{ {string("bicubic"), kInterpolator_bicubic}, 
									{string("lanczos2"), kInterpolator_lanczos2}, 
									{string("lanczos3"), kInterpolator_lanczos3} };



/* ------------------- Function Prototypes ----------------------------- */

void NormalizePSF( double *psfPixels, long nPixels_psf );



// NOTE: (parts of) the following class are used in PyImfit

// ModelObject class

/* ---------------- CONSTRUCTOR ---------------------------------------- */
/// Constructor
ModelObject::ModelObject( )
{
  dataValsSet = weightValsSet = false;

  dataVector = modelVector = weightVector = standardWeightVector = NULL;
  residualVector = maskVector = deviatesVector = NULL;
  outputModelVector = extraCashTermsVector = NULL;
  bootstrapIndices = NULL;
  fsetStartFlags = NULL;

  localPsfPixels = nullptr;
  psfInterpolator = nullptr;
  psfInterpolator_allocated = false;
  
  modelVectorAllocated = false;
  maskVectorAllocated = false;
  weightVectorAllocated = false;
  standardWeightVectorAllocated = false;
  residualVectorAllocated = false;
  outputModelVectorAllocated = false;
  deviatesVectorAllocated = false;
  extraCashTermsVectorAllocated = false;
  localPsfPixels_allocated = false;
  
  fsetStartFlags_allocated = false;
  
  // default setup = use data-based Gaussian errors + chi^2 minimization
  dataErrors = true;
  externalErrorVectorSupplied = false;
  modelErrors = false;
  useCashStatistic = false;
  poissonMLR = false;
  doBootstrap = false;
  bootstrapIndicesAllocated = false;

  modelImageSetupDone = false;
  
  modelImageComputed = false;
  maskExists = false;
  weightValsSet = false;
  doConvolution = false;
  oversampledRegionsExist = false;
  zeroPointSet = false;
  pointSourcesPresent = false;
  
  nFunctions = 0;
  nFunctionSets = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
  nOversampledRegions = 0;
  debugLevel = 0;
  verboseLevel = 0;
  
  // default image characteristics
  gain = 1.0;
  exposureTime = effectiveGain = 1.0;
  nCombined = 1;
  originalSky = readNoise = readNoise_adu_squared = 0.0;
  
  imageOffset_X0 = imageOffset_Y0 = 0;
  
  maxRequestedThreads = 0;   // default value --> use all available processors/cores
  ompChunkSize = DEFAULT_OPENMP_CHUNK_SIZE;
  
  nDataVals = nDataColumns = nDataRows = 0;
  nModelVals = nModelColumns = nModelRows = 0;
  nPSFColumns = nPSFRows = 0;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */
/// Destructor
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
  if (localPsfPixels_allocated)
    free(localPsfPixels);

  if (psfInterpolator_allocated)
    delete psfInterpolator;
  
  if (nFunctions > 0)
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
  
  if (fsetStartFlags_allocated)
    free(fsetStartFlags);
  
  if (doConvolution)
    delete psfConvolver;
  if (oversampledRegionsExist) {
    // since these were originally created with "new", we have to deallocate with "delete"
    for (int i = 0; i < nOversampledRegions; i++)
      delete oversampledRegionsVect[i];
    oversampledRegionsVect.clear();
    nOversampledRegions = 0;
    oversampledRegionsExist = false;
  }
  
  if (bootstrapIndicesAllocated) {
    free(bootstrapIndices);
    bootstrapIndicesAllocated = false;
  }
}



/* ---------------- PUBLIC METHOD: SetVerboseLevel ---------------------- */
/// Set the verbosity level
void ModelObject::SetVerboseLevel( int verbosity )
{
  verboseLevel = verbosity;
}


/* ---------------- PUBLIC METHOD: SetDebugLevel ----------------------- */
/// Set the debugging level (must be 0 [default] or larger).
void ModelObject::SetDebugLevel( int debuggingLevel )
{
  if (debuggingLevel < 0) {
    fprintf(stderr, "ModelObject::SetDebugLevel -- WARNING: debugging level must be > 0");
    fprintf(stderr, " (%d was supplied); debugging level left unchanged.\n", debuggingLevel);
  }
  else
    debugLevel = debuggingLevel;
}


/* ---------------- PUBLIC METHOD: SetMaxThreads ----------------------- */
/// Specify the maximum number of OpenMP threads to use in computations;
/// also sets maximum number FFTW threads for convolutions.
void ModelObject::SetMaxThreads( int maxThreadNumber )
{
  assert( (maxThreadNumber >= 1) );
  maxRequestedThreads = maxThreadNumber;
#ifdef USE_OPENMP
  omp_set_num_threads(maxRequestedThreads);
#endif
}


/* ---------------- PUBLIC METHOD: SetOMPChunkSize --------------------- */
/// Sets the chunk size for OpenMP
void ModelObject::SetOMPChunkSize( int chunkSize )
{
  assert( (chunkSize >= 1) );
  ompChunkSize = chunkSize;
}


/* ---------------- PUBLIC METHOD: AddFunction ------------------------- */
/// Adds a FunctionObject subclass to the model
int ModelObject::AddFunction( FunctionObject *newFunctionObj_ptr )
{
  int  nNewParams, result;
  
  functionObjects.push_back(newFunctionObj_ptr);
  nFunctions += 1;
  nNewParams = newFunctionObj_ptr->GetNParams();
  paramSizes.push_back(nNewParams);
  nFunctionParams += nNewParams;
  
  // handle optional case of PointSource function
  if (newFunctionObj_ptr->IsPointSource()) {
    if (! psfInterpolator_allocated) {
      string interpName = newFunctionObj_ptr->GetInterpolationType();
      int interpType = interpolationMap[interpName];
      result = SetupPsfInterpolation(interpType);
      if (result < 0)
      	return -1;
    }
    newFunctionObj_ptr->AddPsfInterpolator(psfInterpolator);
    pointSourcesPresent = true;
  }
  
  return 0;
}


/* ---------------- PUBLIC METHOD: SetupPsfInterpolation -------------- */
/// Specify that PSF interpolation (by PointSource functions) will be used;
/// causes an internal PsfInterpolator object of the appropriate subclass
/// to be allocated and set up with previously supplied PSF data.
/// This should only be called *after* AddPSFVector has been called to supply the
/// PSF data.
int ModelObject::SetupPsfInterpolation( int interpolationType )
{

  if (! localPsfPixels_allocated) {
    fprintf(stderr, "** ERROR: PointSource image function being used, but no ");
    fprintf(stderr, "PSF image was supplied!\n");
    return -1;
  }
  
  // the following is definitely needed for bicubic interpolation; not sure what
  // the limits for lanczos2 or lanczos3 should be
  if ((nPSFColumns < 4) || (nPSFRows < 4)) {
    fprintf(stderr, "** ERROR: PSF image is too small for interpolation with PointSource functions!\n");
    fprintf(stderr, "   (must be at least 4 x 4 pixels in size for GSL bicubic interpolation)\n");
    return -2;
  }

  switch (interpolationType) {
    case kInterpolator_bicubic:
      psfInterpolator = new PsfInterpolator_bicubic(localPsfPixels, nPSFColumns, nPSFRows);
	  psfInterpolator_allocated = true;
      break;
    case kInterpolator_lanczos2:
      psfInterpolator = new PsfInterpolator_lanczos2(localPsfPixels, nPSFColumns, nPSFRows);
      psfInterpolator_allocated = true;
      break;
    case kInterpolator_lanczos3:
//       printf("ERROR: Lanczos3 interpolation not yet implemented!\n");
//       return -2;
      psfInterpolator = new PsfInterpolator_lanczos3(localPsfPixels, nPSFColumns, nPSFRows);
      psfInterpolator_allocated = true;
      break;
  }
  
  return 0;
}


/* ---------------- PUBLIC METHOD: DefineFunctionSets ---------------- */
/// Given a vector of indices specifying the start of individual function
/// sets, set up the internal fsetStartFlags array; also compute nParamsTot
void ModelObject::DefineFunctionSets( vector<int>& functionStartIndices )
{
  int  nn, i;
  
  nFunctionSets = functionStartIndices.size();
  
  // define array of [false, false, false, ...]
  // WARNING: Possible memory leak (if this function is called more than once)!
  //    If this function *is* called again, nFunctions and/or nFunctionSets could
  //    be different than the first call, in which we'd need to realloc fsetStartFlags
  fsetStartFlags = (bool *)calloc(nFunctions, sizeof(bool));
  fsetStartFlags_allocated = true;

  // just to be extra safe, ensure that the very first parameter is indeed start 
  // of a function set (this will be taken care of during the loop as well)
  fsetStartFlags[0] = true;
  for (i = 0; i < nFunctionSets; i++) {
    nn = functionStartIndices[i];
    // function number nn is start of new function set; change fsetStartFlags[n] to true
    fsetStartFlags[nn] = true;
  }
  
  // total number of parameters = number of parameters for individual functions
  // plus x0/y0 pair for each function set
  // (add computed number to pre-existing value of nParamsTot [= 0 normally]
  // for special case of multimfit, etc., where additional image-description
  // parameters may have been counted)
  nParamsTot += nFunctionParams + 2*nFunctionSets;
}



/* ---------------- PUBLIC METHOD: SetZeroPoint ----------------------- */
/// Set the internal zero-point (for printing flux values)
void ModelObject::SetZeroPoint( double zeroPointValue )
{

  zeroPoint = zeroPointValue;
  zeroPointSet = true;
}


/* ---------------- PUBLIC METHOD: PrintParameterInfoLine ------------- */

void ModelObject::AddParameterInfo( mp_par  *inputParameterInfo )
{
  SimpleParameterInfo paramStruct;

  for (int i = 0; i < nParamsTot; i++) {
    paramStruct.fixed = inputParameterInfo[i].fixed;
    paramStruct.limited[0] = inputParameterInfo[i].limited[0];
    paramStruct.limited[1] = inputParameterInfo[i].limited[1];
    paramStruct.limits[0] = inputParameterInfo[i].limits[0];
    paramStruct.limits[1] = inputParameterInfo[i].limits[1];
    parameterInfoVect.push_back(paramStruct);
  }
}

void ModelObject::AddParameterInfo( vector<mp_par> inputParameterInfo )
{
  SimpleParameterInfo paramStruct;

  for (int i = 0; i < nParamsTot; i++) {
    paramStruct.fixed = inputParameterInfo[i].fixed;
    paramStruct.limited[0] = inputParameterInfo[i].limited[0];
    paramStruct.limited[1] = inputParameterInfo[i].limited[1];
    paramStruct.limits[0] = inputParameterInfo[i].limits[0];
    paramStruct.limits[1] = inputParameterInfo[i].limits[1];
    parameterInfoVect.push_back(paramStruct);
  }
}


/* ---------------- PUBLIC METHOD: AddImageDataVector ------------------ */

int ModelObject::AddImageDataVector( double *pixelVector, int nImageColumns,
                                      int nImageRows )
{
  int  status = 0;
  
  nDataVals = nValidDataVals = (long)nImageColumns * (long)nImageRows;
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
  int  result = 0;
  assert( (nImageColumns >= 1) && (nImageRows >= 1) );
  
  nDataColumns = nImageColumns;
  nDataRows = nImageRows;
  nDataVals = (long)nImageColumns * (long)nImageRows;
  
  if (doConvolution) {
    nModelColumns = nDataColumns + 2*nPSFColumns;
    nModelRows = nDataRows + 2*nPSFRows;
    psfConvolver->SetupImage(nModelColumns, nModelRows);
    result = psfConvolver->DoFullSetup(debugLevel);
    if (result < 0) {
      fprintf(stderr, "*** Error returned from Convolver::DoFullSetup!\n");
      return result;
    }
    nModelVals = (long)nModelColumns * (long)nModelRows;
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
    fprintf(stderr, "    (Requested image size was %d x %d = %ld pixels)\n", nModelRows,
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


// Image offsets are defined as offset_X0 = x_ll - 1, y_ll - offset_Y0,
// where the IRAF-coordinate of the image-subsection's lower-left corner
// is (x_ll,y_ll). When the entire image is being used, (x_ll,y_ll) = (1,1),
// and offset_X0 = offset_Y0 = 0.
/* ---------------- PUBLIC METHOD: AddImageOffsets --------------------- */

void ModelObject::AddImageOffsets( int offset_X0, int offset_Y0 )
{
  imageOffset_X0 = offset_X0;
  imageOffset_Y0 = offset_Y0;
}


/* ---------------- PUBLIC METHOD: GetImageOffsets --------------------- */

std::tuple<int, int> ModelObject::GetImageOffsets( )
{
  return std::make_tuple(imageOffset_X0, imageOffset_Y0);
}


/* ---------------- PUBLIC METHOD: AddErrorVector ---------------------- */

void ModelObject::AddErrorVector( long nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  assert( (nDataValues == nDataVals) && (nImageColumns == nDataColumns) && 
          (nImageRows == nDataRows) );

  long  z;
  
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
      for (z = 0; z < nDataVals; z++) {
        weightVector[z] = 1.0 / weightVector[z];
      }
      break;
    case WEIGHTS_ARE_VARIANCES:
      for (z = 0; z < nDataVals; z++) {
        weightVector[z] = 1.0 / sqrt(weightVector[z]);
      }
      break;
    case WEIGHTS_ARE_WEIGHTS:   // convert external "normal" weights to internal weights
      for (z = 0; z < nDataVals; z++) {
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
      fprintf(stderr, "    (Requested image size was %ld pixels)\n", nDataVals);
      return -1;
    }
    weightVectorAllocated = true;
  }
  
  // Compute noise estimate for each pixel (see above for derivation)
  // Note that we assume a constant sky background (presumably already subtracted)
  for (long z = 0; z < nDataVals; z++) {
    totalFlux = dataVector[z] + originalSky;
    if (totalFlux < 0.0)
      totalFlux = 0.0;
    noise_squared = totalFlux/effectiveGain + nCombined*readNoise_adu_squared;
    // Note that we store 1/sigma instead of 1/sigma^2, since the chi^2 calculation in 
    // ChiSquared() [or the equivalent in mpfit.cpp) will square the individual terms
    weightVector[z] = 1.0 / sqrt(noise_squared);
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
  
  for (long z = 0; z < nDataVals; z++) {
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
int ModelObject::AddMaskVector( long nDataValues, int nImageColumns,
                                int nImageRows, double *pixelVector, int inputType )
{
  int  returnStatus = 0;
  bool  nonFiniteValues = false;
  long  z;
  
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
      for (z = 0; z < nDataVals; z++) {
        // Values of NaN or -infinity will fail > 0 test, but we want them masked, too
        if ( (! isfinite(maskVector[z])) || (maskVector[z] > 0.0) ) {
          if (! isfinite(maskVector[z]))
            nonFiniteValues = true;
          maskVector[z] = 0.0;
        }
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
      for (z = 0; z < nDataVals; z++) {
        // Values of NaN or +infinity will fail < 1 test, but we want them masked, too
        if ( (! isfinite(maskVector[z])) || (maskVector[z] < 1.0) ) {
          if (! isfinite(maskVector[z]))
            nonFiniteValues = true;
          maskVector[z] = 0.0;
        }
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

  if (nonFiniteValues) {
    fprintf(stderr, " ** WARNING: ModelObject::AddMaskVector -- one or more non-finite values in mask image!\n");
    fprintf(stderr, "    (Will be treated as masked)\n");
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
    for (long z = 0; z < nDataVals; z++) {
      newVal = maskVector[z] * weightVector[z];
      // check to make sure that masked non-finite values (e.g. NaN) get zeroed
      // (because if weightVector[z] = NaN, then product will automatically be NaN)
      if ( (! isfinite(newVal)) && (maskVector[z] == 0.0) )
        newVal = 0.0;
      weightVector[z] = newVal;
    }
    if (verboseLevel >= 0) {
      printf("ModelObject: mask vector applied to weight vector. ");
      printf("(%ld valid pixels remain)\n", nValidDataVals);
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
// SetupModelImage() has previously been called.
int ModelObject::AddPSFVector( long nPixels_psf, int nColumns_psf, int nRows_psf,
                         	double *psfPixels, bool normalizePSF )
{
  double  pixVal;
  double  tempSum = 0.0;
  int  returnStatus = 0;
  
  assert( (nPixels_psf >= 1) && (nColumns_psf >= 1) && (nRows_psf >= 1) );

  // store PSF pixels locally (and check for bad pixel values)
  localPsfPixels = (double *) calloc((size_t)nPixels_psf, sizeof(double));
  localPsfPixels_allocated = true;
  for (long i = 0; i < nPixels_psf; i++) {
  	pixVal = psfPixels[i];
    if (! isfinite(pixVal)) {
      fprintf(stderr, "** ERROR: PSF image has one or more non-finite values!\n");
      free(localPsfPixels);
      localPsfPixels_allocated = false;
      return -1;
    }
    localPsfPixels[i] = pixVal;
    tempSum += pixVal;
  }
  if (normalizePSF) {
    if (tempSum <= 0.0) {
      fprintf(stderr, "** ERROR: Sum of PSF pixel values is <= 0 -- PSF cannot be normalized!\n");
      free(localPsfPixels);
      localPsfPixels_allocated = false;
      return -1;
    }
    NormalizePSF(localPsfPixels, nPixels_psf);
  }

  // Finally, set up Convolver object
  nPSFColumns = nColumns_psf;
  nPSFRows = nRows_psf;
  psfConvolver = new Convolver();
  psfConvolver->SetupPSF(psfPixels, nColumns_psf, nRows_psf, normalizePSF);
//   printf("ModelObject::AddPSFVector -- calling psfConvolver->SetMaxThreads()");
//   printf("with maxRequestedThreads = %d\n", maxRequestedThreads);
  psfConvolver->SetMaxThreads(maxRequestedThreads);
  doConvolution = true;
  
  if (modelImageSetupDone) {
    fprintf(stderr, "** ERROR: PSF was added to ModelObject after SetupModelImage() was already called!\n");
    returnStatus = -1;
  }
  
  return returnStatus;
}



/* ---------------- PUBLIC METHOD: AddOversampledPsfInfo --------------- */
int ModelObject::AddOversampledPsfInfo( PsfOversamplingInfo *oversampledPsfInfo )
{
  long  nPixels;
  int  nColumns_psf, nRows_psf, oversampleScale;
  int  x1, x2, y1, y2;
  double  *psfPixels_osamp;
  int  deltaX, deltaY, nCols_osamp, nRows_osamp;
  int  status = 0;
  
  nColumns_psf = oversampledPsfInfo->GetNColumns();
  nRows_psf = oversampledPsfInfo->GetNRows();
  nPixels = (long)nColumns_psf * (long)nRows_psf;
  oversampleScale = oversampledPsfInfo->GetOversamplingScale();
  std::tie(x1, x2, y1, y2) = oversampledPsfInfo->GetCorrectedRegionCoords();
  psfPixels_osamp = oversampledPsfInfo->GetPsfPixels();
  
  assert( (nPixels >= 1) && (nColumns_psf >= 1) && (nRows_psf >= 1) );
  assert( (x1 < x2) && (y1 < y2) );
  assert( (oversampleScale >= 1) );
  // assertion to check that nModelColumns and nModelRows *have* been set to good values
  assert( (nModelColumns > 0) && (nModelRows > 0) );

  for (long i = 0; i < nPixels; i++) {
    if (! isfinite(psfPixels_osamp[i])) {
      fprintf(stderr, "** ERROR: Oversampled PSF image has one or more non-finite values!\n");
      return -1;
    }
  }

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
  deltaX = x2 - x1 + 1;
  deltaY = y2 - y1 + 1;
  assert( deltaX > 0 );
  assert( deltaY > 0 );
  nCols_osamp = oversampleScale * deltaX;
  nRows_osamp = oversampleScale * deltaY;
  
  // oversampled PSF and corresponding Convolver object
  nPSFColumns_osamp = nColumns_psf;
  nPSFRows_osamp = nRows_psf;
  oversampledRegionsExist = true;

  // Size of actual oversampled model sub-image (including padding for PSF conv.)
  nOversampledModelColumns = nCols_osamp + 2*nPSFColumns_osamp;
  nOversampledModelRows = nRows_osamp + 2*nPSFRows_osamp;
  nOversampledModelVals = (long)nOversampledModelColumns * (long)nOversampledModelRows;

  // Allocate OversampledRegion object and give it necessary info
  OversampledRegion *oversampledRegion = new OversampledRegion();
  oversampledRegion->SetDebugLevel(debugLevel);
  oversampledRegion->AddPSFVector(psfPixels_osamp, nPSFColumns_osamp, nPSFRows_osamp,
  									oversampledPsfInfo->GetNormalizationFlag());
  status = oversampledRegion->SetupModelImage(x1, y1, deltaX, deltaY, nModelColumns, nModelRows, 
  									nPSFColumns, nPSFRows, oversampleScale);
  if (status < 0) {
    fprintf(stderr, "*** ERROR: AddOversampledPsfInfo: Call to oversampledRegion->SetupModelImage failed!n");
    return -1;
  }
  oversampledRegionsVect.push_back(oversampledRegion);
  nOversampledRegions++;
  
  return 0;
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
  long  nNonFinitePixels = 0;
  long  nNonFiniteErrorPixels = 0;
  long  z;
  int  returnStatus = 0;
  int  status = 0;
  
  // Create a default all-pixels-valid mask if no mask already exists
  if (! maskExists) {
    maskVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    if (maskVector == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for mask image!\n");
      fprintf(stderr, "    (Requested vector size was %ld pixels)\n", nDataVals);
      // go ahead and return now, otherwise we'll be trying to access a null
      // pointer in the very next step
      return -1;
    }
    for (z = 0; z < nDataVals; z++) {
      maskVector[z] = 1.0;
    }
    maskVectorAllocated = true;
    maskExists = true;
  }

  // Identify currently unmasked data pixels which have non-finite values and 
  // add those pixels to the mask
  for (z = 0; z < nDataVals; z++) {
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
      printf("ModelObject: %ld pixels with non-finite values found (and masked) in data image\n", nNonFinitePixels);
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
    for (z = 0; z < nDataVals; z++) {
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
        printf("ModelObject: %ld pixels with non-finite values found (and masked) in noise/weight image\n", nNonFiniteErrorPixels);
    }
  }

#ifdef DEBUG
  PrintWeights();
#endif

  // Apply mask to weight vector (i.e., set weight -> 0 for masked pixels)
  ApplyMask();
  if (! CheckWeightVector()) {
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
  long  i, j;
  int  n;
  int  offset = 0;
  
  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    fprintf(stderr, "** ModelObject::CreateModelImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g, ", parameterLabels[0].c_str(), params[0]);
    for (int np = 1; np < nParamsTot; np++)
      printf(", %s = %g", parameterLabels[np].c_str(), params[np]);
    printf("\n");
#endif
  }


  // 0. Separate out the individual-component parameters and tell the associated
  // function objects to do setup work.
  // The first component's parameters start at params[0]; the second's start at
  // params[paramSizes[0]], the third at params[paramSizes[0] + paramSizes[1]], and so forth...
  for (n = 0; n < nFunctions; n++) {
    if (fsetStartFlags[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0, y0);
    offset += paramSizes[n];
  }
  
  
  // 1. OK, populate modelVector with the model image -- standard pixel scaling
  double  tempSum, adjVal, storedError;
  
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
  for (long k = 0; k < nModelVals; k++) {
    j = k % nModelColumns;
    i = k / nModelColumns;
    y = (double)(i - nPSFRows + 1);              // Iraf counting: first row = 1
                                                 // (note that nPSFRows = 0 if not doing PSF convolution)
    x = (double)(j - nPSFColumns + 1);           // Iraf counting: first column = 1
                                                 // (note that nPSFColumns = 0 if not doing PSF convolution)
    newValSum = 0.0;
    storedError = 0.0;
    for (n = 0; n < nFunctions; n++) {
      if (! functionObjects[n]->IsPointSource()) {
        // Kahan summation algorithm
        adjVal = functionObjects[n]->GetValue(x, y) - storedError;
        tempSum = newValSum + adjVal;
        storedError = (tempSum - newValSum) - adjVal;
        newValSum = tempSum;
      }
    }
    modelVector[i*nModelColumns + j] = newValSum;
  }
  } // end omp parallel section
  
  
  // 2. Do PSF convolution (using standard pixel scale), if requested
  if (doConvolution)
    psfConvolver->ConvolveImage(modelVector);
  
  
  // 2.B Add flux from PointSource functions, if present 
  // (must be done *after* PSF convolution!)
  if (pointSourcesPresent) {
    // Re-assign psfInterpolator object (bcs. calls made to
    // OversampledRegion::ComputeRegionAndDownsample result in PointSource objects 
    // getting assigned alternate psfInterpolators), so we have to reset PointSource 
    // objects to use the standard-resolution psfInterpolator object held by ModelObject
    for (FunctionObject *funcObj : functionObjects)
      if (funcObj->IsPointSource())
        funcObj->AddPsfInterpolator(psfInterpolator);
    
#pragma omp parallel private(i,j,n,x,y,newValSum,tempSum,adjVal,storedError)
    {
    #pragma omp for schedule (static, ompChunkSize)
    for (long k = 0; k < nModelVals; k++) {
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
        if (functionObjects[n]->IsPointSource()) {
          adjVal = functionObjects[n]->GetValue(x, y) - storedError;
          tempSum = newValSum + adjVal;
          storedError = (tempSum - newValSum) - adjVal;
          newValSum = tempSum;
        }
      }
      modelVector[i*nModelColumns + j] += newValSum;
    }
    } // end omp parallel section
  }
  
  
  // 3. Optional generation of oversampled sub-image and convolution with oversampled PSF
  if (oversampledRegionsExist)
    for (n = 0; n < nOversampledRegions; n++)
      oversampledRegionsVect[n]->ComputeRegionAndDownsample(modelVector, functionObjects, nFunctions);
  
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
// unless you are aware that it will NOT return the full (expanded) model image.)
double * ModelObject::GetSingleFunctionImage( double params[], int functionIndex )
{
  double  x0, y0, x, y, newVal;
  int  offset = 0;
  int  iDataRow, iDataCol;
  long  i, j, z, zModel;
  vector<FunctionObject *> singleFuncObjVector;
  
  assert( (functionIndex >= 0) );
  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    fprintf(stderr, "** ModelObject::SingleFunctionImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g, ", parameterLabels[0].c_str(), params[0]);
    for (int n = 1; n < nParamsTot; n++)
      printf(", %s = %g", parameterLabels[n].c_str(), params[n]);
    printf("\n");
#endif
    fprintf(stderr, "Exiting ...\n\n");
    exit(-1);
  }

  // Separate out individual-component parameters and tell associated function objects 
  // to do setup. The first component's parameters start at params[0]; the second's 
  // at params[paramSizes[0]], the third at params[paramSizes[0] + paramSizes[1]], etc.
  for (int np = 0; np < nFunctions; np++) {
    if (fsetStartFlags[np] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
    }
    functionObjects[np]->Setup(params, offset, x0, y0);
    offset += paramSizes[np];
  }
  
  
  // If requested function is PointSource, re-assign the PsfInterpolator object
  // (see comments in CreateModelImage for why we need to do this)
  if (functionObjects[functionIndex]->IsPointSource())
    functionObjects[functionIndex]->AddPsfInterpolator(psfInterpolator);

  // 1. OK, populate modelVector with the model image -- standard pixel scaling
  // OpenMP Parallel section; see CreateModelImage() for general notes on this
  // Note that since we expect this code to be called only occasionally, we have
  // not converted it to the fast-for-small-images, single-loop version used in
  // CreateModelImages()
#pragma omp parallel private(i,j,x,y,newVal)
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
  
  // 2. Do PSF convolution, if requested and if this is *not* a PointSource function
  if ((doConvolution) && (! functionObjects[functionIndex]->IsPointSource()))
    psfConvolver->ConvolveImage(modelVector);

  // 3. Optional generation of oversampled sub-image and convolution with oversampled PSF
  if (oversampledRegionsExist)
    for (int n = 0; n < nOversampledRegions; n++) {
      // populate FunctionObject vector with just this function
      singleFuncObjVector.push_back(functionObjects[functionIndex]);
      oversampledRegionsVect[n]->ComputeRegionAndDownsample(modelVector, singleFuncObjVector, 1);
    }

  // Return model image (extract correct subimage if PSF convolution was done)
  if (doConvolution) {
    if (! outputModelVectorAllocated) {
      outputModelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
      if (outputModelVector == NULL) {
        fprintf(stderr, "*** ERROR: Unable to allocate memory for output model image!\n");
        fprintf(stderr, "    (Requested image size was %ld pixels)\n", nDataVals);
        return NULL;
      }
      outputModelVectorAllocated = true;
    }
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in output image
    for (z = 0; z < nDataVals; z++) {
      iDataRow = z / nDataColumns;
      iDataCol = z - (long)iDataRow * (long)nDataColumns;
      zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
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
  int  iDataRow, iDataCol;
  long  z, zModel;
  double  totalFlux, noise_squared;
	
  if (doConvolution) {
    for (z = 0; z < nDataVals; z++) {
      if ( (! maskExists) || (maskExists && (maskVector[z] > 0)) ) {
        // only update values that aren't masked out
        // (don't rely on previous weightVector[z] value, since sometimes model flux
        // might be zero for an unmasked pixel)
        iDataRow = z / nDataColumns;
        iDataCol = z - (long)iDataRow * (long)nDataColumns;
        zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
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
double ModelObject::ComputePoissonMLRDeviate( long i, long i_model )
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
  int  iDataRow, iDataCol;
  long  z, zModel, b, bModel;
  
#ifdef DEBUG
  printf("ComputeDeviates: Input parameters: ");
  for (int np = 0; np < nParamsTot; np++)
    printf("p[%d] = %g, ", np, params[np]);
  printf("\n");
#endif

  CreateModelImage(params);
  if (modelErrors)
    UpdateWeightVector();


  // In standard case, z = index into dataVector, weightVector, and yResults; it comes 
  // from linearly stepping through (0, ..., nDataVals).
  // In the bootstrap case, z = index into yResults and bootstrapIndices vector;
  // b = bootstrapIndices[z] = index into dataVector and weightVector
  
  // NOTE: in the doConvolution case, the algorithm is sufficiently complicated that
  // it makes sense to keep it in the current form, with PMLR or chi^2 calculation
  // for each pixel being done with if/else statements
  
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in data and weight images (excluding the outer borders of the model image,
    // which are only for ensuring proper PSF convolution)
    if (doBootstrap) {
      for (z = 0; z < nValidDataVals; z++) {
        b = bootstrapIndices[z];
        iDataRow = b / nDataColumns;
        iDataCol = b - (long)iDataRow * (long)nDataColumns;
        bModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        if (poissonMLR)
          yResults[z] = ComputePoissonMLRDeviate(b, bModel);
        else   // standard chi^2 term
          yResults[z] = weightVector[b] * (dataVector[b] - modelVector[bModel]);
      }
    }
    else {
      for (z = 0; z < nDataVals; z++) {
        iDataRow = z / nDataColumns;
        iDataCol = z - (long)iDataRow * (long)nDataColumns;
        zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        if (poissonMLR)
          yResults[z] = ComputePoissonMLRDeviate(z, zModel);
        else   // standard chi^2 term
          yResults[z] = weightVector[z] * (dataVector[z] - modelVector[zModel]);
      }
    }
  }   // end if convolution case
  
  // NOTE: in the non-convolution case, the algorithm is simple enough that we can
  // write it as separate loops for the PMLR vs chi^2 cases (which allows the
  // compiler to auto-vectorize the loops in the chi^2 case)
  else {
    // No convolution, so model image is same size & shape as data and weight images
    if (doBootstrap) {
      if (poissonMLR)
        for (z = 0; z < nValidDataVals; z++) {
          b = bootstrapIndices[z];
          yResults[z] = ComputePoissonMLRDeviate(b, b);
        } else {   // standard chi^2 term
        for (z = 0; z < nValidDataVals; z++) {
          b = bootstrapIndices[z];
          yResults[z] = weightVector[b] * (dataVector[b] - modelVector[b]);
        }
       }
    }
    else {
      if (poissonMLR)
        for (z = 0; z < nDataVals; z++)
          yResults[z] = ComputePoissonMLRDeviate(z, z);
      else   // standard chi^2 term
        for (z = 0; z < nDataVals; z++)
          yResults[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
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
      fprintf(stderr, "    (Requested image size was %ld pixels)\n", nDataVals);
      return -1;
    }
    weightVectorAllocated = true;
  }
  else {
    fprintf(stderr, "WARNING: ModelImage::UseModelErrors -- weight vector already allocated!\n");
  }

  for (long z = 0; z < nDataVals; z++) {
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
      fprintf(stderr, "    (Requested vector size was %ld pixels)\n", nDataVals);
      return -1;
    }
    extraCashTermsVectorAllocated = true;
  }
  else {
    fprintf(stderr, "WARNING: ModelImage::UseCashStatistic -- extra-terms vector already allocated!\n");
  }

  for (long z = 0; z < nDataVals; z++) {
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
  int  iDataRow, iDataCol;
  long  z, zModel, b, bModel;
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
        iDataCol = b - (long)iDataRow * (long)nDataColumns;
        bModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
        deviatesVector[z] = weightVector[b] * (dataVector[b] - modelVector[bModel]);
      }
    } else {
      for (z = 0; z < nDataVals; z++) {
        iDataRow = z / nDataColumns;
        iDataCol = z - (long)iDataRow * (long)nDataColumns;
        zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
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
      // Note: this loop is auto-vectorized when compiling with -O3 and -sse2 (g++-7)
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
/// Computes and returns the Cash statistic for the current set of model parameters
/// (computes model image and compares it with the data image)
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
  int  iDataRow, iDataCol;
  long  z, zModel, b, bModel;
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
        iDataCol = b - (long)iDataRow * (long)nDataColumns;
        bModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
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
        iDataCol = z - (long)iDataRow * (long)nDataColumns;
        zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
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
/// Prints the number of data values (pixels) in the data image
void ModelObject::PrintDescription( )
{
  // Don't test for verbose level, since we assume user only calls this method
  // if they *want* printed output
  printf("Model Object: %ld data values (pixels)\n", nDataVals);
}


/* ---------------- PUBLIC METHOD: GetFunctionNames ------------------- */
/// Adds the names of image functions in the model (calling each function's 
/// GetShortName method) to the input vector of strings
void ModelObject::GetFunctionNames( vector<string>& functionNames )
{
  functionNames.clear();
  for (FunctionObject *funcObj : functionObjects)
    functionNames.push_back(funcObj->GetShortName());
}


/* ---------------- PUBLIC METHOD: GetFunctionLabels ------------------ */
/// Adds the names of image functions in the model (calling each function's 
/// GetShortName method) to the input vector of strings
void ModelObject::GetFunctionLabels( vector<string>& functionLabels )
{
  functionLabels.clear();
  for (FunctionObject *funcObj : functionObjects)
    functionLabels.push_back(funcObj->GetLabel());
}


/* ---------------- PUBLIC METHOD: PrintModelParamsToStrings ---------- */
/// Prints description of model (e.g., best-fit result) as strings to the input
/// vector of string. 
/// Optionally, the lower and upper limits defined in parameterInfo are also printed, 
/// OR associated lower and upper error bounds in errs can be printed.
///
/// If errs != NULL, then +/- errors are printed as well (only if printLimits is false)
///
/// If prefix != NULL, then the specified character (e.g., '#') is prepended to
/// each output line.
///
/// If printLimits == true, then lower and upper parameter limits will be printed
/// for each parameter (or else "fixed" for fixed parameters)
int ModelObject::PrintModelParamsToStrings( vector<string> &stringVector, double params[], 
									double errs[], const char *prefix, bool printLimits )
{
  double  x0, y0, paramVal;
  int nParamsThisFunc, k;
  int  indexOffset = 0;
  string  funcName, funcLabel, paramName, newLine;
  bool  thisFuncHasUnits = false;
  string  unitsPrefix;
  vector<string>  paramUnits;
  vector<string>  extraParamLines;

  if ((printLimits) && (parameterInfoVect.size() == 0)) {
    fprintf(stderr, "** ERROR: ModelObject::PrintModelParamsToStrings -- printing of parameter limits\n");
    fprintf(stderr, "was requested, but parameterInfoVect is empty!\n");
    return -1;
  }

  for (int n = 0; n < nFunctions; n++) {
    if (fsetStartFlags[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      k = indexOffset;
      x0 = params[k] + imageOffset_X0;
      y0 = params[k + 1] + imageOffset_Y0;
      // blank line (or line with just prefix character) before start of function set
      stringVector.push_back(PrintToString("%s\n", prefix));
      if (printLimits) {
        if (parameterInfoVect[k].fixed == 1)
          newLine = PrintToString(XY_FORMAT_WITH_FIXED, prefix, "X0", x0);
        else
          newLine = PrintToString(XY_FORMAT_WITH_LIMITS, prefix, "X0", x0, 
        				parameterInfoVect[k].limits[0], parameterInfoVect[k].limits[1]);
        stringVector.push_back(newLine);
        if (parameterInfoVect[k + 1].fixed == 1)
          newLine = PrintToString(XY_FORMAT_WITH_FIXED, prefix, "Y0", y0);
        else
          newLine = PrintToString(XY_FORMAT_WITH_LIMITS, prefix, "Y0", y0, 
        				parameterInfoVect[k + 1].limits[0], parameterInfoVect[k + 1].limits[1]);
        stringVector.push_back(newLine);
      } else {
        if (errs != NULL) {
          stringVector.push_back(PrintToString(XY_FORMAT_WITH_ERRS, prefix, "X0", x0, errs[k]));
          stringVector.push_back(PrintToString(XY_FORMAT_WITH_ERRS, prefix, "Y0", y0, errs[k + 1]));
        } else {
          stringVector.push_back(PrintToString(XY_FORMAT, prefix, "X0", x0));
          stringVector.push_back(PrintToString(XY_FORMAT, prefix, "Y0", y0));
        }
      }
      indexOffset += 2;
    }
    
    // Now print the function and its parameters (and units, if they exist)
    funcName = functionObjects[n]->GetShortName();
    funcLabel = functionObjects[n]->GetLabel();
    if (! funcLabel.empty())
      funcLabel = PrintToString("   # LABEL %s", funcLabel.c_str());
    thisFuncHasUnits = false;
    if (functionObjects[n]->HasParameterUnits()) {
      paramUnits.clear();
      functionObjects[n]->GetParameterUnits(paramUnits);
      thisFuncHasUnits = true;
    }
    stringVector.push_back(PrintToString("%sFUNCTION %s%s\n", prefix, funcName.c_str(),
    						funcLabel.c_str()));
    						
    // print optional parameters, if they were set
    if (functionObjects[n]->ExtraParamsSet()) {
      functionObjects[n]->GetExtraParamsDescription(extraParamLines);
      for (int p = 0; p < extraParamLines.size(); p++)
        stringVector.push_back(PrintToString("%s%s\n", prefix, extraParamLines[p].c_str()));
    }
    
    // print the function parameter values
    nParamsThisFunc = paramSizes[n];
    for (int i = 0; i < nParamsThisFunc; i++) {
      paramName = GetParameterName(indexOffset + i);
      paramVal = params[indexOffset + i];
      // the following is placed in front of the unit string *if* we're not printing
      // the parameter errors (since that case we've already put a comment character
      // into the output line)
      unitsPrefix = "\t# ";
      if (printLimits)
        if (parameterInfoVect[indexOffset + i].fixed == 1)
          newLine = PrintToString(PARAM_FORMAT_WITH_FIXED, prefix, paramName.c_str(), 
        						paramVal);
        else
          newLine = PrintToString(PARAM_FORMAT_WITH_LIMITS, prefix, paramName.c_str(), 
        						paramVal, parameterInfoVect[indexOffset + i].limits[0], 
        						parameterInfoVect[indexOffset + i].limits[1]);
      else if (errs != NULL) {
        newLine = PrintToString(PARAM_FORMAT_WITH_ERRS, prefix, paramName.c_str(), 
        						paramVal, errs[indexOffset + i]);
        unitsPrefix = "";
      }
      else
        newLine = PrintToString(PARAM_FORMAT, prefix, paramName.c_str(), paramVal);
      if (thisFuncHasUnits) {
        string  unitsString = paramUnits[i];
        if (unitsString.size() > 0) {
          string  extraString = PrintToString(UNITS_FORMAT, unitsPrefix.c_str(), 
        										unitsString.c_str());
          newLine += extraString;
        }
      }
      stringVector.push_back(newLine + "\n");
    }
    indexOffset += paramSizes[n];
  }
  
  return 0;
}


/* ---------------- PUBLIC METHOD: PrintModelParamsHorizontalString --- */
/// Like PrintModelParamsToString, but prints parameter values all in one line to
/// a string (*without* parameter names or limits or errors), which is returned.
/// Meant to be used in printing results of bootstrap resampling (imfit) or MCMC
/// chains (imfit-mcmc)
string ModelObject::PrintModelParamsHorizontalString( const double params[], const string& separator )
{
  double  x0, y0, paramVal;
  int nParamsThisFunc, k;
  int  indexOffset = 0;
  string  outputString;

  for (int n = 0; n < nFunctions; n++) {
    if (fsetStartFlags[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      k = indexOffset;
      x0 = params[k] + imageOffset_X0;
      y0 = params[k + 1] + imageOffset_Y0;
      // Very first parameter on the line (X0 for first function set) should *not* 
      // be preceded by separator
      if (n == 0)
        outputString += PrintToString("%#.10g%s%#.10g", x0, separator.c_str(), y0);
      else
        outputString += PrintToString("%s%#.10g%s%#.10g", separator.c_str(), x0, separator.c_str(), y0);
      indexOffset += 2;
    }

    // Now print the function and its parameters
    nParamsThisFunc = paramSizes[n];
    for (int i = 0; i < nParamsThisFunc; i++) {
      paramVal = params[indexOffset + i];
      outputString += PrintToString("%s%#.10g", separator.c_str(), paramVal);
    }
    indexOffset += paramSizes[n];
  }

  return outputString;
}


/* ---------------- PUBLIC METHOD: GetImageOffsets --------------------- */
/// Given an input array (all zeros), this method returns the array with
/// locations corresponding to each function set's imageOffset_X0 and 
/// imageOffset_Y0 values properly set (if there are no image offsets, then
/// the corresponding values will remain zero).
/// Intended for use in bootstrap_errors.cpp.
void ModelObject::GetImageOffsets( double params[] )
{
  int  indexOffset = 0;
  int  k;
  
  // assume input array is all zeros to begin with; only
  // modify it if there really are image offsets
  if ((imageOffset_X0 != 0) || (imageOffset_Y0 != 0)) {
    for (int n = 0; n < nFunctions; n++) {
      if (fsetStartFlags[n] == true) {
        k = indexOffset;
        params[k] = imageOffset_X0;
        params[k + 1] = imageOffset_Y0;
        indexOffset += 2;
      }
      indexOffset += paramSizes[n];
    }
  }
}


/* ---------------- PUBLIC METHOD: GetParamHeader ---------------------- */
/// Prints all function and parameter names in order all on one line; e.g., for use as 
/// header in bootstrap-parameters output file.
string ModelObject::GetParamHeader( )
{
  int nParamsThisFunc, nSet;
  int  indexOffset = 0;
  string  paramName, headerLine, newString;

  headerLine = "# ";
  nSet = 0;
  for (int n = 0; n < nFunctions; n++) {
    if (fsetStartFlags[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      nSet += 1;
      newString = PrintToString("X0_%d\t\tY0_%d\t\t", nSet, nSet);
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
/// Returns the status from MakeBootstrapSample(), which will be -1 if memory
/// allocation for the bootstrap-indices vector failed.
int ModelObject::UseBootstrap( )
{
  int  status = 0;
  
  doBootstrap = true;
  // Note that this is slightly inefficient: we don't really *need* to generate
  // a bootstrap sample right now, since we will call MakeBootstrapSample directly
  // later on, every time we need a new sample. But calling this now *does* force
  // allocation of the bootstrapIndices array....
  status = MakeBootstrapSample();
  return status;
}


/* ---------------- PUBLIC METHOD: MakeBootstrapSample ----------------- */
/// Generate a new bootstrap resampling of the data (more precisely, this generate a
/// bootstrap resampling of the data *indices*)
/// Returns -1 if memory allocation for the bootstrap indices vector failed,
/// otherwise returns 0.
int ModelObject::MakeBootstrapSample( )
{
  long  n;
  bool  badIndex;
  
  if (! bootstrapIndicesAllocated) {
    bootstrapIndices = (long *) calloc((size_t)nValidDataVals, sizeof(long));
    if (bootstrapIndices == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for bootstrap-resampling pixel indices!\n");
      fprintf(stderr, "    (Requested vector size was %ld pixels)\n", nValidDataVals);
      return -1;
    }
    bootstrapIndicesAllocated = true;
  }
  for (long i = 0; i < nValidDataVals; i++) {
    // pick random data point between 0 and nDataVals - 1, inclusive;
    // reject masked pixels
    badIndex = true;
    do {
      n = (long)floor( genrand_real2()*nDataVals );
      if (weightVector[n] > 0.0)
        badIndex = false;
    } while (badIndex);
    bootstrapIndices[i] = n;
  }
  return 0;
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
      printf(" %f", pixelVector[(long)i * (long)nColumns + j]);
    printf("\n");
  }
  printf("\n");
}


/* ---------------- PUBLIC METHOD: PrintInputImage -------------------- */
/// Prints the input data image to stdout (for debugging purposes).
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
/// Prints the current computed model image to stdout (for debugging purposes).
void ModelObject::PrintModelImage( )
{

  if (! modelImageComputed) {
    fprintf(stderr, "* ModelObject::PrintModelImage -- Model image has not yet been computed!\n\n");
    return;
  }
  printf("The model image, row by row:\n");
  PrintImage(modelVector, nModelColumns, nModelRows);
}


/* ---------------- PUBLIC METHOD: PrintMask ------------------------- */
/// Prints the input mask image to stdout (for debugging purposes).
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
/// Prints the current weight image to stdout (for debugging purposes).
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
/// Note that this is usually called by AddFunctions() in add_functions.cpp.

void ModelObject::PopulateParameterNames( )
{
  int  n;

  for (n = 0; n < nFunctions; n++) {
    if (fsetStartFlags[n] == true) {
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
/// Prints the total number of image functions (instances of FunctionObject 
/// subclasses) making up the model.

int ModelObject::GetNFunctions( )
{
  return nFunctions;
}


/* ---------------- PUBLIC METHOD: GetNParams -------------------------- */
/// Prints the total number of parameters making up the model.
int ModelObject::GetNParams( )
{
  return nParamsTot;
}


/* ---------------- PUBLIC METHOD: GetNDataValues ---------------------- */
/// Prints the number of data values (pixels, masked or unmasked) in the data image.
long ModelObject::GetNDataValues( )
{
  return nDataVals;
}


/* ---------------- PUBLIC METHOD: GetNValidPixels --------------------- */
/// Prints the number of *valid* (i.e., unmasked) data values (pixels) in the data image.
long ModelObject::GetNValidPixels( )
{
  return nValidDataVals;
}


/* ---------------- PUBLIC METHOD: HasPSF ------------------------------ */
/// Returns true if the model has a PSF image
bool ModelObject::HasPSF( )
{
  return doConvolution;
}

/* ---------------- PUBLIC METHOD: HasOversampledPSF ------------------- */
/// Returns true if the model has one or more oversampled regions (with oversampled PSFs).
bool ModelObject::HasOversampledPSF( )
{
  return oversampledRegionsExist;
}


/* ---------------- PUBLIC METHOD: HasMask ----------------------------- */
/// Returns true if a mask image exists.
bool ModelObject::HasMask( )
{
  return maskExists;
}


/* ---------------- PUBLIC METHOD: GetModelImageVector ----------------- */
/// Returns a pointer to the model image (matching the data image in size if
/// convolution is being done).
/// If the model image has not yet been computed, returns NULL.
double * ModelObject::GetModelImageVector( )
{
  int  iDataRow, iDataCol;
  long  z, zModel;

  if (! modelImageComputed) {
    fprintf(stderr, "* ModelObject::GetModelImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  if (doConvolution) {
    if (! outputModelVectorAllocated) {
      outputModelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
      if (outputModelVector == NULL) {
        fprintf(stderr, "*** ERROR: Unable to allocate memory for output model image!\n");
        fprintf(stderr, "    (Requested image size was %ld pixels)\n", nDataVals);
        return NULL;
      }
      outputModelVectorAllocated = true;
    }
    // Step through (previously computed) model image so that we correctly match 
    // its pixels with corresponding pixels output image
    for (z = 0; z < nDataVals; z++) {
      iDataRow = z / nDataColumns;
      iDataCol = z - (long)iDataRow * (long)nDataColumns;
      zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
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
/// Computes residual image (data - model) using current model image (usually
/// best-fit image), and returns pointer to residual-image vector.
/// Returns NULL if memory allocation failed.
double * ModelObject::GetResidualImageVector( )
{
  int  iDataRow, iDataCol;
  long  z, zModel;

  if (! modelImageComputed) {
    fprintf(stderr, "* ModelObject::GetResidualImageVector -- Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  // WARNING: If we are calling this function for a second or subsequent time,
  // nDataVals *might* have changed; we are currently assuming it hasn't!
  if (! residualVectorAllocated) {
    residualVector = (double *) calloc((size_t)nDataVals, sizeof(double));
    if (residualVector == NULL) {
      fprintf(stderr, "*** ERROR: Unable to allocate memory for output residual image!\n");
      fprintf(stderr, "    (Requested image size was %ld pixels)\n", nDataVals);
      return NULL;
    }
    residualVectorAllocated = true;
  }
    
  if (doConvolution) {
    // Step through model image so that we correctly match its pixels with corresponding
    // pixels in data and weight images
    for (z = 0; z < nDataVals; z++) {
      iDataRow = z / nDataColumns;
      iDataCol = z - (long)iDataRow * (long)nDataColumns;
      zModel = (long)nModelColumns * (long)(nPSFRows + iDataRow) + nPSFColumns + iDataCol;
      residualVector[z] = (dataVector[z] - modelVector[zModel]);
    }
  }
  else {
    // Model image is same size & shape as data and weight images
    // Note: this loop is auto-vectorized when compiling with -O3 and -sse2 (g++-7)
    for (z = 0; z < nDataVals; z++) {
      residualVector[z] = (dataVector[z] - modelVector[z]);
    }
  }

  return residualVector;
}


/* ---------------- PUBLIC METHOD: GetWeightImageVector ---------------- */
/// Returns the weightVector converted to 1/sigma^2 (i.e., 1/variance) form.
/// Returns NULL if memory allocation failed.
double * ModelObject::GetWeightImageVector( )
{
  if (! weightValsSet) {
    fprintf(stderr, "* ModelObject::GetWeightImageVector -- Weight image has not yet been computed!\n\n");
    return NULL;
  }
  
  standardWeightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  if (standardWeightVector == NULL) {
    fprintf(stderr, "*** ERROR: Unable to allocate memory for output weight image!\n");
    fprintf(stderr, "    (Requested image size was %ld pixels)\n", nDataVals);
    return NULL;
  }
  for (long z = 0; z < nDataVals; z++) {
    // Note: this loop is auto-vectorized when compiling with -O3 and -sse2 (g++-7)
    double  w_sqrt = weightVector[z];   // internal weight value (sqrt of formal weight)
  	standardWeightVector[z] = w_sqrt*w_sqrt;
  }
  standardWeightVectorAllocated = true;
  return standardWeightVector;
}


/* ---------------- PUBLIC METHOD: GetDataVector ----------------------- */
/// Returns a pointer to the data image.
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
    if (fsetStartFlags[n] == true) {
      // start of new function set: skip over existing x0,y0 values
      offset += 2;
    }
    functionObjects[n]->Setup(params, offset, x0_all, y0_all);
    offset += paramSizes[n];
  }

  // NOTE: I tried implementing the Kahan summation algorithm for this
  // calculation; it increased the time for a single-component Gaussian
  // model by a factor of 2--2.5, and did *not* produce any noticeable
  // difference in output flux (even for a Gaussian with sigma = 500 pix).
  // So: probably not useful to use Kahan summation here (better accuracy
  // can be achieved by using native TotalFlux() methods of function
  // objects).
  
  totalModelFlux = 0.0;
  // Integrate over the image, once per function
  for (n = 0; n < nFunctions; n++) {
    if (functionObjects[n]->IsBackground())
      totalComponentFlux = 0.0;
    else {
      if (functionObjects[n]->CanCalculateTotalFlux()) {
        totalComponentFlux = functionObjects[n]->TotalFlux();
        if (verboseLevel > 0)
          printf("\tUsing %s.TotalFlux() method...\n", functionObjects[n]->GetShortName().c_str());
      } else {
        totalComponentFlux = 0.0;
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
    } // end else [calculate total flux for non-background component]
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
  
  for (long z = 0; z < nDataVals; z++) {
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
  long  z;
  
  // check individual pixels in weightVector, but only if they aren't masked by maskVector
  if (maskExists) {
    for (z = 0; z < nDataVals; z++) {
      if (maskVector[z] > 0.0) {
        if (! isfinite(weightVector[z]))
          nonFinitePixels = true;
        else if (weightVector[z] < 0.0)
          negativePixels = true;
      }
    }  
  }
  else {
    for (z = 0; z < nDataVals; z++) {
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




// Extra stuff

/// Given an input PSF-image vector, this function normalizes it in place.
void NormalizePSF( double *psfPixels, long nPixels_psf )
{
  // Use Kahan summation to avoid underflow
  long  k;
  double  psfSum = 0.0, storedError = 0.0, adjustedVal = 0.0, tempSum = 0.0;
  for (k = 0; k < nPixels_psf; k++) {
    adjustedVal = psfPixels[k] - storedError;
    tempSum = psfSum + adjustedVal;
    storedError = (tempSum - psfSum) - adjustedVal;
    psfSum = tempSum;
  }
  for (k = 0; k < nPixels_psf; k++)
    psfPixels[k] = psfPixels[k] / psfSum;
}




/* END OF FILE: model_object.cpp --------------------------------------- */
