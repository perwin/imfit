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


/* ------------------------ Include Files (Header Files )--------------- */

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>

#include "definitions.h"
#include "model_object.h"
#include "mp_enorm.h"
#include "param_struct.h"


/* ---------------- Definitions ---------------------------------------- */
static string  UNDEFINED = "<undefined>";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModelObject::ModelObject( )
{
  dataValsSet = weightValsSet = false;
  parameterBoundsSet = false;
  parameterBounds = NULL;
  modelVector = NULL;
  residualVector = NULL;
  modelVectorAllocated = false;
  weightVectorAllocated = false;
  residualVectorAllocated = false;
  setStartFlag_allocated = false;
  modelImageComputed = false;
  maskExists = false;
  doConvolution = false;
  doChisquared = false;
  nFunctions = 0;
  nFunctionSets = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
  debugLevel = 0;
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


/* ---------------- PUBLIC METHOD: AddFunction ------------------------- */

void ModelObject::AddFunction( FunctionObject *newFunctionObj_ptr )
{
  int  nNewParams;
  
  functionObjects.push_back(newFunctionObj_ptr);
  nFunctions += 1;
  nNewParams = newFunctionObj_ptr->GetNParams();
  paramSizes.push_back(nNewParams);
  nFunctionParams += nNewParams;
//  newFunctionObj_ptr->GetParameterNames(parameterLabels);
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



/* ---------------- PUBLIC METHOD: AddDataVectors --------------------- */

void ModelObject::AddDataVectors( int nDataValues, double *xValVector, double *yValVector,
																	bool magnitudeData )
{
  // Just a placeholder for now (needs to be modified & overridden in derived
  // class ModelObject1D
  nDataVals = nDataValues;
}


/* ---------------- PUBLIC METHOD: SetZeroPoint ----------------------- */

void ModelObject::SetZeroPoint( double zeroPointValue )
{
  // Just a placeholder for now (needs to be modified & overridden in derived
  // class ModelObject1D
  ;
}


/* ---------------- PUBLIC METHOD: AddImageDataVector ------------------ */

void ModelObject::AddImageDataVector( double *pixelVector, int nImageColumns,
																			int nImageRows, int nCombinedImages )
{
  nColumns = nImageColumns;
  nRows = nImageRows;
  nDataVals = nValidDataVals = nColumns * nRows;
  dataVector = pixelVector;
  nCombined = nCombinedImages;
  nCombined_sqrt = sqrt(nCombined);
  dataValsSet = true;
  
  SetupModelImage(nDataVals, nColumns, nRows);
}


/* ---------------- PUBLIC METHOD: SetupModelImage -------------------- */
// Called by AddImageDataVector(); can also be used by itself in make-image
// mode. Tells ModelObject to allocate space for the model image.
void ModelObject::SetupModelImage( int nDataValues, int nImageColumns,
                                      int nImageRows )
{
  nDataVals = nDataValues;
  nColumns = nImageColumns;
  nRows = nImageRows;
  
  // Allocate modelimage vector
  modelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  modelVectorAllocated = true;
}


/* ---------------- PUBLIC METHOD: AddErrorVector ---------------------- */

void ModelObject::AddErrorVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  assert( (nDataValues == nDataVals) && (nImageColumns == nColumns) && 
          (nImageRows == nRows) );
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


/* ---------------- PUBLIC METHOD: AddErrorVector1D -------------------- */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
void ModelObject::AddErrorVector1D( int nDataValues, double *inputVector,
                                      int inputType )
{
  ;
}



/* ---------------- PUBLIC METHOD: AddMaskVector1D --------------------- */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
void ModelObject::AddMaskVector1D( int nDataValues, double *inputVector,
                                      int inputType )
{
  ;
}



/* ---------------- PUBLIC METHOD: GenerateErrorVector ----------------- */
// Generate an error vector based on Poisson statistics.
//    noise^2 = object_flux + sky + rdnoise^2
//
// Since sigma_adu = sigma_e/gain, we can go from
//    noise(e-)^2 = object_flux(e-) + sky(e-) + rdnoise^2
// to
//    noise(adu)^2 = object_flux(adu)/gain + sky(adu)/gain + rdnoise^2/gain^2
// (assuming that read noise is in units of e-, as is usual)
void ModelObject::GenerateErrorVector( double gain, double readNoise, double skyValue )
{
  double  sky_plus_readNoise, noise_squared, objectFlux;
  
  assert( (gain > 0.0) && (readNoise >= 0.0) );
  if ( skyValue <= 0.0 ) {
    printf("ModelObject::GenerateErrorVector -- original sky value (%g) is zero or negative.\n",
           skyValue);
  }

  // Allocate storage for weight image:
  weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  weightVectorAllocated = true;
  
  // Compute noise estimate for each pixel (see above for derivation)
  // Note that we assume a constant sky background (presumably already
  // subtracted)
  sky_plus_readNoise = skyValue/gain + readNoise*readNoise/(gain*gain);
  for (int z = 0; z < nDataVals; z++) {
    objectFlux = dataVector[z];
    if (objectFlux < 0.0)
      objectFlux = 0.0;
    noise_squared = objectFlux/gain + sky_plus_readNoise;
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
  assert( (nDataValues == nDataVals) && (nImageColumns == nColumns) && 
          (nImageRows == nRows) );

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
// This function must be called *after* SetupModelImage() is called (to ensure
// that we can tell the Convolver object the model image dimensions).
// [Note that SetupModelImage() is called by AddImageDataVector().]
void ModelObject::AddPSFVector(int nPixels_psf, int nColumns_psf, int nRows_psf,
                         double *psfPixels)
{

  if (modelVectorAllocated) {
    psfConvolver = new Convolver();
    psfConvolver->SetupPSF(psfPixels, nColumns_psf, nRows_psf);
    psfConvolver->SetupImage(nColumns, nRows);
    psfConvolver->DoFullSetup(debugLevel);
    doConvolution = true;
  }
  else {
    printf("WARNING: PSF image was supplied before main image dimensions!\n");
    printf("\tUnable to do convolution ...\n");
    doConvolution = false;
  }
}


/* ---------------- PUBLIC METHOD: AddPSFVector1D ---------------------- */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
void ModelObject::AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector )
{
  ;
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
  double  x0, y0, x, y, newVal;
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
  for (int n = 0; n < nFunctions; n++) {
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
  int  chunk = 100;
// Note that we cannot specify modelVector as shared [or private] bcs it is part
// of a class (not an independent variable); happily, by default all references in
// an omp-parallel section are shared unless specificied otherwise
#pragma omp parallel private(i,j,n,x,y,newVal)
  {
  #pragma omp for schedule (static, chunk)
  for (i = 0; i < nRows; i++) {   // step by row number = y
    y = (double)(i + 1);              // Iraf counting: first row = 1
    for (j = 0; j < nColumns; j++) {   // step by column number = x
      x = (double)(j + 1);                 // Iraf counting: first column = 1
      newVal = 0.0;
      for (n = 0; n < nFunctions; n++)
        newVal += functionObjects[n]->GetValue(x, y);
      modelVector[i*nColumns + j] = newVal;
    }
  }
  
  } // end omp parallel section
  
  
  // Do PSF convolution, if requested
  if (doConvolution)
    psfConvolver->ConvolveImage(modelVector);
  
  modelImageComputed = true;
}


/* ---------------- PUBLIC METHOD: ComputeDeviates --------------------- */
/* This function computes the vector of weighted deviates (differences between
 * model and data pixel values).  Note that a proper chi^2 sum requires *squaring*
 * each deviate value before summing them; we assume this is done elsewhere, by 
 * whatever function calls ComputeDeviates().
 */
void ModelObject::ComputeDeviates( double yResults[], double params[] )
{

#ifdef DEBUG
  printf("ComputeDeviates: Input parameters: ");
  for (int z = 0; z < nParamsTot; z++)
    printf("p[%d] = %g, ", z, params[z]);
  printf("\n");
#endif

  CreateModelImage(params);
  
  for (int z = 0; z < nDataVals; z++) {
    yResults[z] = nCombined_sqrt * weightVector[z] * (dataVector[z] - modelVector[z]);
// #ifdef DEBUG
//     printf("yResults[%d] = %g = %g * (%g - %g)  ", z, yResults[z], weightVector[z],
//     				dataVector[z], modelVector[z]);
// #endif
  }
}


/* ---------------- PUBLIC METHOD: SetupChisquaredCalcs ---------------- */
/* Function which tells object to prepare for making chi-square calculations
 * (i.e., allocate memory for self-stored deviates vector).
 */
void ModelObject::SetupChisquaredCalcs( )
{
  deviatesVector = (double *) malloc(nDataVals * sizeof(double));
  doChisquared = true;
}



/* ---------------- PUBLIC METHOD: ChiSquared -------------------------- */
/* Function for calculating chi^2 value for a model.
 *
 * IMPORTANT: SetupChisquaredCalcs() should be called (once) prior to any calls
 * to this function!
 */
double ModelObject::ChiSquared( double params[] )
{
  double  chi;
  
  CreateModelImage(params);
  
  for (int z = 0; z < nDataVals; z++)
    deviatesVector[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
  chi = mp_enorm(nDataVals, deviatesVector);
  
  return (nCombined*chi*chi);
}


/* ---------------- PUBLIC METHOD: PrintDescription ------------------- */

void ModelObject::PrintDescription( )
{
  printf("Model Object: %d data values (pixels)\n", nDataVals);
}


/* ---------------- PUBLIC METHOD: PrintModelParams --------=---------- */
// Basic function which prints to a file (or, e.g. stdout) a summary of the
// best-fitting model, in form suitable for future use as an input config file.
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
      x0 = params[k] + parameterInfo[k].offset;
      y0 = params[k + 1] + parameterInfo[k + 1].offset;
      if (errs != NULL) {
        fprintf(output_ptr, "\nX0\t\t%f # +/- %f\n", x0, errs[k]);
        fprintf(output_ptr, "Y0\t\t%f # +/- %f\n", y0, errs[k + 1]);
      } else {
        fprintf(output_ptr, "\nX0\t\t%f\n", x0);
        fprintf(output_ptr, "Y0\t\t%f\n", y0);
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
        fprintf(output_ptr, "%s\t\t%f # +/- %f\n", paramName.c_str(), paramVal, errs[indexOffset + i]);
      else
        fprintf(output_ptr, "%s\t\t%f\n", paramName.c_str(), paramVal);
    }
    indexOffset += paramSizes[n];
  }
}


/* ---------------- PUBLIC METHOD: PrintImage ------------------------- */
// Basic function which prints an image to stdout.  Mainly meant to be
// called by PrintInputImage, PrintModelImage, and PrintWeights.

void ModelObject::PrintImage( double *pixelVector )
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
    printf("* ModelObject: No image data supplied!\n\n");
    return;
  }
  printf("The whole input image, row by row:\n");
  PrintImage(dataVector);
}



/* ---------------- PUBLIC METHOD: PrintModelImage -------------------- */

void ModelObject::PrintModelImage( )
{

  if (! modelImageComputed) {
    printf("* ModelObject: Model image has not yet been computed!\n\n");
    return;
  }
  printf("The model image, row by row:\n");
  PrintImage(modelVector);
}


/* ---------------- PUBLIC METHOD: PrintWeights ----------------------- */

void ModelObject::PrintWeights( )
{

  if (! weightValsSet) {
    printf("* ModelObject: Weight vector has not yet been computed!\n\n");
    return;
  }
  printf("The weight image, row by row:\n");
  PrintImage(weightVector);
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
  if (! modelImageComputed) {
    printf("* ModelObject: Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  return modelVector;
}


/* ---------------- PUBLIC METHOD: GetResidualImageVector -------------- */

double * ModelObject::GetResidualImageVector( )
{
  if (! modelImageComputed) {
    printf("* ModelObject: Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  residualVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  residualVectorAllocated = true;
  for (int z = 0; z < nDataVals; z++)
    residualVector[z] = (dataVector[z] - modelVector[z]);
  return residualVector;
}


/* ---------------- PUBLIC METHOD: GetWeightImageVector ---------------- */

double * ModelObject::GetWeightImageVector( )
{
  if (! weightValsSet) {
    printf("* ModelObject: Weight image has not yet been computed!\n\n");
    return NULL;
  }
  
  return weightVector;
}


/* ---------------- PUBLIC METHOD: GetModelVector ---------------------- */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
int ModelObject::GetModelVector( double *profileVector )
{
  ;
}



/* ---------------- PUBLIC METHOD: UseBootstrap ------------------------ */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
void ModelObject::UseBootstrap( )
{
  ;
}



/* ---------------- PUBLIC METHOD: MakeBootstrapSample ----------------- */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
void ModelObject::MakeBootstrapSample( )
{
  ;
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
  if (doChisquared)
    free(deviatesVector);
  if (residualVectorAllocated)
    free(residualVector);
  
  if (nFunctions > 0)
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
  
  if (setStartFlag_allocated)
    free(setStartFlag);
  
  if (doConvolution)
    delete psfConvolver;
}



/* END OF FILE: model_object.cpp --------------------------------------- */
