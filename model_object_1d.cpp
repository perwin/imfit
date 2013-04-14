/* FILE: model_object1d.cpp -------------------------------------------- */
/* VERSION 0.1
 *
 * This is intended to be an abstract base class for the various
 * "model" objects (e.g., image data + fitting functions).
 * 
 *
 * length of: dataVector, weightVector = nDataVals
 * length of: modelVector = nModelVals = nDataVals IF NO PSF
 * length of: modelVector = nModelVals = nDataVals + 2*nPSFVals IF PSF USED
 *
 *     [v0.1]: 27 Nov 2009: Created; initial development.
 *
 */


/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "definitions.h"
#include "model_object_1d.h"
#include "mp_enorm.h"
#include "convolver1d.h"
#include "mersenne_twister.h"


/* ---------------- Definitions ---------------------------------------- */
static string  UNDEFINED = "<undefined>";

#define OPENMP_CHUNK_SIZE  10


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModelObject1d::ModelObject1d( )
{
  dataValsSet = weightValsSet = false;
  parameterBoundsSet = false;
  parameterBounds = NULL;
  modelVector = NULL;
  modelVectorAllocated = false;
  modelImageComputed = false;
  maskExists = false;
  dataAreMagnitudes = true;
  doBootstrap = false;
  bootstrapIndicesAllocated = false;
  zeroPointSet = false;
  nFunctions = 0;
  nFunctionSets = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
  dataStartOffset = 0;
  debugLevel = 0;
  nCombined = 1;
  zeroPoint = 0.0;
}


/* ---------------- PUBLIC METHOD: DefineFunctionSets ------------------ */
// We have to redefine this function from the ModelObject base function because
// nParamsTot is calculated differently
void ModelObject1d::DefineFunctionSets( vector<int>& functionStartIndices )
{
  int  nn, i;
  
  nFunctionSets = functionStartIndices.size();
    // define array of [false, false, false, ...]
  setStartFlag = (bool *)calloc(nFunctions, sizeof(bool));
  for (i = 0; i < nFunctionSets; i++) {
    nn = functionStartIndices[i];
    // function number n is start of new function set; 
    // change setStartFlag[n] to true
    setStartFlag[nn] = true;
  }
  
  // total number of parameters = number of parameters for individual functions
  // plus x0 for each function set
  nParamsTot = nFunctionParams + nFunctionSets;
}



/* ---------------- PUBLIC METHOD: SetZeroPoint ----------------------- */

void ModelObject1d::SetZeroPoint( double zeroPointValue )
{
  zeroPoint = zeroPointValue;
  zeroPointSet = true;
  if (nFunctions < 1) {
    fprintf(stderr, "ModelObject1d: WARNING: zero point added to model object");
    fprintf(stderr, " before any functions were added!\n");
    return;
  }
  else {
    for (int n = 0; n < nFunctions; n++)
      functionObjects[n]->SetZeroPoint(zeroPoint);
  }
}


/* ---------------- PUBLIC METHOD: AddDataVectors --------------------- */

void ModelObject1d::AddDataVectors( int nDataValues, double *xValVector, 
										double *yValVector, bool magnitudeData )
{
  nModelVals = nDataVals = nValidDataVals = nDataValues;
  modelXValues = dataXValues = xValVector;
  dataVector = yValVector;
  dataValsSet = true;
  dataAreMagnitudes = magnitudeData;  // are yValVector data magnitudes?

  modelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  modelVectorAllocated = true;
}


/* ---------------- PUBLIC METHOD: AddErrorVector1D -------------------- */

void ModelObject1d::AddErrorVector1D( int nDataValues, double *inputVector,
                                      int inputType )
{
  assert (nDataValues == nDataVals);
  weightVector = inputVector;
  
  // Convert noise values into weights, if needed
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
      // do nothing, since input values *are* weights
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


/* ---------------- PUBLIC METHOD: AddMaskVector1D --------------------- */
// Code for adding and processing a vector containing the 1-D mask.
// Note that although our default *input* format is "0 = good pixel, > 0 =
// bad pixel", internally we convert all bad pixels to 0 and all good pixels
// to 1, so that we can multiply the weight vector by the (internal) mask values.
void ModelObject1d::AddMaskVector1D( int nDataValues, double *inputVector,
                                      int inputType )
{
  assert (nDataValues == nDataVals);

  maskVector = inputVector;
  nValidDataVals = 0;   // Since there's a mask, not all pixels from the original
                        // profile will be valid
    
  // We need to convert the mask values so that good pixels = 1 and bad
  // pixels = 0.
  switch (inputType) {
    case MASK_ZERO_IS_GOOD:
      // This is our "standard" input mask: good pixels are zero, bad pixels
      // are positive integers
      printf("ModelObject1D::AddMaskVector -- treating zero-valued pixels as good ...\n");
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
      printf("ModelObject1D::AddMaskVector -- WARNING: unknown inputType detected!\n\n");
      exit(-1);
  }
      
  maskExists = true;
}


/* ---------------- PUBLIC METHOD: AddPSFVector1D ---------------------- */
// This function needs to be redefined because the base function in ModelObject
// assumes a 2-D PSF.
// NOTE: PSF vector y-values are assumed to be intensities, *not* magnitudes!
void ModelObject1d::AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector )
{
  nPSFVals = nPixels_psf;
  
  // Do full setup for convolution
  // 1. Figure out extra size for model profile (PSF dimensions added to each end)
  nModelVals = nDataVals + 2*nPSFVals;
  dataStartOffset = nPSFVals;
  // 2. Create new model vector and set dataStartOffset to nPSFVals
  if (modelVectorAllocated)
    free(modelVector);
  modelVector = (double *) calloc((size_t)nModelVals, sizeof(double));
  modelVectorAllocated = true;
  // 3. Create new xVals vector for model
  modelXValues = (double *) calloc((size_t)nModelVals, sizeof(double));
  double  deltaX = dataXValues[1] - dataXValues[0];
  double  newXStart = dataXValues[0] - nPSFVals*deltaX;
  for (int i = 0; i < nPSFVals; i++)
    modelXValues[i] = newXStart + deltaX*i;
  for (int i = 0; i < nDataVals; i++)
    modelXValues[i + nPSFVals] = dataXValues[i];
  for (int i = 0; i < nPSFVals; i ++)
    modelXValues[i + nPSFVals + nDataVals] = dataXValues[nDataVals - 1] + deltaX*(i + 1);
  // 4. Create and setup Convolver1D object
  psfConvolver = new Convolver1D();
  psfConvolver->SetupPSF(yValVector, nPSFVals);
  psfConvolver->SetupProfile(nModelVals);
  psfConvolver->DoFullSetup(debugLevel);
  doConvolution = true;
}


/* ---------------- PUBLIC METHOD: CreateModelImage -------------------- */

void ModelObject1d::CreateModelImage( double params[] )
{
  double  x0, x, newVal;
  int  i, n, z;
  int  offset = 0;

  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    printf("** ModelObject1d::CreateModelImage -- non-finite values detected in parameter vector!\n");
#ifdef DEBUG
    printf("   Parameter values: %s = %g", parameterLabels[0].c_str(), params[0]);
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
      // start of new function set: extract x0 and then skip over them
      x0 = params[offset];
      offset += 1;
    }
    functionObjects[n]->Setup(params, offset, x0);
    offset += paramSizes[n];
  }
  
  // populate modelVector with the model
  // OpenMP Parallel Section
  int  chunk = OPENMP_CHUNK_SIZE;
// Note that we cannot specify modelVector as shared [or private] bcs it is part
// of a class (not an independent variable); happily, by default all references in
// an omp-parallel section are shared unless specified otherwise
#pragma omp parallel private(i,n,x,newVal)
  {
  #pragma omp for schedule (static, chunk)
  for (i = 0; i < nModelVals; i++) {
    x = modelXValues[i];
    newVal = 0.0;
    for (n = 0; n < nFunctions; n++)
      newVal += functionObjects[n]->GetValue(x);
    modelVector[i] = newVal;
#ifdef DEBUG
    printf("x = %g, newVal = %g,  ", x, newVal);
#endif
  }
  
  } // end omp parallel section

#ifdef DEBUG
    printf("   Parameter values: %s = %g", parameterLabels[0].c_str(), params[0]);
    for (z = 1; z < nParamsTot; z++)
      printf(", %s = %g", parameterLabels[z].c_str(), params[z]);
    printf("\n");
#endif

  // Do PSF convolution, if requested
  if (doConvolution) {
    psfConvolver->ConvolveProfile(modelVector);
  }
  
  // Convert to magnitudes, if required
  if (dataAreMagnitudes) {
    for (i = 0; i < nModelVals; i++) {
      modelVector[i] = zeroPoint - 2.5 * log10(modelVector[i]);
    }
  }
  
  modelImageComputed = true;
}


/* ---------------- PUBLIC METHOD: ComputeDeviates --------------------- */

void ModelObject1d::ComputeDeviates( double yResults[], double params[] )
{

#ifdef DEBUG
  printf("ComputeDeviates: Input parameters: ");
  for (int z = 0; z < nParamsTot; z++)
    printf("p[%d] = %g, ", z, params[z]);
  printf("\n");
#endif

  CreateModelImage(params);
  
  if (doBootstrap) {
    for (int z = 0; z < nDataVals; z++) {
      int i = bootstrapIndices[z];
      yResults[z] = weightVector[i] * (dataVector[i] - modelVector[dataStartOffset + i]);
    }
  } else {
    for (int z = 0; z < nDataVals; z++) {
      yResults[z] = weightVector[z] * (dataVector[z] - modelVector[dataStartOffset + z]);
#ifdef DEBUG
      printf("weight = %g, data = %g, model = %g ==> yResults = %g\n", weightVector[z], dataVector[z], modelVector[dataStartOffset + z], yResults[z]);
#endif
    }
  }

//   for (int z = 0; z < nDataVals; z++) {
//     yResults[z] = weightVector[z] * (dataVector[z] - modelVector[dataStartOffset + z]);
// #ifdef DEBUG
//     printf("weight = %g, data = %g, model = %g ==> yResults = %g\n", weightVector[z], dataVector[z], modelVector[dataStartOffset + z], yResults[z]);
// #endif
//   }
}


/* ---------------- PUBLIC METHOD: PrintDescription -------------------- */

void ModelObject1d::PrintDescription( )
{
  printf("ModelObject(1d): %d data values\n", nDataVals);
}


/* ---------------- PUBLIC METHOD: PrintModelParams --------=---------- */
// Basic function which prints to a file a summary of the best-fitting model,
// in form suitable for future use as an input config file.

void ModelObject1d::PrintModelParams( FILE *output_ptr, double params[], mp_par *parameterInfo,
																		double errs[] )
{
  double  x0, paramVal;
  int nParamsThisFunc, k;
  int  indexOffset = 0;
  string  funcName, paramName;

  for (int n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      k = indexOffset;
      x0 = params[k] + parameterInfo[k].offset;
      if (errs != NULL) {
        fprintf(output_ptr, "\nX0\t\t%f # +/- %f\n", x0, errs[k]);
      } else {
        fprintf(output_ptr, "\nX0\t\t%f\n", x0);
      }
      indexOffset += 1;
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


/* ---------------- PUBLIC METHOD: PopulateParameterNames -------------- */
// This function is redefined because the base function in ModelObject assumes
// two positional paramters ("X0").

void ModelObject1d::PopulateParameterNames( )
{
  int  n;

  for (n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0
      parameterLabels.push_back("X0");
    }
    functionObjects[n]->GetParameterNames(parameterLabels);
  }
}


/* ---------------- PUBLIC METHOD: UseBootstrap ------------------------ */
// Tells ModelObject1d object that from now on we'll operate in bootstrap
// resampling mode, so that bootstrapIndices vector is used to access the
// data and model values (and weight values, if any).

void ModelObject1d::UseBootstrap( )
{
  doBootstrap = true;
  MakeBootstrapSample();
}


/* ---------------- PUBLIC METHOD: MakeBootstrapSample ----------------- */
// Generate a new bootstrap resampling of the data (actually, we generate a
// bootstrap resampling of the data *indices*)

void ModelObject1d::MakeBootstrapSample( )
{
  
  if (! bootstrapIndicesAllocated) {
    bootstrapIndices = (int *) calloc((size_t)nDataVals, sizeof(int));
    bootstrapIndicesAllocated = true;
  }
  for (int i = 0; i < nDataVals; i++) {
    /* pick random data point between 0 and nDataVals - 1, inclusive */
    //n = round( (random()/MAX_RANDF)*(nDataVals - 1) );
    int n = (int)floor( genrand_real2()*nDataVals );
    bootstrapIndices[i] = n;
  }
}


/* ---------------- PUBLIC METHOD: GetModelProfile --------------------- */

int ModelObject1d::GetModelVector( double *profileVector )
{
  if (! modelImageComputed) {
    printf("* ModelObject: Model profile has not yet been computed!\n\n");
    return -1;
  }
  
  for (int z = 0; z < nDataVals; z++)
    profileVector[z] = modelVector[dataStartOffset + z];
  return nDataVals;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */
// Note that we have to turn various bool variables off and set nFunctions = 0,
// else we have problems when the *base* class (ModelObject) destructor is called!
ModelObject1d::~ModelObject1d()
{
  if (modelVectorAllocated) {
    free(modelVector);
    modelVectorAllocated = false;
  }
  if (doConvolution) {
    free(psfConvolver);
    free(modelXValues);
    doConvolution = false;
  }
  
  if (nFunctions > 0) {
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
    nFunctions = 0;
  }
  if (setStartFlag_allocated) {
    free(setStartFlag);
    setStartFlag_allocated = false;
  }
  
  if (bootstrapIndicesAllocated) {
    free(bootstrapIndices);
    bootstrapIndicesAllocated = false;
  }
}



/* END OF FILE: model_object1d.cpp ------------------------------------- */
