/* FILE: model_object1d.cpp -------------------------------------------- */
/* VERSION 0.1
 *
 * This is intended to be an abstract base class for the various
 * "model" objects (e.g., image data + fitting functions).
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


/* ---------------- Definitions ---------------------------------------- */
static string  UNDEFINED = "<undefined>";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModelObject1d::ModelObject1d( )
{
  dataValsSet = weightValsSet = false;
  parameterBoundsSet = false;
  parameterBounds = NULL;
  modelVector = NULL;
  modelVectorAllocated = false;
  modelImageComputed = false;
  dataAreMagnitudes = true;
  nFunctions = 0;
  nFunctionSets = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
  dataStartOffset = 0;
  debugLevel = 0;
  nCombined = 1;
  nCombined_sqrt = 1.0;
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
      
  if (inputType == WEIGHTS_ARE_SIGMAS) {
  }

  weightValsSet = true;
}


/* ---------------- PUBLIC METHOD: AddPSFVector1D ---------------------- */
// This function needs to be redefined because the base function in ModelObject
// assumes a 2-D PSF.
// Still mostly a stub function at this point!
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
  int  offset = 0;

  // Check parameter values for sanity
  if (! CheckParamVector(nParamsTot, params)) {
    printf("** ModelObject1d::CreateModelImage -- non-finite values detected in parameter vector!\n");
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
      // start of new function set: extract x0 and then skip over them
      x0 = params[offset];
      offset += 1;
// #ifdef DEBUG
//       printf("  Function %d = start of new set; x0 = %g\n",
//               n, x0);
// #endif
    }
    functionObjects[n]->Setup(params, offset, x0);
    offset += paramSizes[n];
  }
  
  // populate modelVector with the model
  for (int i = 0; i < nModelVals; i++) {
    x = modelXValues[i];
    newVal = 0.0;
    for (int n = 0; n < nFunctions; n++)
      newVal += functionObjects[n]->GetValue(x);
    modelVector[i] = newVal;
//    printf("newVal = %g  ", newVal);
  }


  // Do PSF convolution, if requested
  if (doConvolution) {
    psfConvolver->ConvolveProfile(modelVector);
  }
  
  // Convert to magnitudes, if required
  if (dataAreMagnitudes) {
    for (int i = 0; i < nModelVals; i++) {
      modelVector[i] = -2.5 * log10(modelVector[i]);
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
  
  for (int z = 0; z < nDataVals; z++) {
    yResults[z] = weightVector[z] * (dataVector[z] - modelVector[dataStartOffset + z]);
    //printf("weight = %g, data = %g, model = %g\n", weightVector[z], dataVector[z], modelVector[z]);
  }
  //for (int z = 0; z < nDataVals; z++) {
  //  yResults[z] = yWeights[z] * (yVals[z] - GetValue(xVals[z], params));
  //}
}


/* ---------------- PUBLIC METHOD: PrintDescription -------------------- */

void ModelObject1d::PrintDescription( )
{
  printf("ModelObject(1d): %d data values\n", nDataVals);
}


/* ---------------- PUBLIC METHOD: PrintModelParams --------=---------- */
// Basic function which prints to a file a summary of the best-fitting model,
// in form suitable for future use as an input config file.

void ModelObject1d::PrintModelParams( double params[], FILE *output_ptr )
{
  double  x0, paramVal;
  int nParamsThisFunc;
  int  offset = 0;
  string  funcName, paramName;

  for (int n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0 and then skip over it
      x0 = params[offset];
      offset += 1;
      fprintf(output_ptr, "\nX0\t\t%g\n", x0);
    }
    
    // Now print the function and its parameters
    nParamsThisFunc = paramSizes[n];
    funcName = functionObjects[n]->GetShortName();
    fprintf(output_ptr, "FUNCTION %s\n", funcName.c_str());
    for (int i = 0; i < nParamsThisFunc; i++) {
      paramName = GetParameterName(offset + i);
      paramVal = params[offset + i];
      fprintf(output_ptr, "%s\t\t%g\n", paramName.c_str(), paramVal);
    }
    offset += paramSizes[n];
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
}



/* END OF FILE: model_object1d.cpp ------------------------------------- */
