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


/* ---------------- Definitions ---------------------------------------- */


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
  nParamsTot = 0;
}


/* ---------------- PUBLIC METHOD: AddDataVectors --------------------- */

void ModelObject1d::AddDataVectors( int nDataValues, double *xValVector, 
										double *yValVector, bool magnitudeData )
{
  nDataVals = nValidDataVals = nDataValues;
  xValuesVector = xValVector;
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


/* ---------------- PUBLIC METHOD: CreateModelImage -------------------- */

void ModelObject1d::CreateModelImage( double params[] )
{
  double  x, newVal;
  int  offset = 0;

//  printf("CreateModelImage: nFunctions = %d\n", nFunctions);
  // Separate out the individual-component parameters and tell the
  // associated function objects to do setup work.
  // The first component's parameters start at params[0]; the second's
  // start at params[paramSizes[0]], the third at 
  // params[paramSizes[0] + paramSizes[1]], and so forth...
  for (int n = 0; n < nFunctions; n++) {
    //printf("ModelObject1d: component %d setup: ", n);
//    printf("Setting up function object %d with offset = %d; ", n, offset);
//    printf("params = (%g, %g)\n", params[0], params[1]);
    functionObjects[n]->Setup(params, offset);
    offset += paramSizes[n];
  }
  
  // populate modelVector with the model
  for (int i = 0; i < nDataVals; i++) {
    x = xValuesVector[i];
    newVal = 0.0;
    for (int n = 0; n < nFunctions; n++)
      newVal += functionObjects[n]->GetValue(x);
    if (dataAreMagnitudes)
      newVal = -2.5 * log10(newVal);
    modelVector[i] = newVal;
//    printf("newVal = %g  ", newVal);
  }
  modelImageComputed = true;
}


/* ---------------- PUBLIC METHOD: ComputeDeviates --------------------- */

void ModelObject1d::ComputeDeviates( double yResults[], double params[] )
{

  CreateModelImage(params);
  
  for (int z = 0; z < nDataVals; z++) {
    yResults[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
    //printf("weight = %g, data = %g, model = %g\n", weightVector[z], dataVector[z], modelVector[z]);
  }
  //for (int z = 0; z < nDataVals; z++) {
  //  yResults[z] = yWeights[z] * (yVals[z] - GetValue(xVals[z], params));
  //}
}


/* ---------------- PUBLIC METHOD: GetDescription ---------------------- */

void ModelObject1d::GetDescription( )
{
  printf("Model Object(1d): %ld data values\n", nDataVals);
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

ModelObject1d::~ModelObject1d()
{
  if (modelVectorAllocated)
    free(modelVector);
  
  if (nFunctions > 0)
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
}



/* END OF FILE: model_object1d.cpp ------------------------------------- */
