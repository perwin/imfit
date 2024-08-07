// NOTES FOR MULTIMFI-PARAMETERS CLASS
// 
// 
// Stores:
// 	Current image-description parameter values, for all images
// 	[Current image-specific function parameter values, when we get around
// 	to adding these]
// 	Current general model parameter values
// 	Parameter labels?

// The following are double-precision arrays containing image-description parameters.
// In most or all cases, we *ignore* the first element in each, which would
// correspond to the reference image.
//   pixelScales[nImages], rotationAngles[nImages], intensityScales[nImages]
//   imageX0s[nImages], imageY0x[nImages]
// 
// The following are copies of the variables in ModelObject and ModelObjectMultImage
//   nImages = same as nModelObjs in ModelObjectMultImage
//   nModelParameters = same as nFuncParamsTot in ModelObjectMultImage
//   nImageParameters = N_IMAGE_PARAMS*(nImages - 1)
//      -- this is the number of image-description parameters in the externally-provided
//      parameter vector (e.g., generated by the minimization algorithm)
//   functionSetStartFlags[nFunctions] : array of bool = whether a function is
//                                              start of new function set
//   parameterSizes[nFunctions] : vector<int> = how many parameters each function has


// Methods
// 	[] Generate initial 1D parameter array for input to solvers (in main)
// 	[] Take 1D parameter array from solvers and use it to update internal
// 	parameter values
// 	[] Generate image-description parameter values for a ModelObject
// 		E.g., GetImageDescriptionParams(imageDescriptionParams)
// 		
// 	[] Generate general-model-only parameter values for a ModelObject
// 	
// 	[] Print current parameter values to screen or file
// 



/* ------------------------ Include Files (Header Files )--------------- */

#include <stdlib.h>
#include <vector>

#include "definitions_multimage.h"
#include "param_holder.h"
#include "paramvector_processing.h"


/* ---------------- Definitions ---------------------------------------- */



/* ---------------- CONSTRUCTOR ---------------------------------------- */

ParamHolder::ParamHolder( )
{
  arraysAllocated = false;
  perImageModelParamsVect_allocated = false;
  
  nLocalParams = 0;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

ParamHolder::~ParamHolder( )
{
  if (arraysAllocated) {
    free(pixelScales);
    free(rotationAngles);
    free(intensityScales);
    free(imageX0s);
    free(imageY0s);
    free(functionSetStartFlags);
  }
  
  if (perImageModelParamsVect_allocated) {
    for (int i = 0; i < nImages; i++)
      free(perImageModelParameters[i]);
    perImageModelParameters.clear();
    modelParamArray = nullptr;
  }
}



// Tell ParamHolder object details about the model we're working with: number of
// images (ModelObject instances); number of functions, function sets, and
// number of parameters
//    nFuncs = number of functions in *global* model
//    nParamsTot = total number of parameters in ModelObjectMultImage model
//               = global *and* local function parameters *and* image-description parameters
//    nFuncSets = number of function sets in the *global* model
//    paramSizes = vector of # of individual-function parameter values, for each function
//               in *global* model
//    fsetStartFlags = flags for whether functions in *global* model are at start of
//               new function set
void ParamHolder::GeneralSetup( int nImagesTot, int nFuncs, int nParamsTot, int nFuncSets,
								const vector<int> paramSizes, bool fsetStartFlags[],
								const vector<int> localParamCounts )
{
  nImages = nImagesTot;
  nImageParameters = N_IMAGE_PARAMS * (nImages - 1);
  nFunctions = nFuncs;
  nParametersTot = nParamsTot;
  nFunctionSets = nFuncSets;
  
  // figure out how many per-image-function parameters there are
  nLocalParamsForImages = localParamCounts;
  for (int i = 0; i < nImages; i++)
    nLocalParams += nLocalParamsForImages[i];
  nModelParameters = nParametersTot - nImageParameters - nLocalParams;

  pixelScales = (double *)calloc(nImages, sizeof(double));
  rotationAngles = (double *)calloc(nImages, sizeof(double));
  intensityScales = (double *)calloc(nImages, sizeof(double));
  imageX0s = (double *)calloc(nImages, sizeof(double));
  imageY0s = (double *)calloc(nImages, sizeof(double));
  functionSetStartFlags = (bool *)calloc(nFunctions, sizeof(bool));
  arraysAllocated = true;
  
  for (int i = 0; i < nFunctions; i++)
    functionSetStartFlags[i] = fsetStartFlags[i];
  parameterSizes = paramSizes;
  
  // initial setup for reference image [not normally used]
  pixelScales[0] = 1.0;
  rotationAngles[0] = 0.0;
  intensityScales[0] = 1.0;
  imageX0s[0] = 0.0;
  imageY0s[0] = 0.0;
  
  for (int nImage = 0; nImage < nImages; nImage++) {
    int  nParamsThisImage = nModelParameters + nLocalParamsForImages[nImage];
    modelParamArray = (double *)calloc(nParamsThisImage, sizeof(double));
    perImageModelParameters.push_back(modelParamArray);
    nParamsPerModelObj.push_back(nParamsThisImage);
  }
  perImageModelParamsVect_allocated = true;
}


// Pass in the updated global parameter vector; extract image-description
// parameters and calculate per-ModelObject parameter arrays
void ParamHolder::AddNewParameterVector( double params[] )
{
  double  pixScale, rotation, intensityScale;
  double  x0, y0;
  int  status;
  
  // Extract and store image-description parameters
  for (int nImage = 0; nImage < nImages; nImage++) {
    if (nImage == 0) {
      // default values for reference image
      pixelScales[0] = 1.0;
      rotationAngles[0] = 0.0;
      intensityScales[0] = 1.0;
      imageX0s[0] = 0.0;
      imageY0s[0] = 0.0;
    } else {
      status = ExtractImageParams(params, nImage, pixScale, rotation, intensityScale, 
      	 							x0, y0);
      pixelScales[nImage] = pixScale;
      rotationAngles[nImage] = rotation;
      intensityScales[nImage] = intensityScale;
      imageX0s[nImage] = x0;
      imageY0s[nImage] = y0;      
    }
  }
  
  // Generate the per-image parameter vectors we'll use later
  for (int nImage = 0; nImage < nImages; nImage++)
    AssembleParametersForImage(params, nImage, perImageModelParameters[nImage]);
}


/// Returns the current set of image-description parameters for image number nImage,
/// storing them in imageParamsArray
/// Called by external functions (e.g., inside ModelObjectMultImage)
void ParamHolder::GetImageParams( int nImage, double imageParamsArray[] )
{
  imageParamsArray[0] = pixelScales[nImage];
  imageParamsArray[1] = rotationAngles[nImage];
  imageParamsArray[2] = intensityScales[nImage];
}


/// Returns the current set of function parameters usable by a ModelObject instance,
/// with X0,Y0 values adjusted for image number nImage (and including appropriate
/// per-image parameters for this instance, if any), storing them in 
/// paramsForModelObject
/// Called by external functions (e.g., inside ModelObjectMultImage)
void ParamHolder::GetParamsForModelObject( int nImage, double paramsForModelObject[] )
{
  for (int i = 0; i < nParamsPerModelObj[nImage]; i++)
    paramsForModelObject[i] = perImageModelParameters[nImage][i];
}






// NOTE: The following is intended to be a PRIVATE function (not called from outside);
// we declared it public so we could do unit tests on it.
// Python prototype: parameter_munging.AssembleParams
// Given a raw parameter vector (e.g., as generated by an external minimization
// algorithm) which consists of nImageParameters values followed by nModelParameters
// values, this function generates a vector of parameters for the ModelObject
// instance indexed by imageNumber, including appropriate offset/rotated/scaled
// X0,Y0 values for individual function sets.
void ParamHolder::AssembleParametersForImage( double externalInputParamsVect[], int imageNumber, 
											double outputModelParams[] )
{
  double pixScale_im, rot_im, iScale;
  double X0_0_im, Y0_0_im;

  // Copy input *model* parameters (skipping over image-description parameters)
  // Note that model parameters start at externalInputParamsVect[nImageParameters]
  // (i.e., skipping over the initial image-description parameters)
  for (int i = 0; i < nModelParameters; i++)
    outputModelParams[i] = externalInputParamsVect[nImageParameters + i];

  // If this is for 2nd or subsequent image, update X0,Y0 for first function set
  // and work out transformed X0,Y0 values for 2nd and subsequent function sets
  if (imageNumber > 0) {
    ExtractImageParams(externalInputParamsVect, imageNumber, pixScale_im, rot_im, iScale, 
    					X0_0_im, Y0_0_im);
    // start updating copy of main vector with proper X0,Y0 for this image
    outputModelParams[0] = X0_0_im;
    outputModelParams[1] = Y0_0_im;
    // now compute offsets and second and subsequent function-set X0,Y0 values
    if (nFunctionSets > 1) {
      double X0_0_ref = externalInputParamsVect[nImageParameters];
      double Y0_0_ref = externalInputParamsVect[nImageParameters + 1];
      // skip the first function (X0,Y0 + parameterSizes[0])
      int offset = parameterSizes[0] + 2;
      for (int n = 1; n < nFunctions; n++) {
        if (functionSetStartFlags[n] == true) {
          // new function set, so we calculate and store offset/rotated/scaled
          // X0,Y0 for this function set in this image
          double X0_n_ref = externalInputParamsVect[nImageParameters + offset];
          double Y0_n_ref = externalInputParamsVect[nImageParameters + offset + 1];
          double X0_n_im, Y0_n_im;
          std::tie(X0_n_im, Y0_n_im, std::ignore) = CalculateOffset_X0Y0(X0_0_ref, Y0_0_ref, 
          											X0_n_ref, Y0_n_ref, X0_0_im, Y0_0_im,
          											pixScale_im, rot_im);
          outputModelParams[offset] = X0_n_im;
          outputModelParams[offset + 1] = Y0_n_im;
          offset += 2;
        }
        offset += parameterSizes[n];
      }
    }
  }

  // Extract and store of local/per-image-function parameters (if they exist for this image)
  int  nLocalParams_thisImage = nLocalParamsForImages[imageNumber];
  if (nLocalParams_thisImage > 0) {
    // figure out offset from end of global-model parameters due to any previous
    // local/per-image-function parameters
    int  start = nImageParameters + nModelParameters;
    for (int nIm = 0; nIm < imageNumber; nIm++)
      start += nLocalParamsForImages[nIm];
    for (int i = 0; i < nLocalParams_thisImage; i++)
      outputModelParams[nModelParameters + i] = externalInputParamsVect[start + i];
  }
  
}





// How variables are defined/set in ModelObjectMultImage
// 	nParamsTot -- built up by repeated calls to DefineFunctionSets *and* AddModelObject
// 	nFunctionParams -- built up by repeated calls to AddFunction
// 	paramSizes -- built up by repeated calls to AddFunction
// 	nFunctions -- built up by repeated calls to AddFunction
// 	nModelObjects -- built up by repeated calls to AddModelObject
