/* FILE: model_object_multimage.cpp ------------------------------------ */
/* 
 * This is intended to be the main class for the "model" object (i.e., image data + 
 * fitting functions); it can also serve as a base class for derived versions thereof 
 *(e.g., fitting 1D models to profiles).
 * 
 *
 *
 * March 2017: Created as variation of model_object.cpp
 */

// Copyright 2010--2022 by Peter Erwin.
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

// logging
#ifdef USE_LOGGING
#include "loguru/loguru.hpp"
#endif

#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object_multimage.h"
//#include "oversampled_region.h"
#include "add_functions.h"
#include "mp_enorm.h"
#include "param_struct.h"
#include "utilities_pub.h"
#include "paramvector_processing.h"
#include "param_holder.h"
#include "definitions_multimage.h"


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


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModelObjectMultImage::ModelObjectMultImage( )
{
  // Extra stuff not done by ModelObject constructor
  nModelObjects = nGlobalFuncParams = 0;
  paramsForModelObjects_allocated = false;
  imageDescriptionParams_allocated = false;
  parameterHolderSet = false;
}





/* ---------------- DESTRUCTOR ----------------------------------------- */

// Note that we have to turn various bool variables off and set nFunctions = 0,
// else we have problems when the *base* class (ModelObject) destructor is called
// (which happens automatically just after *this* destructor is called) -- we
// can end up trying to free vectors that have already been freed, because the
// associated bool variables are still = true...
ModelObjectMultImage::~ModelObjectMultImage()
{

  // NEW
  // deallocate individual arrays of per-ModelObject parameters
  if (paramsForModelObjects_allocated) {
    for (int i = 0; i < nModelObjects; i++)
      free(paramsForModelObjects[i]);
    paramsForModelObjects.clear();
    paramsForModelObjects_allocated = false;
  }

  if (imageDescriptionParams_allocated) {
    for (int i = 0; i < nModelObjects; i++)
      free(imageDescriptionParams[i]);
    imageDescriptionParams.clear();
    imageDescriptionParams_allocated = false;
  }

  // deallocate individual ModelObject instances in modelObjectsVect
  // (since these were originally created with "new", we have to deallocate with "delete")
  if (modelObjectsVect.size() > 0) {
    for (int i = 0; i < (int)modelObjectsVect.size(); i++)
      delete modelObjectsVect[i];
    modelObjectsVect.clear();
    nModelObjects = 0;
  }
  
  // Now the base class (ModelObject) destructor is automatically called...
}


int ModelObjectMultImage::SetupGeneralOptions( shared_ptr<OptionsBase> options )
{
  if (options->maxThreadsSet)
    SetMaxThreads(options->maxThreads);
  SetDebugLevel(options->debugLevel);
  
  return 0;
}


void ModelObjectMultImage::AddModelObject( ModelObject* newModelObj )
{
  double *currentModelParamsVect;
  double *imageDescription;
  int nCurrentModelParams;
  
  if (nModelObjects > 0)
    nParamsTot += N_IMAGE_PARAMS;
  modelObjectsVect.push_back(newModelObj);
  nModelObjects++;
  
  // Other setup stuff ...
  
  // NOTE: we assume that the input ModelObject instance *already* has its set
  // of FunctionObjects, so we can query it for the number of model parameters
  nCurrentModelParams = newModelObj->GetNParams();
  currentModelParamsVect = (double *) calloc((size_t)nCurrentModelParams, sizeof(double));
  if (currentModelParamsVect == nullptr) {
    fprintf(stderr, "*** ERROR: Unable to allocate memory for input ModelObject parameter vector!\n");
    fprintf(stderr, "    (Requested vector size was %d)\n", nCurrentModelParams);
#ifdef USE_LOGGING
    LOG_F(ERROR, "Unable to allocate memory for input ModelObject parameter vector!");
#endif
    //return -1;
  }
  paramsForModelObjects.push_back(currentModelParamsVect);
  nParamsForModelObjects.push_back(nCurrentModelParams);
  paramsForModelObjects_allocated = true;
  
  imageDescription = (double *) calloc((size_t)N_IMAGE_DESCRIPTION_PARAMS, sizeof(double));
  imageDescriptionParams_allocated = true;  
  // default values
  imageDescription[0] = 1.0;  // pixelScale
  imageDescription[1] = 0.0;  // rotation
  imageDescription[2] = 1.0;  // intensityScale
  imageDescriptionParams.push_back(imageDescription);
  
  // add image-description parameter names
  if (nModelObjects > 1) {
    parameterLabels.push_back("PIXEL_SCALE");
    parameterLabels.push_back("ROTATION");
    parameterLabels.push_back("FLUX_SCALE");
    parameterLabels.push_back("X0");
    parameterLabels.push_back("Y0");
  }
    
  nDataVals += newModelObj->GetNDataValues();
}


// This is currently called from main() -- via add_function.cpp's AddFunction() -- *after* 
// all the individual ModelObject instances have had their AddFunction methods called
// When this is finished, ModelObjectMultImage will know about *global* functions, but
// not any of the per-image functions
int ModelObjectMultImage::AddFunction( FunctionObject *newFunctionObj_ptr, bool isGlobalFunc )
{
  int  nNewParams;

  if (nModelObjects < 1) {
    fprintf(stderr, "WARNING: ModelObjectMultImage::AddFunction called before any");
    fprintf(stderr, " ModelObject instances were added!\n");
    return -1;
  }
  functionObjects.push_back(newFunctionObj_ptr);
  nFunctions += 1;
  nNewParams = newFunctionObj_ptr->GetNParams();
  paramSizes.push_back(nNewParams);
  nFunctionParams += nNewParams;
  
  return 0;
}


// We redefine this to give an error message, since it *shouldn't* be called
// for this subclass!
void ModelObjectMultImage::CreateModelImage( double params[] )
{
  fprintf(stderr, "ERROR: CreateModelImage should not be called on ModelObjectMultImage instance!\n");
#ifdef USE_LOGGING
  LOG_F(ERROR, "CreateModelImage was erroneously called on ModelObjectMultImage instance!");
#endif
}



int ModelObjectMultImage::GetNImages( )
{
   return nModelObjects;
}



long ModelObjectMultImage::GetNValidPixels( )
{
  long  nPix = 0;
  
  for (int i = 0; i < nModelObjects; i++)
    nPix += modelObjectsVect[i]->GetNValidPixels();
  return nPix;
}



std::tuple<int, int> ModelObjectMultImage::GetImageDimensions( int imageNumber )
{
  int nColumns, nRows;
  modelObjectsVect[imageNumber]->GetDataImageDimensions(&nColumns, &nRows);
  return std::make_tuple(nColumns, nRows);
}



void ModelObjectMultImage::SetupParamHolder( )
{
  parameterHolder.GeneralSetup(nModelObjects, nFunctions, nParamsTot, nFunctionSets, 
    							paramSizes, fsetStartFlags, nPerImageParamsForModelObjects);
  parameterHolderSet = true;
}



/* ---------------- PUBLIC METHOD: FinalModelSetup -------------------- */
/// Do any final setup involving determination of parameter numbers (including
/// accounting for global vs per-image functions in multimfit mode)
//
// We assume that the following has already occurred:
//   1. all functions have been added via repeated calls to AddFunction()
//   2. DefineFunctionSets(functionSetIndices) has been called;
//   3. PopulateParameterNames() has been called
int ModelObjectMultImage::FinalModelSetup( )
{
  int  nParamsThisImage, nExtra;

  nGlobalFuncParams = nFunctionParams + 2*nFunctionSets;

  // Compute nPerImageParams, for each ModelObject and running total
  for (int i = 0; i < nModelObjects; i++) {
    nParamsThisImage = nParamsForModelObjects[i];
    nExtra = nParamsThisImage - nGlobalFuncParams;
    nPerImageParamsForModelObjects.push_back(nExtra);
    if (nExtra > 0) {
      nPerImageParams += nExtra;
    }
  }
  nParamsTot += nPerImageParams;
  
  // This must be done *after* we've computed the final values of nParamsTot, etc.
  SetupParamHolder();
  
  return 0;
}


int ModelObjectMultImage::FinalSetupForFitting( )
{
  int  status;
  
  for (int i = 0; i < nModelObjects; i++) {
      status = modelObjectsVect[i]->FinalSetupForFitting();
      if (status < 0)
        return status;
  }
  
  if (! parameterHolderSet)
    SetupParamHolder();
  
  return 0;
}

// parameter vector should be of the form:
//    [ pixScale_1,rotation_1,intensityScale_1,X0_1,Y0_1, ..., 
//       <params for first ModelObject> ]
// where <params for first ModelObject> = X0_0,Y0_0, <first function params>, ...
void ModelObjectMultImage::CreateAllModelImages( double params[] )
{
  string textLine;
  
  if (! parameterHolderSet)
    SetupParamHolder();

#ifdef USE_LOGGING
  LOG_F(2, "In CreateAllModelImages:");
  textLine = "   params[] = ";
  for (int i = 0; i < nParamsTot; i++)
    textLine += PrintToString("\t%.4f", params[i]);
  LOG_F(2, "%s\n", textLine.c_str());
#endif

  parameterHolder.AddNewParameterVector(params);

  for (int i = 0; i < nModelObjects; i++) {
    parameterHolder.GetImageParams(i, imageDescriptionParams[i]);
    modelObjectsVect[i]->SetImageParameters(imageDescriptionParams[i]);
    parameterHolder.GetParamsForModelObject(i, paramsForModelObjects[i]);

#ifdef USE_LOGGING
    LOG_F(2, "  image i = %d: pixScale, rot, iScale = ", i);
    textLine = "   ";
    for (int j = 0; j < 3; j++)
      textLine += "   " + PrintToString("%.1f", imageDescriptionParams[i][j]);
    LOG_F(2, "%s", textLine.c_str());
    textLine = "  params for ModelObj:";
    for (int j = 0; j < nGlobalFuncParams + nPerImageParamsForModelObjects[i]; j++)
      textLine += PrintToString("\t%.4f", paramsForModelObjects[i][j]);
    LOG_F(2, "%s", textLine.c_str());
#endif

    modelObjectsVect[i]->CreateModelImage(paramsForModelObjects[i]);
  }
}



double ModelObjectMultImage::GetFitStatistic( double params[] )
{
  double  cumulativeFitStatistic;
  double  fitstat;
  
  if (! parameterHolderSet)
    SetupParamHolder();

  cumulativeFitStatistic = 0.0;
  
#ifdef USE_LOGGING
  LOG_F(2, "In GetFitStatistic:");
  string  textLine = "   params[] = ";
  for (int i = 0; i < nParamsTot; i++)
    textLine += PrintToString("\t%.4f", params[i]);
  LOG_F(2, "%s\n", textLine.c_str());
#endif

  parameterHolder.AddNewParameterVector(params);
  
  for (int i = 0; i < nModelObjects; i++) {
    parameterHolder.GetImageParams(i, imageDescriptionParams[i]);
    modelObjectsVect[i]->SetImageParameters(imageDescriptionParams[i]);
    parameterHolder.GetParamsForModelObject(i, paramsForModelObjects[i]);

#ifdef USE_LOGGING
    LOG_F(2, "image i = %d, pixScale, rot, iScale = %.1f, %.1f, %.1f", i,
    	imageDescriptionParams[i][0], imageDescriptionParams[i][1], imageDescriptionParams[i][2]);
    textLine = "  params for ModelObj:";
    for (int j = 0; j < nGlobalFuncParams + nPerImageParamsForModelObjects[i]; j++)
      textLine += PrintToString("\t%.4f", paramsForModelObjects[i][j]);
    LOG_F(2, "%s", textLine.c_str());
#endif

    if (useCashStatistic)
      cumulativeFitStatistic += modelObjectsVect[i]->CashStatistic(paramsForModelObjects[i]);
    else {
      fitstat = modelObjectsVect[i]->ChiSquared(paramsForModelObjects[i]);
      cumulativeFitStatistic += fitstat;
#ifdef USE_LOGGING
      LOG_F(2, "   fitstat = %f", fitstat);
#endif
    }
  }

  return cumulativeFitStatistic;
}


void ModelObjectMultImage::ComputeDeviates( double yResults[], double params[] )
{
  long  startIndex = 0;

#ifdef USE_LOGGING
  LOG_F(2, "In ComputeDeviates:");
#endif

  if (! parameterHolderSet)
    SetupParamHolder();

  parameterHolder.AddNewParameterVector(params);

  for (int i = 0; i < nModelObjects; i++) {
    parameterHolder.GetImageParams(i, imageDescriptionParams[i]);
    modelObjectsVect[i]->SetImageParameters(imageDescriptionParams[i]);
    parameterHolder.GetParamsForModelObject(i, paramsForModelObjects[i]);

#ifdef USE_LOGGING
    string  textLine = PrintToString("   image i = %d: pixScale, rot, iScale =", i);
    for (int j = 0; j < 3; j++)
      textLine += "   " + PrintToString("%.1f", imageDescriptionParams[i][j]);
    LOG_F(2, "%s", textLine.c_str());
    textLine = "  params for ModelObj:";
    for (int j = 0; j < nGlobalFuncParams + nPerImageParamsForModelObjects[i]; j++)
      textLine += PrintToString("\t%.4f", paramsForModelObjects[i][j]);
    LOG_F(2, "%s", textLine.c_str());
#endif

    modelObjectsVect[i]->ComputeDeviates(yResults + startIndex, 
    									paramsForModelObjects[i]);
    startIndex += modelObjectsVect[i]->GetNDataValues();
  }
}


double * ModelObjectMultImage::GetModelImageVector( int imageNumber )
{
  return modelObjectsVect[imageNumber]->GetModelImageVector();
}


/* ---------------- PUBLIC METHOD: PrintDescription ------------------- */

void ModelObjectMultImage::PrintDescription( )
{
  // Don't test for verbose level, since we assume user only calls this method
  // if they *want* printed output
  printf("ModelObjectMultImage: %ld total data values\n", nDataVals);
  if (nModelObjects > 0)
    for (int i = 0; i < nModelObjects; i++) {
      printf("   ModelObject %d: ", i + 1);
      modelObjectsVect[i]->PrintDescription();
    }
}



/* ---------------- PUBLIC METHOD: PrintModelParamsToStrings ---------- */
/// Basic function which prints to a vector of strings a summary of the
/// best-fitting model, in form suitable for future use as an input config file.
///
/// If errs != NULL, then +/- errors are printed as well
///
/// If prefix != NULL, then the specified character (e.g., '#') is prepended to
/// each output line.
///
/// Note that this version (in ModelObjectMultImage) prints the image-description
/// parameters and the function parameters in their default form (i.e., scaled for the
/// reference image). To print function parameters scaled/transformed for *other*
/// images than the reference image, use GetParameterStringsForOneImage
///
/// Order of output:
///    1. "Main model parameters (for reference image"
///    2. Model parameters [*no* image-description parameters], with errors if provided
///    3. Loop over additional images:
///       A. "Image-description parameters for image %d"
///       B. Image-description parameters for this image, with errors if provided
int ModelObjectMultImage::PrintModelParamsToStrings( vector<string> &stringVector, double params[], 
									double errs[], const char *prefix, bool printLimits )
{
//  double  x0, y0, paramVal;
  double  pixScale, rot, iScale, imageX0, imageY0;
  double  err_pixScale, err_rot, err_iScale, err_imageX0, err_imageY0;
  int  indexOffset = 0;
  int  offsetX0 = 0;
  int  offsetY0 = 0;
  string  funcName, paramName, newLine, imageName;

  if ((printLimits) && (parameterInfoVect.size() == 0)) {
    fprintf(stderr, "** ERROR: ModelObjectMultImage::PrintModelParamsToStrings -- printing of parameter limits\n");
    fprintf(stderr, "was requested, but parameterInfoVect is empty!\n");
    return -1;
  }

  imageName = modelObjectsVect[0]->GetDataFilename();
  newLine = PrintToString("\n%s Main model parameters (for reference image = image 1 (%s))", 
  							"#", imageName.c_str());  
  stringVector.push_back(newLine);

  GetParameterStringsForOneImage(stringVector, params, 0, errs);
  // TODO: Get local, per-image functions and parameters for this image, if they exist
  GetPerImageFuncStrings(stringVector, params, 0, errs);

  
  // Now print image-description parameters for second and subsequent images
  for (int i = 1; i < nModelObjects; i++) {
    imageName = modelObjectsVect[i]->GetDataFilename();
    newLine = PrintToString("\n%s Image-description parameters for image %d (%s)\n", "#", 
    						i + 1, imageName.c_str());
    stringVector.push_back(newLine);
    ExtractImageParams(params, i, pixScale, rot, iScale, imageX0, imageY0);
    // correct imageX0,imageY0 for image-section-based offsets
    std::tie(offsetX0, offsetY0) = modelObjectsVect[i]->GetImageOffsets();
    imageX0 += offsetX0;
    imageY0 += offsetY0;
    if (errs != NULL) {
      ExtractImageParams(errs, i, err_pixScale, err_rot, err_iScale, err_imageX0, err_imageY0);
      stringVector.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, prefix, "PIXEL_SCALE", 
      							pixScale, err_pixScale) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, prefix, "ROTATION", 
      							rot, err_rot) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, prefix, "FLUX_SCALE", 
      							iScale, err_iScale) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, prefix, "X0", 
      							imageX0, err_imageX0) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, prefix, "Y0", 
      							imageY0, err_imageY0) + "\n");
      GetPerImageFuncStrings(stringVector, params, i, errs);
    }
    else {
      stringVector.push_back(PrintToString(PARAM_FORMAT, prefix, "PIXEL_SCALE", pixScale) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT, prefix, "ROTATION", rot) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT, prefix, "FLUX_SCALE", iScale) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT, prefix, "X0", imageX0) + "\n");
      stringVector.push_back(PrintToString(PARAM_FORMAT, prefix, "Y0", imageY0) + "\n");
      GetPerImageFuncStrings(stringVector, params, i);
    }
	indexOffset += N_IMAGE_PARAMS;
  }
  // TODO: Get local, per-image functions and parameters for this image, if they exist

  return 0;
}



/* ---------------- PUBLIC METHOD: GetParameterStringsForOneImage ----- */
/// Basic function which generates a list of strings specifying the function
/// parameters for a single-image output file, corresponding to imageNumber.
void ModelObjectMultImage::GetParameterStringsForOneImage( vector<string> &strings, 
									double params[], int imageNumber, double errors[] )
{
  double *paramErrors = NULL;
  
  if (! parameterHolderSet)
    SetupParamHolder();

  // put full parameter vector into parameterHolder
  parameterHolder.AddNewParameterVector(params);
  // Get the set of function parameters, with X0,Y0 adjusted for image #imageNumber
  // along with per-image function parameters (if any)
  parameterHolder.GetParamsForModelObject(imageNumber, paramsForModelObjects[imageNumber]);

  // Get image-description parameters for image #imageNumber
  parameterHolder.GetImageParams(imageNumber, imageDescriptionParams[imageNumber]);
  // tell ModelObject #imageNumber to pass image-description parameters to its
  // internal FunctionObjects, so they can adjust intensities, sizes, etc.
  modelObjectsVect[imageNumber]->SetImageParameters(imageDescriptionParams[imageNumber]);
  
  // If errors were supplied, compute offset into the error vector so we skip over
  // the image-description-parameter errors and get just the function-parameter errors, 
  // not the image-description-parameter errors
  if (errors != NULL) {
    int nImageParams = (nModelObjects - 1)*N_IMAGE_PARAMS;
    paramErrors = errors + nImageParams;
  }

  // Finally, tell ModelObject #imageNumber to print the (adjusted) function parameter
  // values to strings
  string prefix = "";
  modelObjectsVect[imageNumber]->PrintModelParamsToStrings(strings, 
  										paramsForModelObjects[imageNumber], paramErrors, 
  										"", false);
}



/* ---------------- PUBLIC METHOD: GetPerImageFuncStrings ------------- */
/// Basic function which adds to a vector of strings specifying the function
/// parameters for the *per-image functions* of a single-image ModelObject, 
/// instance, corresponding to imageNumber.
/// If that particular ModelObject instance has no per-image functions, then
/// no extra strings are added.
void ModelObjectMultImage::GetPerImageFuncStrings( vector<string> &strings, 
									double params[], int imageNumber, double errors[] )
{
  ;
}



/* ---------------- PUBLIC METHOD: GetAllImageDescriptionParameters --- */
/// Given a set of global parameter values (e.g., as returned by minimization process),
/// Returns just the image-description "prefix", storing it in outputImageParams
/// Used by SaveImageInfoParameters in print_results_multi.cpp
void ModelObjectMultImage::GetAllImageDescriptionParameters( double params[], 
															double outputImageParams[] )
{
  for (int i = 0; i < (nModelObjects - 1)*N_IMAGE_PARAMS; i++)
    outputImageParams[i] = params[i];
}



/* ---------------- PUBLIC METHOD: GetOffsetsForImage ----------------- */
/// Returns the x and y offsets for the specified data image.
std::tuple<int, int> ModelObjectMultImage::GetOffsetsForImage( int imageNumber )
{
  return modelObjectsVect[imageNumber]->GetImageOffsets();
}



/* ---------------- PUBLIC METHOD: GetDataFilename -------------------- */
/// Returns the x and y offsets for the specified data image.
string ModelObjectMultImage::GetDataFilename( int imageNumber )
{
  return modelObjectsVect[imageNumber]->GetDataFilename();
}



/* ---------------- PUBLIC METHOD: PrintModelParamsHorizontalString --- */
/// Like PrintModelParamsToString, but prints parameter values all in one line to
/// a string (*without* parameter names or limits or errors), which is returned.
/// Meant to be used in printing results of bootstrap resampling (imfit) or MCMC
/// chains (imfit-mcmc)
string ModelObjectMultImage::PrintModelParamsHorizontalString( const double params[], 
																const string& separator )
{
  double  x0, y0, paramVal;
  int  image1Offset_x0, image1Offset_y0;
  int  nParamsThisFunc, k, indexOffset, nImageParams;
  string  outputString = "";

  // pixel offsets for reference image
  std::tie(image1Offset_x0, image1Offset_y0) = modelObjectsVect[0]->GetImageOffsets();
  // number of image-description parameters
  nImageParams = (nModelObjects - 1)*N_IMAGE_PARAMS;
  
  // Line should start with parameter value (no initial separator)
  outputString = PrintToString("%#.10g", params[0]);
  // now we precede each subsequent parameter value with separator
  for (int i = 1; i < nImageParams; i++)
    outputString += PrintToString("%s%#.10g", separator.c_str(), params[i]);

  indexOffset = nImageParams;
  for (int n = 0; n < nFunctions; n++) {
    if (fsetStartFlags[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      k = indexOffset;
      x0 = params[k] + image1Offset_x0;
      y0 = params[k + 1] + image1Offset_y0;
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


// FIXME: Update this to generate base model + image params headers
// i = 0 up to (nModelObjects - 1)*N_IMAGE_PARAMS ==> image parameters
// i = (nModelObjects - 1)*N_IMAGE_PARAMS up to nParamsTot ==> model parameters
//
/* ---------------- PUBLIC METHOD: GetParamHeader ---------------------- */
/// Prints all function and parameter names in order all on one line; e.g., for use as 
/// header in bootstrap-parameters output file.
string ModelObjectMultImage::GetParamHeader( )
{
  int  imageNum;
  string  headerLine, baseModelParamHeader;

  // construct image-parameters part of header
  headerLine = "# ";
  for (int n = 1; n < nModelObjects; n++) {
    imageNum = n + 1;
    headerLine += PrintToString("PIXEL_SCALE_IM_%d\t", imageNum);
    headerLine += PrintToString("IMAGE_PA_IM_%d\t", imageNum);
    headerLine += PrintToString("FLUX_SCALE_IM_%d\t", imageNum);
    headerLine += PrintToString("X0_IM_%d\t", imageNum);
    headerLine += PrintToString("Y0_IM_%d\t", imageNum);
  }
  
  // get model-parameters part of header, chop off leading "# "
  baseModelParamHeader = modelObjectsVect[0]->GetParamHeader();
  headerLine += baseModelParamHeader.substr(2, string::npos);
  
  return headerLine;
}


/* ---------------- PUBLIC METHOD: UseBootstrap ------------------------ */
/// Tells ModelObjectMultImage object that from now on we'll operate in bootstrap
/// resampling mode. 
/// This (overridden) version simply calls the UseBootstrap() method on each of the 
/// ModelObject instances in modelObjectsVect.
/// Returns the "cumulative" status from the individual UseBootstrap() calls, which 
/// will be -1 if any of the memory allocation for the bootstrap-indices vectors failed.
int ModelObjectMultImage::UseBootstrap( )
{
  int  status = 0;
  int  cumulativeStatus = 0;
  
  for (int i = 0; i < nModelObjects; i++) {
    status = modelObjectsVect[i]->UseBootstrap();
    if (status < 0)
      cumulativeStatus = -1;
  }
  return cumulativeStatus;
}


/* ---------------- PUBLIC METHOD: MakeBootstrapSample ----------------- */
/// Generate a new bootstrap resampling of the data (more precisely, this generate a
/// bootstrap resampling of the data *indices*).
/// This (overridden) version simply calls the MakeBootstrapSample() method on each 
/// of the ModelObject instances in modelObjectsVect.
/// Returns the "cumulative" status from the individual MakeBootstrapSample() calls, 
/// which will be -1 if any of the memory allocation for the bootstrap-indices 
/// vectors failed.
int ModelObjectMultImage::MakeBootstrapSample( )
{
  int  status = 0;
  int  cumulativeStatus = 0;
  
  for (int i = 0; i < nModelObjects; i++) {
    status = modelObjectsVect[i]->MakeBootstrapSample();
    if (status < 0)
      cumulativeStatus = -1;
  }
  return cumulativeStatus;
}




/* END OF FILE: model_object_multimage.cpp ----------------------------- */
