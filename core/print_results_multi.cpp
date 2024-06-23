/** @file
 * \brief Utility functions for interpreting and printing results from fits by multimfit.
 */

/* FILE: print_results_multi.cpp ----------------------------------- */

// Copyright 2017-2024 by Peter Erwin.
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


#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

#include "definitions.h"
#include "definitions_multimage.h"
#include "print_results_multi.h"
#include "model_object_multimage.h"
#include "utilities_pub.h"
#include "solver_results.h"



/* ------------------- Function Prototypes ----------------------------- */
void PrintSingleImageInfoToStrings( const ImageInfo& imInfo, vector<string>& strings,
									bool refImageFlag, double *imageDescParams, int offset,
									double *imageDescParams_errs );


/* ------------------------ Global Variables --------------------------- */
#define  FILE_OPEN_ERR_STRING "\n   Couldn't open file \"%s\"\n\n"





void SaveMultImageParameters( double *params, ModelObjectMultImage *multModel, 
								string& outputFilenameRoot, vector<string>& outputHeader )
{
  FILE  *file_ptr;
  vector<string> stringsForFile;
  string  newFilename, dataFilename;

  // NOTE: simplest for now to write errors *only* for the first (reference) image
  for (int i = 0; i < multModel->GetNImages(); i++) {
    stringsForFile.clear();
    if (i == 0) {
      newFilename = PrintToString("%s_refimage.dat", outputFilenameRoot.c_str(), i + 1);
      printf("\tSaving %s (for reference image)...\n", newFilename.c_str());
    }
    else {
      newFilename = PrintToString("%s_image%d.dat", outputFilenameRoot.c_str(), i + 1);
      printf("\tSaving %s...\n", newFilename.c_str());
    }
    file_ptr = fopen(newFilename.c_str(), "w");
    for (auto line: outputHeader)
      fprintf(file_ptr, "%s\n", line.c_str());
    dataFilename = multModel->GetDataFilename(i);
    fprintf(file_ptr, "\n# MODEL PARAMETERS FOR DATA IMAGE: %s\n", dataFilename.c_str());
    multModel->GetParameterStringsForOneImage(stringsForFile, params, i);
    for (auto line: stringsForFile)
      fprintf(file_ptr, "%s", line.c_str());
    fclose(file_ptr);
  }
}



void SaveImageInfoParameters( double *params, ModelObjectMultImage *multModel, 
								vector<ImageInfo>& imageInfoVect, 
								SolverResults& solverResults, string& outputFilename,
								vector<string>& outputHeader )
{
  int  nImages = multModel->GetNImages();
  FILE  *file_ptr;
  string  simpleHeaderLine;
  vector<string>  stringsForFile;
  bool  isReferenceImage;
  bool  errorsPresent = false;
  int  whichSolver;
  double  *imageDescriptionParams = (double *)calloc((size_t)nImages*N_IMAGE_PARAMS, 
  														sizeof(double));
  double  *imageDescriptionParams_errs = NULL;
  double  *paramErrs = NULL;

  multModel->GetAllImageDescriptionParameters(params, imageDescriptionParams);

  whichSolver = solverResults.GetSolverType();
  if (whichSolver == MPFIT_SOLVER) {
    errorsPresent = true;
    paramErrs = (double *)calloc(multModel->GetNParams(), sizeof(double));
    solverResults.GetErrors(paramErrs);
    imageDescriptionParams_errs = (double *)calloc((size_t)nImages*N_IMAGE_PARAMS, 
  														sizeof(double));
    multModel->GetAllImageDescriptionParameters(paramErrs, imageDescriptionParams_errs);
  }

  
  for (int i = 0; i < nImages; i++) {
    if (i == 0) {
      isReferenceImage = true;
      stringsForFile.push_back("# Reference image\n");
    }
    else {
      isReferenceImage = false;
      if (i == 1)
        stringsForFile.push_back("# Additional images\n");
    }
    
    // offset index = (i - 1)*N_IMAGE_PARAMS): skip over previous single-image
    // parameters in units of N_IMAGE_PARAMS; first set of params is for *second*
    // image (i = 1)
    int  offset = (i - 1)*N_IMAGE_PARAMS;
    PrintSingleImageInfoToStrings(imageInfoVect[i], stringsForFile, isReferenceImage,
    							imageDescriptionParams, offset, imageDescriptionParams_errs);
  }
  
  if (stringsForFile.size() > 0) {
    file_ptr = fopen(outputFilename.c_str(), "w");
    for (auto line: outputHeader)
      fprintf(file_ptr, "%s\n", line.c_str());
    fprintf(file_ptr, "\n# IMAGE PARAMETERS\n\n");
    for (auto line: stringsForFile)
      fprintf(file_ptr, "%s", line.c_str());
    fclose(file_ptr);
  }
  
  free(imageDescriptionParams);
  if (errorsPresent) {
    free(paramErrs);
    free(imageDescriptionParams_errs);
  }
}



/// Utility function which prints a block of text with image-info data for a single
/// image, using the data contained in the ImageInfo object; output is appended as
/// sequence of lines to the in/out parameter strings.
void PrintSingleImageInfoToStrings( const ImageInfo& imInfo, vector<string>& strings,
									bool refImageFlag, double *imageDescParams,
									int offset, double *imageDescParams_errs )
{
  strings.push_back("IMAGE_START\n");
  // data image and auxiliary info
  if (imInfo.dataImage_present)
    strings.push_back(PrintToString("DATA    %s\n", imInfo.dataImageFileName.c_str()));
  if (imInfo.originalSky != 0.0)
    strings.push_back(PrintToString("ORIGINAL_SKY\t\t%g\n", imInfo.originalSky));
  if (imInfo.gain != 1.0)
    strings.push_back(PrintToString("GAIN\t\t%g\n", imInfo.gain));
  if (imInfo.readNoise != 0.0)
    strings.push_back(PrintToString("READNOISE\t\t%g\n", imInfo.readNoise));
  if (imInfo.nCombined != 1)
    strings.push_back(PrintToString("NCOMBINED\t\t%g\n", imInfo.nCombined));
  // mask, error, and PSF images
  if (imInfo.maskImage_present)
    strings.push_back(PrintToString("MASK    %s\n", imInfo.maskImageFileName.c_str()));
  if (imInfo.errorImage_present)
    strings.push_back(PrintToString("ERROR    %s\n", imInfo.errorImageFileName.c_str()));
  if (imInfo.psfImage_present)
    strings.push_back(PrintToString("PSF    %s\n", imInfo.psfImageFileName.c_str()));

  // oversampled PSF images and info
  if (imInfo.psfOversampledImage_present) {
    int  nOversampledImages = imInfo.nOversampledFileNames;
    int  nOversampledRegions = imInfo.nOversampleRegions;
    for (int i = 0; i < nOversampledRegions; i++) {
      if (i < nOversampledImages) {
        strings.push_back(PrintToString("OVERSAMPLED_PSF         %s\n", 
        			imInfo.psfOversampledFileNames[i].c_str()));
        strings.push_back(PrintToString("OVERSAMPLE_SCALE        %d\n", 
        			imInfo.psfOversamplingScales[i]));
      }
      strings.push_back(PrintToString("OVERSAMPLED_REGION      %s\n", 
        			imInfo.psfOversampleRegions[i].c_str()));
    }
  }
  
  if (! refImageFlag) {
    // save pixel_scale, image_pa, flux_scale, and initial X0,Y0
    if (imageDescParams_errs == NULL) {
      // no errors on these parameters from the fitting process
      // note that we have to add line-feed at end of each PARAM_FORMAT-based string
      strings.push_back(PrintToString(PARAM_FORMAT, "","PIXEL_SCALE", imageDescParams[0 + offset]) + "\n");
      strings.push_back(PrintToString(PARAM_FORMAT, "","IMAGE_PA", imageDescParams[1 + offset]) + "\n");
      strings.push_back(PrintToString(PARAM_FORMAT, "","FLUX_SCALE", imageDescParams[2 + offset]) + "\n");
      strings.push_back(PrintToString(XY_FORMAT, "","X0", imageDescParams[3 + offset]));
      strings.push_back(PrintToString(XY_FORMAT, "","Y0", imageDescParams[4 + offset]));
    }
    else {
      // note that we have to add line-feed at end of each PARAM_FORMAT-based string
      strings.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, "","PIXEL_SCALE", imageDescParams[0 + offset], 
      									imageDescParams_errs[0 + offset]) + "\n");
      strings.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, "","IMAGE_PA", imageDescParams[1 + offset], 
      									imageDescParams_errs[1 + offset]) + "\n");
      strings.push_back(PrintToString(PARAM_FORMAT_WITH_ERRS, "","FLUX_SCALE", imageDescParams[2 + offset], 
      									imageDescParams_errs[2 + offset]) + "\n");
      strings.push_back(PrintToString(XY_FORMAT_WITH_ERRS, "","X0", imageDescParams[3 + offset], 
      									imageDescParams_errs[3 + offset]));
      strings.push_back(PrintToString(XY_FORMAT_WITH_ERRS, "","Y0", imageDescParams[4 + offset], 
      									imageDescParams_errs[4 + offset]));
    }
  }
// PIXEL_SCALE             1.0			fixed   # pixels in this image are PIXEL_SCALE * reference image
// IMAGE_PA                0.0			fixed
// FLUX_SCALE              1.0		 	fixed
// X0                      54			50,60
// Y0                      64			60,70

  strings.push_back("\n");
}


// from definitions_multimage.h
//     ImageInfo( )
//     {
//       dataImage_present = false;
//       dataImageFileName = "";
//       weight = 1.0;
//       gain = expTime = 1.0;
//       readNoise = 0.0;
//       originalSky = 0.0;
//       nCombined = 1;
//       maskImage_present = false;
//       maskImageFileName = "";
//       errorImage_present = false;
//       errorImageFileName = "";
//       psfImage_present = false;
//       psfImageFileName = "";
// 
//       nOversampledFileNames = 0;
//       psfOversampledImage_present = false;
//       nOversamplingScales = 0;
//       //psfOversampledFileName = "";
//       //psfOversamplingScale = 0;
//       oversampleRegionSet = false;
//       nOversampleRegions = 0;
//       
//       nColumns = nRows = 0;
//       x0 = y0 = 0.0;
//       
//       pixelScale = intensityScale = 1.0;
//       rotation = 0.0;

//     vector<string>  psfOversampledFileNames;
//     int  nOversampledFileNames;
//     bool  psfOversampledImage_present;
//     vector<int>  psfOversamplingScales;
//     int  nOversamplingScales;
//     bool  oversampleRegionSet;
//     int  nOversampleRegions;
//     vector<string>  psfOversampleRegions;



/* END OF FILE: print_results_multi.cpp ---------------------------- */
