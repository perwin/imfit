// Copyright 2017-2020 by Peter Erwin.
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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>  // for strtol
#include <strings.h>  // for memset
#include <tuple>

using namespace std;

#include "utilities_pub.h"
#include "read_simple_params.h"
#include "imageparams_file_parser.h"
#include "config_file_parser.h"   // for ValidParameterLine()
#include "definitions_multimage.h"


// DATA	n100_hst-acs-wfc_ss.fits[250:350,200:300]	1.0
// GAIN		1.0
// READNOISE	0.6
// ORIGINAL_SKY	2.359
// MASK	n100_hst_mask.fits[250:350,200:300]
// PSF		n100_hst_tinytim_psf.fits
// PIXEL_SCALE		1.0		fixed   # pixels in this image are PIXEL_SCALE * reference image
// IMAGE_PA		0.0		fixed
// FLUX_SCALE		1.0		fixed
// X0				40
// X0				40

const string PARAMETER_NAME_NCOLS = "NCOLS";
const string PARAMETER_NAME_NROWS = "NROWS";
const string PARAMETER_NAME_NCOLS_NROWS = "NCOLS,NROWS";
const string PARAMETER_NAME_NROWS_NCOLS = "NROWS,NCOLS";
const string PARAMETER_NAME_WEIGHT = "WEIGHT";
const string PARAMETER_NAME_GAIN = "GAIN";
const string PARAMETER_NAME_READNOISE = "READNOISE";
const string PARAMETER_NAME_NCOMBINED = "NCOMBINED";
const string PARAMETER_NAME_EXPTIME = "EXPTIME";
const string PARAMETER_NAME_ORIGINAL_SKY = "ORIGINAL_SKY";
const string PARAMETER_NAME_DATA = "DATA";
const string PARAMETER_NAME_MASK = "MASK";
const string PARAMETER_NAME_ERROR = "ERROR";
const string PARAMETER_NAME_PSF = "PSF";

const string PARAMETER_FLAG_MASK_ZERO_IS_GOOD = "MASK_ZERO_IS_GOOD";  // default
const string PARAMETER_FLAG_MASK_ZERO_IS_BAD = "MASK_ZERO_IS_BAD";
const string PARAMETER_FLAG_ERRORS_ARE_SIGMAS = "ERRORS_ARE_SIGMAS";
const string PARAMETER_FLAG_ERRORS_ARE_VARIANCES = "ERRORS_ARE_VARIANCES";
const string PARAMETER_FLAG_ERRORS_ARE_WEIGHTS = "ERRORS_ARE_WEIGHTS";

const string PARAMETER_NAME_OVERPSF = "OVERSAMPLED_PSF";
const string PARAMETER_NAME_OVERPSF_SCALE = "OVERSAMPLE_SCALE";
const string PARAMETER_NAME_OVERPSF_REGION = "OVERSAMPLED_REGION";

const string PARAMETER_NAME_PIXEL_SCALE = "PIXEL_SCALE";
const string PARAMETER_NAME_IMAGE_PA = "IMAGE_PA";
const string PARAMETER_NAME_FLUX_SCALE = "FLUX_SCALE";
const string PARAMETER_NAME_X0 = "X0";
const string PARAMETER_NAME_Y0 = "Y0";

const string fixedIndicatorString = "fixed";



void ReportConfigError( const vector<int> &errorCodes, 
						const vector<int> &origLineNumbers );

/* ---------------- FUNCTION: ReportConfigError ------------------------ */
void ReportConfigError( const vector<int> &errorCodes, 
						const vector<int> &origLineNumbers )
{
  fprintf(stderr, "One or more errors in the input image-info file!\n");
  if (errorCodes[0] == IMAGEPARAMS_FILE_ERROR_NOLINES) {
    fprintf(stderr, "   No non-empty/non-comment lines in image-info file!\n");
    return;
  }
  if (errorCodes[0] == IMAGEPARAMS_FILE_ERROR_NOLINENUMBERS) {
    fprintf(stderr, "   [internal error] Empty line-numbers vector in call to VetInputLines!\n");
    return;
  }
  int  nErrors = (int)errorCodes.size();
  for (int i = 0; i < nErrors; i++) {
    if (errorCodes[i] == IMAGEPARAMS_FILE_ERROR_NOIMAGESTART)
      fprintf(stderr, "   No IMAGE_START lines found!\n");
    if (errorCodes[i] == IMAGEPARAMS_FILE_ERROR_BADLINE)
      fprintf(stderr, "   Invalid line at line #%d\n", origLineNumbers[i]);
  }
}



bool VetInputLines( const vector<string> &inputLines, const vector<int> &origLineNumbers,
					vector<int> &errorCodes, vector<int> &badLineNumbers )
{
  int  i, nInputLines;
  int  nErrors = 0;
  bool  imageSectionFound = false;
  bool  dataFound = false;

  // check for erroneous input (no input lines, no original line numbers)
  if (inputLines.empty()) {
    errorCodes.push_back(IMAGEPARAMS_FILE_ERROR_NOLINES);
    badLineNumbers.push_back(0);
    return false;
  }
  if (origLineNumbers.empty()) {
    errorCodes.push_back(IMAGEPARAMS_FILE_ERROR_NOLINENUMBERS);
    badLineNumbers.push_back(0);
    return false;
  }
  
  // check for required stuff
  i = 0;
  nInputLines = inputLines.size();
  while (i < nInputLines) {
    if (inputLines[i].find("IMAGE_START", 0) != string::npos) {
      imageSectionFound = true;
      break;
    }
    i++;
  }
  
  // check to see that individual lines are well-formatted
  for (i = 0; i < nInputLines; i++) {
    // need to copy string so ValidParameterLine() can accept & operate on it.
    string currentLine = inputLines[i];
    if (currentLine.find("IMAGE_START", 0) != string::npos) {
      imageSectionFound = true;
      continue;
    }
    // Currently, we have two possibilities: FUNCTIONS_START or a parameter line
    if (! (currentLine.find("FUNCTIONS_START", 0) != string::npos)) {
      if (! ValidParameterLine(currentLine, true)) {
        badLineNumbers.push_back(origLineNumbers[i]);
        errorCodes.push_back(IMAGEPARAMS_FILE_ERROR_BADLINE);
        nErrors += 1;
      }
    }
  }

  if (! imageSectionFound) {
    errorCodes.push_back(IMAGEPARAMS_FILE_ERROR_NOIMAGESTART);
    badLineNumbers.push_back(0);
    nErrors += 1;
  }

  if (nErrors > 0)
    return false;
  else
    return true;
}


/* ---------------- FUNCTION: GetImageDescriptionParams ---------------- */
// Returns info about "image-description parameters" = image-description info 
// for each image (pixelScale, rotation, intensityScale, X0, Y0).
// This is the simple version for use with makemultimages
// Returns tuple of (array of image-description parameter values for all images, status)
std::tuple<double *, int> GetImageDescriptionParams( const vector<ImageInfo> imageInfoList )
{
  double *imageDescParams;
  int  nImages = imageInfoList.size();

  int  nImageDescParams = N_IMAGE_PARAMS * (nImages - 1);
  imageDescParams = (double *)calloc(nImageDescParams, sizeof(double));

  for (int i = 0; i < (nImages - 1); i++) {
    imageDescParams[i*N_IMAGE_PARAMS] = imageInfoList[i + 1].pixelScale;
    imageDescParams[i*N_IMAGE_PARAMS + 1] = imageInfoList[i + 1].rotation;
    imageDescParams[i*N_IMAGE_PARAMS + 2] = imageInfoList[i + 1].intensityScale;
    imageDescParams[i*N_IMAGE_PARAMS + 3] = imageInfoList[i + 1].x0;
    imageDescParams[i*N_IMAGE_PARAMS + 4] = imageInfoList[i + 1].y0;
  }
  
  return std::make_tuple(imageDescParams, nImageDescParams);
}



/* ---------------- FUNCTION: FindIfParamLimitsExist ------------------- */
// Determine if parameter limits exist: either fixed = true, or limited[0]
// or limited[1] = true
bool FindIfParamLimitsExist( mp_par& parameterLimits )
{
  if ( (parameterLimits.fixed == 1) || (parameterLimits.limited[0] == 1) ||
  		(parameterLimits.limited[1] == 1) )
    return true;
  else
    return false;
}


/* ---------------- FUNCTION: GetImageDescriptionParams ---------------- */
// Returns info about "image-description parameters" = image-description info 
// for each image (pixelScale, rotation, intensityScale, X0, Y0).
// This is the more complicated version for use with multimfit
// Here, we simultaneously populate the parameter-value array (returned in the
// tuple) *and* a vector of mp_par values with corresponding fixed/limits info
// for each image-description parameter.
// Returns tuple of (array of image-description parameter values for all images, status);
// also returns parameter-limits info (vector of mp_par) in parameter imageParamLimits.
std::tuple<double *, int> GetImageDescriptionParams( const vector<ImageInfo> imageInfoList,
											vector<mp_par>& imageParamLimits, 
											bool& paramLimitsExist )
{
  double *imageDescParams;
  int  nImages = imageInfoList.size();
  mp_par  paramLimit;

//  bzero(&paramLimit, sizeof(mp_par));
  memset(&paramLimit, 0, sizeof(mp_par));
  paramLimitsExist = false;
  int  nImageDescParams = N_IMAGE_PARAMS * (nImages - 1);
  imageDescParams = (double *)calloc(nImageDescParams, sizeof(double));

  for (int i = 0; i < (nImages - 1); i++) {
    // pixelScale value and limits info
    imageDescParams[i*N_IMAGE_PARAMS] = imageInfoList[i + 1].pixelScale;
    paramLimit.fixed = imageInfoList[i + 1].pixelScale_limitInfo.fixed;
    paramLimit.limited[0] = imageInfoList[i + 1].pixelScale_limitInfo.limited[0];
    paramLimit.limited[1] = imageInfoList[i + 1].pixelScale_limitInfo.limited[1];
    paramLimit.limits[0] = imageInfoList[i + 1].pixelScale_limitInfo.limits[0];
    paramLimit.limits[1] = imageInfoList[i + 1].pixelScale_limitInfo.limits[1];
    if (FindIfParamLimitsExist(paramLimit))
      paramLimitsExist = true;
    imageParamLimits.push_back(paramLimit);
    
    // rotation value and limits info
    imageDescParams[i*N_IMAGE_PARAMS + 1] = imageInfoList[i + 1].rotation;
    paramLimit.fixed = imageInfoList[i + 1].rotation_limitInfo.fixed;
    paramLimit.limited[0] = imageInfoList[i + 1].rotation_limitInfo.limited[0];
    paramLimit.limited[1] = imageInfoList[i + 1].rotation_limitInfo.limited[1];
    paramLimit.limits[0] = imageInfoList[i + 1].rotation_limitInfo.limits[0];
    paramLimit.limits[1] = imageInfoList[i + 1].rotation_limitInfo.limits[1];
    if (FindIfParamLimitsExist(paramLimit))
      paramLimitsExist = true;
    imageParamLimits.push_back(paramLimit);

    // intensityScale value and limits info
    imageDescParams[i*N_IMAGE_PARAMS + 2] = imageInfoList[i + 1].intensityScale;
    paramLimit.fixed = imageInfoList[i + 1].intensityScale_limitInfo.fixed;
    paramLimit.limited[0] = imageInfoList[i + 1].intensityScale_limitInfo.limited[0];
    paramLimit.limited[1] = imageInfoList[i + 1].intensityScale_limitInfo.limited[1];
    paramLimit.limits[0] = imageInfoList[i + 1].intensityScale_limitInfo.limits[0];
    paramLimit.limits[1] = imageInfoList[i + 1].intensityScale_limitInfo.limits[1];
    if (FindIfParamLimitsExist(paramLimit))
      paramLimitsExist = true;
    imageParamLimits.push_back(paramLimit);

    // X0 value and limits info
    imageDescParams[i*N_IMAGE_PARAMS + 3] = imageInfoList[i + 1].x0;
    paramLimit.fixed = imageInfoList[i + 1].x0_limitInfo.fixed;
    paramLimit.limited[0] = imageInfoList[i + 1].x0_limitInfo.limited[0];
    paramLimit.limited[1] = imageInfoList[i + 1].x0_limitInfo.limited[1];
    paramLimit.limits[0] = imageInfoList[i + 1].x0_limitInfo.limits[0];
    paramLimit.limits[1] = imageInfoList[i + 1].x0_limitInfo.limits[1];
    imageParamLimits.push_back(paramLimit);
    if (FindIfParamLimitsExist(paramLimit))
      paramLimitsExist = true;

    // Y0 value and limits info
    imageDescParams[i*N_IMAGE_PARAMS + 4] = imageInfoList[i + 1].y0;
    paramLimit.fixed = imageInfoList[i + 1].y0_limitInfo.fixed;
    paramLimit.limited[0] = imageInfoList[i + 1].y0_limitInfo.limited[0];
    paramLimit.limited[1] = imageInfoList[i + 1].y0_limitInfo.limited[1];
    paramLimit.limits[0] = imageInfoList[i + 1].y0_limitInfo.limits[0];
    paramLimit.limits[1] = imageInfoList[i + 1].y0_limitInfo.limits[1];
    if (FindIfParamLimitsExist(paramLimit))
      paramLimitsExist = true;
    imageParamLimits.push_back(paramLimit);
  }

  return std::make_tuple(imageDescParams, nImageDescParams);
}



// status = 0 for just a parameter value, 1 for fixed parameter,
// 2 for parameter value and limit(s), or -1 for problem
std::tuple<double, mp_par, int> GetParameterAndLimits( const string theLine )
{
  double  paramVal;
  double  lowerLimit, upperLimit;
  string  extraPiece;
  vector<string>  stringPieces, newPieces;
  mp_par  newParamLimit;
  int  status = 0;

  stringPieces.clear();
  SplitString(theLine, stringPieces);
  // first piece is parameter name, which we ignore; second piece is initial value
  paramVal = strtod(stringPieces[1].c_str(), NULL);

  // OK, now we create a new mp_par structure and check for possible parameter limits
//   bzero(&newParamLimit, sizeof(mp_par));
  memset(&newParamLimit, 0, sizeof(mp_par));
  if (stringPieces.size() > 2) {
    // parse and store parameter limits, if any
    extraPiece = stringPieces[2];
    if (extraPiece == fixedIndicatorString) {
      newParamLimit.fixed = 1;
      status = 1;
    } else {
      if (extraPiece.find(',', 0) != string::npos) {
        status = 2;
        newPieces.clear();
        SplitString(extraPiece, newPieces, ",");
        newParamLimit.limited[0] = 1;
        newParamLimit.limited[1] = 1;
        lowerLimit = strtod(newPieces[0].c_str(), NULL);
        upperLimit = strtod(newPieces[1].c_str(), NULL);
        if (lowerLimit >= upperLimit) {
          fprintf(stderr, "*** WARNING: first parameter limit for \"%s\" (%g) must be < second limit (%g)!\n",
          				stringPieces[0].c_str(), lowerLimit, upperLimit);
          if (lowerLimit == upperLimit) {
             fprintf(stderr, "    (To specify fixed parameters, use \"fixed\" keyword, not identical upper & lower parameter limits)\n");
          }
          status = -1;
        }
        if ((paramVal < lowerLimit) || (paramVal > upperLimit)) {
          fprintf(stderr, "*** WARNING: initial parameter value for \"%s\" (%g) must lie between the limits (%g,%g)!\n",
          				stringPieces[0].c_str(), paramVal, lowerLimit, upperLimit);
          status = -1;
        }
        newParamLimit.limits[0] = lowerLimit;
        newParamLimit.limits[1] = upperLimit;
      }
    }
  }

  return std::make_tuple(paramVal, newParamLimit, status);
}



/* ---------------- FUNCTION: StoreInfoFromLine ------------------------ */
// Given a line from an image-info file, this function parses the line and stores
// appropriate info in the input ImageInfo object
// We assume that theLine has had leading and trailing whitespace (including
// newlines) removed before this function is called.
int StoreInfoFromLine( const string theLine, ImageInfo& imageInfo )
{
  double  numericVal;
  int  status;
  vector<string>  stringPieces;
  mp_par  paramLimits;
  string  parameterName, stringVal;
  
  SplitString(theLine, stringPieces);
  parameterName = stringPieces[0];
    
  // floating-point values
  if (parameterName == PARAMETER_NAME_GAIN) 
    imageInfo.gain = strtod(stringPieces[1].c_str(), NULL);
      
  else if (parameterName == PARAMETER_NAME_READNOISE) 
    imageInfo.readNoise = strtod(stringPieces[1].c_str(), NULL);
      
  else if (parameterName == PARAMETER_NAME_ORIGINAL_SKY) 
    imageInfo.originalSky = strtod(stringPieces[1].c_str(), NULL);
      
  else if (parameterName == PARAMETER_NAME_EXPTIME) 
    imageInfo.expTime = strtod(stringPieces[1].c_str(), NULL);
      
  else if (parameterName == PARAMETER_NAME_WEIGHT) 
    imageInfo.weight = strtod(stringPieces[1].c_str(), NULL);


  // Special cases of parameters with optional limits or fixed status
  
  else if (parameterName == PARAMETER_NAME_PIXEL_SCALE) {
    std::tie(numericVal, paramLimits, status) = GetParameterAndLimits(theLine);
    imageInfo.pixelScale = numericVal;
    // special rule for PIXEL_SCALE: fixed unless limits were explicitly set
    if (status == 2)
      imageInfo.pixelScale_limitInfo = paramLimits;
    else
      imageInfo.pixelScale_limitInfo.fixed = 1;
  }
      
  else if (parameterName == PARAMETER_NAME_IMAGE_PA) {
    std::tie(numericVal, paramLimits, status) = GetParameterAndLimits(theLine);
    imageInfo.rotation = numericVal;
    // special rule for IMAGE_PA: fixed unless limits were explicitly set
    if (status == 2)
      imageInfo.rotation_limitInfo = paramLimits;
    else
      imageInfo.rotation_limitInfo.fixed = 1;
  }
      
  else if (parameterName == PARAMETER_NAME_FLUX_SCALE) {
    std::tie(numericVal, paramLimits, status) = GetParameterAndLimits(theLine);
    imageInfo.intensityScale = numericVal;
    imageInfo.intensityScale_limitInfo = paramLimits;
  }
      
  else if (parameterName == PARAMETER_NAME_X0) {
    std::tie(numericVal, paramLimits, status) = GetParameterAndLimits(theLine);
    imageInfo.x0 = numericVal;
    imageInfo.x0_limitInfo = paramLimits;
  }
      
  else if (parameterName == PARAMETER_NAME_Y0) {
    std::tie(numericVal, paramLimits, status) = GetParameterAndLimits(theLine);
    imageInfo.y0 = numericVal;
    imageInfo.y0_limitInfo = paramLimits;
  }
  
  
  // integer values
  else if (parameterName == PARAMETER_NAME_NCOLS)
    imageInfo.nColumns = (int)strtol(stringPieces[1].c_str(), NULL, 0);

  else if (parameterName == PARAMETER_NAME_NROWS) 
    imageInfo.nRows = (int)strtol(stringPieces[1].c_str(), NULL, 0);

  else if (parameterName == PARAMETER_NAME_NROWS_NCOLS) {
    vector<string> subPieces;
    SplitString(stringPieces[1], subPieces, ",");
    imageInfo.nRows = (int)strtol(subPieces[0].c_str(), NULL, 0);
    imageInfo.nColumns = (int)strtol(subPieces[1].c_str(), NULL, 0);
  }

  else if (parameterName == PARAMETER_NAME_NCOLS_NROWS) {
    vector<string> subPieces;
    SplitString(stringPieces[1], subPieces, ",");
    imageInfo.nColumns = (int)strtol(subPieces[0].c_str(), NULL, 0);
    imageInfo.nRows = (int)strtol(subPieces[1].c_str(), NULL, 0);
  }

  else if (parameterName == PARAMETER_NAME_NCOMBINED)
    imageInfo.nCombined = (int)strtol(stringPieces[1].c_str(), NULL, 0);

  
  // string values
  else if (parameterName == PARAMETER_NAME_DATA) {
    imageInfo.dataImageFileName = stringPieces[1];
    imageInfo.dataImage_present = true;
  }

  else if (parameterName == PARAMETER_NAME_MASK) {
    imageInfo.maskImageFileName = stringPieces[1];
    imageInfo.maskImage_present = true;
  }

  else if (parameterName == PARAMETER_NAME_ERROR) {
    imageInfo.errorImageFileName = stringPieces[1];
    imageInfo.errorImage_present = true;
  }

  else if (parameterName == PARAMETER_NAME_PSF) {
    imageInfo.psfImageFileName = stringPieces[1];
    imageInfo.psfImage_present = true;
      }

  // Special cases of oversampled PSF specifications
  else if (parameterName == PARAMETER_NAME_OVERPSF_SCALE) {
    int osampleScale = (int)strtol(stringPieces[1].c_str(), NULL, 0);
    imageInfo.psfOversamplingScales.push_back(osampleScale);
    imageInfo.nOversamplingScales += 1;
  }

  else if (parameterName == PARAMETER_NAME_OVERPSF) {
    imageInfo.psfOversampledFileNames.push_back(stringPieces[1]);
    imageInfo.nOversampledFileNames += 1;
    imageInfo.psfOversampledImage_present = true;
  }
  
  else if (parameterName == PARAMETER_NAME_OVERPSF_REGION) {
    imageInfo.psfOversampleRegions.push_back(stringPieces[1]);
    imageInfo.nOversampleRegions += 1;
    imageInfo.oversampleRegionSet = true;
  }
  
  // Miscellaneous flags
  else if (parameterName == PARAMETER_FLAG_MASK_ZERO_IS_GOOD)   // default
    imageInfo.maskFormat = MASK_ZERO_IS_GOOD;
  else if (parameterName == PARAMETER_FLAG_MASK_ZERO_IS_BAD)
    imageInfo.maskFormat = MASK_ZERO_IS_BAD;
  else if (parameterName == PARAMETER_FLAG_ERRORS_ARE_SIGMAS)   // default
    imageInfo.errorType = WEIGHTS_ARE_SIGMAS;
  else if (parameterName == PARAMETER_FLAG_ERRORS_ARE_VARIANCES)
    imageInfo.errorType = WEIGHTS_ARE_VARIANCES;
  else if (parameterName == PARAMETER_FLAG_ERRORS_ARE_WEIGHTS)
    imageInfo.errorType = WEIGHTS_ARE_WEIGHTS;


  // FIXME: MISSING: 
  // optional codes for mask and error files (e.g., equivalents for 
  // "--errors-are-variances" and "--mask-zero-is-bad";
  // [x]exposure time; [x]nCombined
  else if (parameterName == PARAMETER_NAME_DATA) 
    imageInfo.dataImageFileName = stringPieces[1];

  else {
    ; // handle error
  }
  
  return 0;
}


/* ---------------- FUNCTION: AddInfoFromLines ------------------------- */
// Given a list of (valid) lines (comments and blank lines removed, etc.) from 
// an image-info file, this function parses each line (using StoreInfoFromLine) 
// and populates the input vector of ImageInfo objects, one for each specified
// image.
// Returns 0 for standard operation (no per-image functions found); returns 
// index (> 0) for first line of per-image-function section if it exists
int AddInfoFromLines( const vector<string>& lines, vector<ImageInfo>& imageInfoList )
{
  ImageInfo newImageInfo = ImageInfo();
  int  nLines = lines.size();
  int  returnVal = 0;
  
  for (int i = 0; i < nLines; i++) {
    if (lines[i].find("FUNCTIONS_START", 0) != string::npos) {
      returnVal = i;
      break;
    }
    StoreInfoFromLine(lines[i], newImageInfo);
  }
  imageInfoList.push_back(newImageInfo);
  
  return returnVal;
}


/* ---------------- FUNCTION: ReadImageParamsFile ---------------------- */
// Reads and parses an image-info configuration file.
// Returns (in input parameters):
//    Number of images specified in the file (nImages)
//    Vector of ImageInfo objects (one per specified image) 
int ReadImageParamsFile( const string& configFileName, int& nImages,
						vector<ImageInfo>& imageInfoList )
{
  ifstream  inputFileStream;
  string  inputLine;
  vector<string>  inputLines;
  vector<string>  stringPieces;
  vector<int>  origLineNumbers;
  vector<int>  imageStartIndices;
  vector<string>  funcNameList, funcLabelList;
  vector<double>  funcParams;
  vector<mp_par>  funcParamLimits;
  vector<int>  funcSetStartIndices;
  int  imageFunctionsOffset, imageFunctionsStartIndex;
  int  i, k, nInputLines;

  // Read in and store non-empty lines (also store 1-based line numbers from
  // original file, accounting for empty lines)
  inputFileStream.open(configFileName.c_str());
  if( ! inputFileStream ) {
     cerr << "Error opening input stream for file " << configFileName.c_str() << endl;
  }
  int  currentLine = 0;   // 1-based for user
  while ( getline(inputFileStream, inputLine) ) {
    currentLine++;
    // strip off leading & trailing spaces; turns a blank line with spaces/tabs
    // into an empty string
    TrimWhitespace(inputLine);
    // store non-empty, non-comment lines in a vector of strings
    if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
      inputLines.push_back(inputLine);
      origLineNumbers.push_back(k);
    }
  }
  inputFileStream.close();
  nInputLines = inputLines.size();

  vector<int> badLineNumbers;
  vector<int> errorCodes;
  bool inputOK = VetInputLines(inputLines, origLineNumbers, errorCodes, badLineNumbers);
  if (! inputOK) {
    ReportConfigError(errorCodes, badLineNumbers);
    return -1;
  }
  
  // Clear the input vectors before we start appending things to them
  imageInfoList.clear();

  // look for IMAGE_START to chop up input into per-image sections
  for (i = 0; i < nInputLines; i++) {
    if (inputLines[i].find("IMAGE_START", 0) != string::npos)
      imageStartIndices.push_back(i);
  }
  nImages = imageStartIndices.size();
  
  // Step through each block of per-image lines and analyze them
  vector<string> currentImageLines;
  int  currentStartIndex, currentEndIndex;
  vector<string> imageFunctionsLines;
  for (k = 0; k < nImages; k++) {
    currentStartIndex = imageStartIndices[k];
    if (k < (nImages - 1))
      currentEndIndex = imageStartIndices[k + 1];
    else
      currentEndIndex = nInputLines;
    currentImageLines.clear();
    for (i = currentStartIndex; i < currentEndIndex; i++)
      currentImageLines.push_back(inputLines[i]);
    // start of per-image functions section is relative to currentImageLines (0 = no such section)
    imageFunctionsOffset = AddInfoFromLines(currentImageLines, imageInfoList);
    // process per-image function lines, if they exist
    if (imageFunctionsOffset > 0) {
      // start with global location of current per-image functions section
      imageFunctionsStartIndex = currentStartIndex + imageFunctionsOffset;
      vector<int>  functionSectionLineNumbers;
      imageFunctionsLines.clear();
      for (i = imageFunctionsStartIndex + 1; i < currentEndIndex; i++) {
        imageFunctionsLines.push_back(inputLines[i]);
        functionSectionLineNumbers.push_back(i);
      }
      bool  paramLimitsExist = false;
      int status = ParseFunctionSection(imageFunctionsLines, true, funcNameList,
						funcLabelList, funcParams, funcParamLimits, funcSetStartIndices,
						paramLimitsExist, functionSectionLineNumbers);
	  if (status < 0) {
	    fprintf(stderr, "ERROR MESSAGE ABOUT FAILURE IN ParseFunctionSection\n");
	    return -1;
	  }
	  imageInfoList[k].perImageFuncNames = funcNameList;
	  imageInfoList[k].perImageFuncLabels = funcLabelList;
	  imageInfoList[k].perImageFuncSetIndices = funcSetStartIndices;
	  imageInfoList[k].perImageParamVals = funcParams;
	  imageInfoList[k].perImageParamLimitsExist = paramLimitsExist;
	  imageInfoList[k].perImageParamLimits = funcParamLimits;
	  imageInfoList[k].perImageFunctionsExist = true;
    }
 }
  
  return 0;
}
