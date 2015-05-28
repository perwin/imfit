/*! \file
    \brief Struct containing options info for makeimage, along with function for
           setting up default values within the struct.
 *
 */

#ifndef _MAKEIMAGE_OPTION_STRUCT_H_
#define _MAKEIMAGE_OPTION_STRUCT_H_

#include <string>

#include "definitions.h"


//! struct for holding various makeimage options (set by command-line flags & options)
typedef struct {
  bool  noConfigFile;
  std::string  configFileName;   // []
  
  bool  noImageName;
  std::string  outputImageName;
  std::string  functionRootName;   // root name for individual-image functions  []
  bool  noRefImage;
  std::string  referenceImageName;   // []

  bool  subsamplingFlag;

  bool  noImageDimensions;
  bool  nColumnsSet;
  bool  nRowsSet;
  int  nColumns;
  int  nRows;

  bool  psfImagePresent;
  std::string  psfFileName;   // []
  bool  psfOversampledImagePresent;
  std::string  psfOversampledFileName;     // []
  int  psfOversamplingScale;
  bool  oversampleRegionSet;
  std::string  psfOversampleRegion;     // []

  double  magZeroPoint;
  
  bool  printImages;
  bool  saveImage;
  bool  saveExpandedImage;
  bool  saveAllFunctions;  // save individual-function images

  bool  printFluxes;
  int  estimationImageSize;
  
  bool  maxThreadsSet;
  int  maxThreads;
  
  int  debugLevel;
} makeimageCommandOptions;




void SetDefaultMakeimageOptions( makeimageCommandOptions *theOptions )
{
  theOptions->noConfigFile = true;
  theOptions->configFileName = "";

  theOptions->noImageName = true;
  theOptions->outputImageName = DEFAULT_MAKEIMAGE_OUTPUT_FILENAME;
  theOptions->functionRootName = "";
  theOptions->noRefImage = true;
  theOptions->referenceImageName = "";
  
  theOptions->subsamplingFlag = true;
  
  theOptions->noImageDimensions = true;
  theOptions->nColumnsSet = false;
  theOptions->nRowsSet = false;
  theOptions->nColumns = 0;
  theOptions->nRows = 0;
  
  theOptions->psfImagePresent = false;
  theOptions->psfFileName = "";
  theOptions->psfOversampledImagePresent = false;
  theOptions->psfOversampledFileName = "";
  theOptions->psfOversamplingScale = 0;
  theOptions->oversampleRegionSet = false;
  theOptions->psfOversampleRegion = "";
  
  theOptions->magZeroPoint = NO_MAGNITUDES;

  theOptions->printImages = false;
  theOptions->saveImage = true;
  theOptions->saveExpandedImage = false;
  theOptions->saveAllFunctions = false;

  theOptions->printFluxes = false;
  theOptions->estimationImageSize = DEFAULT_ESTIMATION_IMAGE_SIZE;

  theOptions->maxThreadsSet = false;
  theOptions->maxThreads = 0;

  theOptions->debugLevel = 0;
}



#endif  // _MAKEIMAGE_OPTION_STRUCT_H_
