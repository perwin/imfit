/** @file
    \brief Base class for program options; meant to be subclassed by imfit,
    imfit-mcmc, and makeimage
 *
 */
#ifndef _OPTIONS_MAKEIMAGE_H_
#define _OPTIONS_MAKEIMAGE_H_

#include <string>
#include <vector>

#include "definitions.h"
#include "options_base.h"


//! Derived class for holding various options for makeimage (set by command-line flags & options)
class MakeimageOptions : public OptionsBase
{

  public:
    // Constructors:
    MakeimageOptions( )
    {
      noConfigFile = true;
      configFileName = "";

      noOutputImageName = true;
      outputImageName = DEFAULT_MAKEIMAGE_OUTPUT_FILENAME;
      functionRootName = "";
      noRefImage = true;
      referenceImageName = "";
  
      subsamplingFlag = true;
  
      noImageDimensions = true;
      nColumnsSet = false;
      nRowsSet = false;
      nColumns = 0;
      nRows = 0;
  
      psfImagePresent = false;
      psfFileName = "";
      psfOversampledImagePresent = false;
      psfOversampledFileName = "";
      psfOversamplingScale = 0;
      oversampleRegionSet = false;
      nOversampleRegions = 0;
  
      magZeroPoint = NO_MAGNITUDES;

      printImages = false;
      saveImage = true;
      saveExpandedImage = false;
      saveAllFunctions = false;

      printFluxes = false;
      estimationImageSize = DEFAULT_ESTIMATION_IMAGE_SIZE;
      timingIterations = 0;

      maxThreadsSet = false;
      maxThreads = 0;

      debugLevel = 0;
    };

    // Extra data members (in addition to those in options_base.h):
    bool  noOutputImageName;
    std::string  outputImageName;   // [] = assign default value in main?

    std::string  functionRootName;   // root name for individual-image functions  []
    bool  noRefImage;
    std::string  referenceImageName;   // []

    bool  noImageDimensions;
    bool  nColumnsSet;
    bool  nRowsSet;
    int  nColumns;
    int  nRows;

    double  magZeroPoint;

    bool  printImages;
    bool  saveImage;
    bool  saveExpandedImage;
    bool  saveAllFunctions;  // save individual-function images

    bool  printFluxes;
    int  estimationImageSize;
  
    int  timingIterations;
};


#endif  // _OPTIONS_MAKEIMAGE_H_
