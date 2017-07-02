/** @file
    \brief Subclass of OptionsBase blass for makeimage-specific program options.
 *
 */
#ifndef _OPTIONS_MAKEIMAGE_H_
#define _OPTIONS_MAKEIMAGE_H_

#include <string>
#include <vector>

using namespace std;

#include "definitions.h"
#include "options_base.h"


//! Derived class for holding various options for makeimage (set by command-line flags & options)
class MakeimageOptions : public OptionsBase
{

  public:
    // Constructor:
    MakeimageOptions( )
    {
      configFileName = "";

      noOutputImageName = true;
      outputImageName = DEFAULT_MAKEIMAGE_OUTPUT_FILENAME;
      functionRootName = "";
      noRefImage = true;
      referenceImageName = "";
  
      noImageDimensions = true;
      nColumnsSet = false;
      nRowsSet = false;
      nColumns = 0;
      nRows = 0;
  
      magZeroPoint = NO_MAGNITUDES;

      printImages = false;
      saveImage = true;
      saveExpandedImage = false;
      saveAllFunctions = false;

      printFluxes = false;
      estimationImageSize = DEFAULT_ESTIMATION_IMAGE_SIZE;
      timingIterations = 0;
    };

    // Extra data members (in addition to those in options_base.h):
    bool  noOutputImageName;
    string  outputImageName;

    string  functionRootName;
    bool  noRefImage;
    string  referenceImageName;

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
