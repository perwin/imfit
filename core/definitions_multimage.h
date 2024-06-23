/** @file
    \brief Generally useful definitions & constants for multimfit, etc.

 */

#ifndef _DEFINITIONS_MULTIMAGE_H_
#define _DEFINITIONS_MULTIMAGE_H_

#include <string>
#include <vector>
#include <strings.h>

#include "definitions.h"
#include "param_struct.h"

using namespace std;

const int N_IMAGE_DESCRIPTION_PARAMS = 3;  // pixScale, rotation, intensityScale
const int N_IMAGE_PARAMS = 5;  // pixScale, rotation, intensityScale, X0_0, Y0_0

#define DEFAULT_MAKEMULTIMAGE_OUTPUT_FILENAME   "modelimage_multi"

#define DEFAULT_BESTFIT_ROOTNAME   "bestfit_parameters_multimfit"




class ImageInfo
{

  public:
    // Constructor:
    ImageInfo( )
    {
      dataImage_present = false;
      dataImageFileName = "";
      weight = 1.0;
      gain = expTime = 1.0;
      readNoise = 0.0;
      originalSky = 0.0;
      nCombined = 1;
      maskImage_present = false;
      maskImageFileName = "";
      maskFormat = MASK_ZERO_IS_GOOD;
      errorImage_present = false;
      errorImageFileName = "";
      errorType = WEIGHTS_ARE_SIGMAS;
      psfImage_present = false;
      psfImageFileName = "";

      nOversampledFileNames = 0;
      psfOversampledImage_present = false;
      nOversamplingScales = 0;
      //psfOversampledFileName = "";
      //psfOversamplingScale = 0;
      oversampleRegionSet = false;
      nOversampleRegions = 0;
      
      nColumns = nRows = 0;
      x0 = y0 = 0.0;
      
      pixelScale = intensityScale = 1.0;
      rotation = 0.0;
      
      bzero(&pixelScale_limitInfo, sizeof(mp_par));
      bzero(&rotation_limitInfo, sizeof(mp_par));
      bzero(&intensityScale_limitInfo, sizeof(mp_par));
      bzero(&x0_limitInfo, sizeof(mp_par));
      bzero(&y0_limitInfo, sizeof(mp_par));
      
      perImageFunctionsExist = false;
    }
    
    string  dataImageFileName;
    string  maskImageFileName;
    int  maskFormat;
    string  errorImageFileName;
    int  errorType;
    string  psfImageFileName;
    double  weight, gain, readNoise, originalSky, expTime;
    int  nCombined;
    double  x0, y0;
    int  nColumns, nRows;
    double  pixelScale, intensityScale, rotation;

    bool  dataImage_present, maskImage_present, errorImage_present;
    bool  psfImage_present;

    vector<string>  psfOversampledFileNames;
    int  nOversampledFileNames;
    bool  psfOversampledImage_present;
    vector<int>  psfOversamplingScales;
    int  nOversamplingScales;
    bool  oversampleRegionSet;
    int  nOversampleRegions;
    vector<string>  psfOversampleRegions;
    
    // stuff for per-image functions, if any
    bool perImageFunctionsExist;
    vector<string> perImageFuncNames;
    vector<string> perImageFuncLabels;
    vector<int> perImageFuncSetIndices;
    vector<double> perImageParamVals;
    bool perImageParamLimitsExist;
    vector<mp_par> perImageParamLimits;
    

    mp_par  pixelScale_limitInfo;
    mp_par  rotation_limitInfo;
    mp_par  intensityScale_limitInfo;
    mp_par  x0_limitInfo;
    mp_par  y0_limitInfo;
};


#endif /* _DEFINITIONS_MULTIMAGE_H_ */
