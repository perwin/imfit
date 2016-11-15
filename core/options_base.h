/** @file
    \brief Base class for program options; meant to be subclassed by imfit,
    imfit-mcmc, and makeimage
 *
 */
#ifndef _OPTIONS_BASE_H_
#define _OPTIONS_BASE_H_

#include <string>
#include <vector>

#include "definitions.h"


//! base class for holding various imfit options (set by command-line flags & options)
class OptionsBase
{

  public:
    // Constructors:
    OptionsBase( ) {;};

    // Data members:
    bool  noConfigFile;
    std::string  configFileName;   // []

    bool  psfImagePresent;
    std::string  psfFileName;     // []
  
    bool  psfOversampledImagePresent;
    std::string  psfOversampledFileName;     // []
    int  psfOversamplingScale;
    bool  oversampleRegionSet;
    int  nOversampleRegions;
    std::string  psfOversampleRegion;     // []
    std::vector<std::string>  psfOversampleRegions;     // []
  
    bool  noiseImagePresent;
    std::string  noiseFileName;   // []
    int  errorType;
  
    std::string  maskFileName;   //  []
    bool  maskImagePresent;
    int  maskFormat;
  
    bool  subsamplingFlag;

    bool  gainSet;
    double  gain;
    bool  readNoiseSet;
    double  readNoise;
    bool  expTimeSet;
    double  expTime;
    bool  nCombinedSet;
    int  nCombined;
    bool  originalSkySet;
    double  originalSky;

    bool  useModelForErrors;
    bool  useCashStatistic;
    bool  usePoissonMLR;

    int  maxThreads;
    bool  maxThreadsSet;
  
    int  debugLevel;
    int  verbose;
};


#endif  // _OPTIONS_BASE_H_
