/** @file
    \brief Base class for program options; meant to be subclassed by imfit,
    imfit-mcmc, and makeimage
 *
 *  Fundamentally, this is a list of parameter values and settings. For ease
 *  of access, all data members are public. (The only reason we don't use a C-style
 *  struct is so we can use the constructor to set default values automatically,
 *  and so we can have derived classes with the same base set of parameters.)
 *
 *  In addition to serving as a base class for expanded, program-specific option
 *  listing, this interface is meant for use by SetupModelImage(), so that any
 *  instance of a derived class can be used by that function.
 */
#ifndef _OPTIONS_BASE_H_
#define _OPTIONS_BASE_H_

#include <string>
#include <vector>

using namespace std;

#include "definitions.h"


//! base class for holding various imfit options (set by command-line flags & options)
class OptionsBase
{

  public:
    // Constructor:
    OptionsBase( )
    {
      noConfigFile = true;
      
      psfImagePresent = false;
      psfFileName = "";
      psfOversampledImagePresent = false;
      psfOversampledFileName = "";
      psfOversamplingScale = 0;
      oversampleRegionSet = false;
      nOversampleRegions = 0;
      psfOversampleRegion = "";
      
      noiseImagePresent = false;
      noiseFileName = "";
      errorType = WEIGHTS_ARE_SIGMAS;
  
      maskImagePresent = false;
      maskFileName = "";
      maskFormat = MASK_ZERO_IS_GOOD;

      gainSet = false;
      gain = 1.0;
      readNoiseSet = false;
      readNoise = 0.0;
      expTimeSet = false;
      expTime = 1.0;
      nCombinedSet = false;
      nCombined = 1;
      originalSkySet = false;
      originalSky = 0.0;

      printFitStatisticOnly = false;
      solver = MPFIT_SOLVER;

      subsamplingFlag = true;

      rngSeed = 0;           // 0 = get seed value from system clock
  
      maxThreads = 0;
      maxThreadsSet = false;

      verbose = 1;
      debugLevel = 0;
    };

    // Data members:
    bool  noConfigFile;
    string  configFileName;   // []

    bool  psfImagePresent;
    string  psfFileName;     // []
  
    bool  psfOversampledImagePresent;
    string  psfOversampledFileName;     // []
    vector<string>  psfOversampledFileNames;     // []
    int  psfOversamplingScale;
    vector<int>  psfOversamplingScales;
    bool  oversampleRegionSet;
    int  nOversampleRegions;
    string  psfOversampleRegion;     // []
    vector<string>  psfOversampleRegions;     // []
  
    bool  noiseImagePresent;
    string  noiseFileName;   // []
    int  errorType;
  
    string  maskFileName;   //  []
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

    bool printFitStatisticOnly;
    int  solver;
    bool  useModelForErrors;
    bool  useCashStatistic;
    bool  usePoissonMLR;

    int  maxThreads;
    bool  maxThreadsSet;
  
    unsigned long  rngSeed;

    int  debugLevel;
    int  verbose;
};


#endif  // _OPTIONS_BASE_H_
