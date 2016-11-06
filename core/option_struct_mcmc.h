/** @file
    \brief Struct containing options info for imfit-mcmc, along with function for
           setting up default values within the struct.
 *
 */
#ifndef _MCMC_OPTION_STRUCT_H_
#define _MCMC_OPTION_STRUCT_H_

#include <string>
#include <vector>

#include "definitions.h"


//! struct for holding various imfit-mcmc options (set by command-line flags & options)
typedef struct {
  std::string  configFileName;
  
  bool  noImage;
  std::string  imageFileName;   // [] = assign default value in main?
  
  bool  psfImagePresent;
  std::string  psfFileName;     // []
  
  bool  psfOversampledImagePresent;
  std::string  psfOversampledFileName;     // []
  int  psfOversamplingScale;
  bool  oversampleRegionSet;
  int  nOversampleRegions;
  std::string  psfOversampleRegion;     // []
  
  bool  noiseImagePresent;
  std::string  noiseFileName;   // []
  int  errorType;
  
  std::string  maskFileName;   //  []
  bool  maskImagePresent;
  int  maskFormat;
  
  bool  subsamplingFlag;
  
  // NOTE: the following are necessary as inputs for EstimateMemoryUse(), even
  // though we don't use them otherwise
  bool  saveModel;
  bool  saveResidualImage;
  
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

  bool  noParamLimits;
  
  // MCMC-related stuff
  bool  appendToOutput;
  std::string  outputFileRoot;
  int  nChains;
  int  maxEvals;
  int  nBurnIn;
  int  nGelmanEvals;
  double  GRScaleReductionLimit;
  double  mcmcNoise;
  double  mcmc_bstar;

  unsigned long  rngSeed;
  
  int  maxThreads;
  bool  maxThreadsSet;
  
  int  verbose;
} mcmcCommandOptions;


void SetDefaultMCMCOptions( mcmcCommandOptions *theOptions )
{
  theOptions->configFileName = DEFAULT_IMFIT_CONFIG_FILE;
  
  theOptions->noImage = true;
  theOptions->imageFileName = "";
  
  theOptions->psfImagePresent = false;
  theOptions->psfFileName = "";
  
  theOptions->psfOversampledImagePresent = false;
  theOptions->psfOversampledFileName = "";
  theOptions->psfOversamplingScale = 0;
  theOptions->oversampleRegionSet = false;
  theOptions->nOversampleRegions = 0;
  theOptions->psfOversampleRegion = "";
  
  theOptions->noiseImagePresent = false;
  theOptions->noiseFileName = "";
  theOptions->errorType = WEIGHTS_ARE_SIGMAS;
  
  theOptions->maskImagePresent = false;
  theOptions->maskFileName = "";
  theOptions->maskFormat = MASK_ZERO_IS_GOOD;
  
  theOptions->subsamplingFlag = true;
  
  theOptions->saveModel = false;
  theOptions->saveResidualImage = false;

  theOptions->gainSet = false;
  theOptions->gain = 1.0;
  theOptions->readNoiseSet = false;
  theOptions->readNoise = 0.0;
  theOptions->expTimeSet = false;
  theOptions->expTime = 1.0;
  theOptions->nCombinedSet = false;
  theOptions->nCombined = 1;
  theOptions->originalSkySet = false;
  theOptions->originalSky = 0.0;

  theOptions->useModelForErrors = false;
  theOptions->useCashStatistic = false;
  theOptions->usePoissonMLR = false;

  theOptions->noParamLimits = true;

  theOptions->appendToOutput = false;
  theOptions->outputFileRoot = "mcmc_out";
  theOptions->nChains = -1;          // -1 = use default, which is nChains = nFreeParams
  theOptions->maxEvals = 100000;
  theOptions->nBurnIn = 5000;
  theOptions->nGelmanEvals = 1000;
  theOptions->GRScaleReductionLimit = 1.01;
  theOptions->mcmcNoise = 0.01;      // b parameter in DREAM (uniform scaling w/in 1 +/- b)
                                     // 0.01 seems to work well for (small) image fits
  theOptions->mcmc_bstar = 1.0e-6;   // b^star parameter in DREAM (sigma for epsilon)
                                     // 1.0e-6 to 1.0e-3 seem to work ~ equally well; 0.01 is worse
  theOptions->rngSeed = 0;           // 0 = get seed value from system clock
  
  theOptions->maxThreads = 0;
  theOptions->maxThreadsSet = false;

  theOptions->verbose = 1;
}



#endif  // _MCMC_OPTION_STRUCT_H_
