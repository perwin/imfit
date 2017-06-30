/** @file
    \brief Subclass of OptionsBase class for imfit-specific program options.
 *
 */
#ifndef _OPTIONS_MCMC_H_
#define _OPTIONS_MCMC_H_

#include <string>
#include <vector>

using namespace std;

#include "definitions.h"
#include "options_base.h"


//! Derived class for holding various options for imfit (set by command-line flags & options)
class MCMCOptions : public OptionsBase
{

  public:
    // Constructor:
    MCMCOptions( )
    {
      configFileName = DEFAULT_IMFIT_CONFIG_FILE;
  
      noImage = true;
      imageFileName = "";
  
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
  
      subsamplingFlag = true;
  
      saveModel = false;
      saveResidualImage = false;

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

      useModelForErrors = false;
      useCashStatistic = false;
      usePoissonMLR = false;

      noParamLimits = true;

      appendToOutput = false;
      outputFileRoot = "mcmc_out";
      nChains = -1;          // -1 = use default, which is nChains = nFreeParams
      maxEvals = 100000;
      nBurnIn = 5000;
      nGelmanEvals = 5000;
      GRScaleReductionLimit = 1.01;
      mcmcNoise = 0.01;      // b parameter in DREAM (uniform scaling w/in 1 +/- b)
                        	 // 0.01 seems to work well for (small) image fits
      mcmc_bstar = 1.0e-6;   // b^star parameter in DREAM (sigma for epsilon)
                             // 1.0e-6 to 1.0e-3 seem to work ~ equally well; 0.01 is worse  
    };

    // Extra data members (in addition to those in options_base.h):  
    bool  noImage;
    string  imageFileName;
  
    string  psfOversampleRegion;
  
    // NOTE: the following are necessary as inputs for EstimateMemoryUse(), even
    // though we don't use them otherwise
    bool  saveModel;
    bool  saveResidualImage;

    bool  noParamLimits;
  
    // MCMC-related stuff
    bool  appendToOutput;
    string  outputFileRoot;
    int  nChains;
    int  maxEvals;
    int  nBurnIn;
    int  nGelmanEvals;
    double  GRScaleReductionLimit;
    double  mcmcNoise;
    double  mcmc_bstar;

};


#endif  // _OPTIONS_MCMC_H_
