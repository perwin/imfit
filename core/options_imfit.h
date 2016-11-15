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
class ImfitOptions : public OptionsBase
{

  public:
    // Constructors:
    ImfitOptions( )
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
      outputModelFileName = "";
      saveResidualImage = false;
      outputResidualFileName = "";
      saveWeightImage = false;
      outputWeightFileName = "";
      saveBestFitParams = true;
      outputParameterFileName = DEFAULT_OUTPUT_PARAMETER_FILE;

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

      printFitStatisticOnly = false;
      noParamLimits = true;
      ftolSet = false;
      ftol = DEFAULT_FTOL;
      solver = MPFIT_SOLVER;
      nloptSolverName = "NM";   // default value = Nelder-Mead Simplex

      magZeroPoint = NO_MAGNITUDES;
  
      printImages = false;

      doBootstrap = false;
      bootstrapIterations = 0;
      saveBootstrap = false;
      outputBootstrapFileName = "";

      rngSeed = 0;
  
      maxThreads = 0;
      maxThreadsSet = false;

      verbose = 1;
    };

    // Extra data members (in addition to those in options_base.h):  
    bool  noImage;
    std::string  imageFileName;   // [] = assign default value in main?
  
    std::string  psfOversampleRegion;     // []
  
    bool  saveModel;
    std::string  outputModelFileName;   // []
    bool  saveResidualImage;
    std::string  outputResidualFileName;   // []
    bool  saveWeightImage;
    std::string  outputWeightFileName;    // []
    bool  saveBestFitParams;
    std::string  outputParameterFileName;
  
    bool printFitStatisticOnly;
    bool  noParamLimits;
    bool  ftolSet;
    double  ftol;
    int  solver;
    std::string  nloptSolverName;
  
    double  magZeroPoint;
  
    bool  printImages;
  
    bool  doBootstrap;
    int  bootstrapIterations;
    bool  saveBootstrap;
    std::string  outputBootstrapFileName;
  
    unsigned long  rngSeed;
  
};


#endif  // _OPTIONS_MAKEIMAGE_H_
