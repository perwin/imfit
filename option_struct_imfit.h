/** @file
    \brief Struct containing options info for imfit, along with function for
           setting up default values within the struct.
 *
 */
#ifndef _IMFIT_OPTION_STRUCT_H_
#define _IMFIT_OPTION_STRUCT_H_

#include <string>

#include "definitions.h"


//! struct for holding various imfit options (set by command-line flags & options)
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
  std::string  psfOversampleRegion;     // []
  
  bool  noiseImagePresent;
  std::string  noiseFileName;   // []
  int  errorType;
  
  std::string  maskFileName;   //  []
  bool  maskImagePresent;
  int  maskFormat;
  
  bool  subsamplingFlag;
  
  bool  saveModel;
  std::string  outputModelFileName;   // []
  bool  saveResidualImage;
  std::string  outputResidualFileName;   // []
  bool  saveWeightImage;
  std::string  outputWeightFileName;    // []
  bool  saveBestFitParams;
  std::string  outputParameterFileName;
  
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
  
  int  maxThreads;
  bool  maxThreadsSet;
  
  int  verbose;
} imfitCommandOptions;


void SetDefaultImfitOptions( imfitCommandOptions *theOptions )
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
  theOptions->psfOversampleRegion = "";
  
  theOptions->noiseImagePresent = false;
  theOptions->noiseFileName = "";
  theOptions->errorType = WEIGHTS_ARE_SIGMAS;
  
  theOptions->maskImagePresent = false;
  theOptions->maskFileName = "";
  theOptions->maskFormat = MASK_ZERO_IS_GOOD;
  
  theOptions->subsamplingFlag = true;
  
  theOptions->saveModel = false;
  theOptions->outputModelFileName = "";
  theOptions->saveResidualImage = false;
  theOptions->outputResidualFileName = "";
  theOptions->saveWeightImage = false;
  theOptions->outputWeightFileName = "";
  theOptions->saveBestFitParams = true;
  theOptions->outputParameterFileName = DEFAULT_OUTPUT_PARAMETER_FILE;

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

  theOptions->printFitStatisticOnly = false;
  theOptions->noParamLimits = true;
  theOptions->ftolSet = false;
  theOptions->ftol = DEFAULT_FTOL;
  theOptions->solver = MPFIT_SOLVER;
  theOptions->nloptSolverName = "NM";   // default value = Nelder-Mead Simplex

  theOptions->magZeroPoint = NO_MAGNITUDES;
  
  theOptions->printImages = false;

  theOptions->doBootstrap = false;
  theOptions->bootstrapIterations = 0;
  theOptions->saveBootstrap = false;
  theOptions->outputBootstrapFileName = "";

  theOptions->maxThreads = 0;
  theOptions->maxThreadsSet = false;

  theOptions->verbose = 1;
}



#endif  // _IMFIT_OPTION_STRUCT_H_
