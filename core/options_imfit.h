/** @file
    \brief Subclass of OptionsBase class for imfit-specific program options.
 *
 */
#ifndef _OPTIONS_IMFIT_H_
#define _OPTIONS_IMFIT_H_

#include <string>
#include <vector>

#include "definitions.h"
#include "options_base.h"


//! Derived class for holding various options for imfit (set by command-line flags & options)
class ImfitOptions : public OptionsBase
{

  public:
    // Constructor:
    ImfitOptions( )
    {
      configFileName = DEFAULT_IMFIT_CONFIG_FILE;
  
      noImage = true;
      imageFileName = "";
  
      saveModel = false;
      outputModelFileName = "";
      saveResidualImage = false;
      outputResidualFileName = "";
      saveWeightImage = false;
      outputWeightFileName = "";
      saveBestFitParams = true;
      outputParameterFileName = DEFAULT_OUTPUT_PARAMETER_FILE;

      useModelForErrors = false;
      useCashStatistic = false;
      usePoissonMLR = false;

      noParamLimits = true;

      ftolSet = false;
      ftol = DEFAULT_FTOL;
      nloptSolverName = "NM";   // default value = Nelder-Mead Simplex

      magZeroPoint = NO_MAGNITUDES;
  
      printImages = false;

      doBootstrap = false;
      bootstrapIterations = 0;
      saveBootstrap = false;
      outputBootstrapFileName = "";
    };

    // Extra data members (in addition to those in options_base.h):  
    bool  noImage;
    std::string  imageFileName;
  
    std::string  psfOversampleRegion;
  
    bool  saveModel;
    std::string  outputModelFileName;
    bool  saveResidualImage;
    std::string  outputResidualFileName;
    bool  saveWeightImage;
    std::string  outputWeightFileName;
    bool  saveBestFitParams;
    std::string  outputParameterFileName;
  
    bool  noParamLimits;
    bool  ftolSet;
    double  ftol;
    std::string  nloptSolverName;
  
    double  magZeroPoint;
  
    bool  printImages;
  
    bool  doBootstrap;
    int  bootstrapIterations;
    bool  saveBootstrap;
    std::string  outputBootstrapFileName;
    
};


#endif  // _OPTIONS_IMFIT_H_
