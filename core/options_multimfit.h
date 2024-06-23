/** @file
    \brief Subclass of OptionsBase class for multimfit-specific program options.
 *
 */
#ifndef _OPTIONS_MULTIMFIT_H_
#define _OPTIONS_MULTIMFIT_H_

#include <string>
#include <vector>

using namespace std;

#include "definitions.h"
#include "options_base.h"


//! Derived class for holding various options for multimfit (set by command-line flags & options)
class MultimfitOptions : public OptionsBase
{

  public:
    // Constructor:
    MultimfitOptions( )
    {
      configFileName = DEFAULT_IMFIT_CONFIG_FILE;
    
      saveModel = false;
      outputModelFileName = "";
      saveResidualImage = false;
      outputResidualFileName = "";
      saveWeightImage = false;
      outputWeightFileName = "";
      saveBestFitParams = true;
      outputParameterFileName = DEFAULT_OUTPUT_PARAMETER_FILE;
      outputParameterRootName = DEFAULT_BESTFIT_ROOTNAME;
      imageInfoFile_exists = false;
      imageInfoFile = "";

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
      
      loggingOn = false;
    };

    // Extra data members (in addition to those in options_base.h):  
  
    string  psfOversampleRegion;
  
    bool  saveModel;
    string  outputModelFileName;
    bool  saveResidualImage;
    string  outputResidualFileName;
    bool  saveWeightImage;
    string  outputWeightFileName;
    bool  saveBestFitParams;
    string  outputParameterFileName;
    string  outputParameterRootName;
  
    bool  noParamLimits;
    bool  ftolSet;
    double  ftol;
    string  nloptSolverName;
  
    double  magZeroPoint;
  
    bool  printImages;
  
    bool  doBootstrap;
    int  bootstrapIterations;
    bool  saveBootstrap;
    string  outputBootstrapFileName;
    
    bool  loggingOn;
    
    std::string  imageInfoFile;
    bool  imageInfoFile_exists;
};


#endif  // _OPTIONS_MULTIMFIT_H_
